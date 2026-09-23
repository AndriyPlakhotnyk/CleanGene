from __future__ import annotations

import copy
import importlib.util
import io
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from cleangene.cli import apply_cli_overrides, local, main
from cleangene.defaults import DEFAULTS
from cleangene.resistance import (align_sequences, assign_families, cluster_sequences, controller,
    definitions, extract_loci, features_from_gff, global_identity, merge, origins, parse_hits,
    pileup_counts, read_support, revcomp, root_dir, validate_config, variant_events)
from cleangene.util import atomic_json, load_json, read_tsv, write_tsv
from cleangene.workers import _index_done, _stage_limit

HAS_EDLIB = importlib.util.find_spec("edlib") is not None
HAS_PLOTS = importlib.util.find_spec("matplotlib") is not None


def fixture(isolate="a", strand="+", mutation=False):
    rng = random.Random(501)
    sequence = "".join(rng.choice("ACGT") for _ in range(400))
    if mutation:
        sequence = sequence[:200] + ("T" if sequence[200] != "T" else "C") + sequence[201:]
    features = []
    for number, (start, end) in enumerate(((10, 39), (80, 179), (220, 249))):
        features.append({"id": f"{isolate}_{number}", "contig": "c", "start": start, "end": end,
            "strand": strand, "gene": "", "product": "hypothetical protein", "sequence": sequence[start-1:end]})
    hit = {"symbol": "tet(M)", "contig": "c", "start": 80, "end": 179, "strand": strand, "method": "EXACTX"}
    result = extract_loci(isolate, "Efaecium", {"c": sequence}, features, [hit], {}, 10000)
    for locus in result:
        locus.update(evidence="assembly_only", read_supported_fraction="NA", mean_depth="NA", discordant_positions="NA", read_note="fixture")
    return sequence, features, hit, result


class ResistanceTests(unittest.TestCase):
    def test_requested_typo_alias_and_canonical_flag(self):
        for flag in ("-ressitanec-operon", "--resistance-operon", "-resistance-operon", "--ressitanec-operon"):
            with patch("cleangene.cli.run_command", return_value=0) as run:
                self.assertEqual(main(["run", "--manifest", "samples.tsv", flag, "--skip-downsampling", "--ignore-checkm2"]), 0)
                cfg = apply_cli_overrides(DEFAULTS, run.call_args.args[0])
                self.assertEqual(cfg["RESISTANCE_OPERON"], "true")
                self.assertEqual(cfg["SKIP_DOWNSAMPLING"], "true")
                self.assertEqual(cfg["CHECKM2_MODE"], "off")

    def test_local_stage_order_and_no_checkm2(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp)
            cfg = {**DEFAULTS, "CHECKM2_MODE": "off", "TAXONOMY_MODE": "off", "RESISTANCE_OPERON": "true"}
            atomic_json(run / "provenance/resolved_config.json", cfg)
            rows = [["a", "g"], ["b", "g"]]
            for path in ("provenance/manifest.tsv", "state/isolate_tasks.tsv"):
                write_tsv(run / path, ["isolate_id", "group_id"], rows)
            write_tsv(run / "state/group_tasks.tsv", ["group_id"], [["g"]])
            with patch("cleangene.cli.dispatch") as dispatch, patch("cleangene.resistance.prepare"):
                local(run)
            calls = [c.args[0] for c in dispatch.call_args_list]
            self.assertNotIn("checkm2_db_setup", calls)
            self.assertEqual(calls[-4:], ["resistance_scan", "resistance_scan", "resistance_merge", "summary"])
            self.assertLess(calls.index("plot"), calls.index("resistance_scan"))

    def test_slurm_thousand_isolate_indices_and_resources(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp)
            write_tsv(run / "state/isolate_tasks.tsv", ["isolate_id", "group_id"], [[f"i{i}", "g"] for i in range(999)])
            with patch("cleangene.resistance.prepare"), patch("cleangene.workers._run_index_stage") as scan, patch("cleangene.workers._run_single_job") as finish:
                controller(run, DEFAULTS)
            self.assertEqual(scan.call_args.args[2], "resistance_scan")
            self.assertEqual(scan.call_args.args[3], list(range(999)))
            self.assertEqual(finish.call_args.args[2], "resistance_merge")
            atomic_json(run / "state/resistance_scan/i998.done.json", {"status": "complete"})
            self.assertTrue(_index_done(run, "resistance_scan", 998))
            self.assertEqual(_stage_limit(DEFAULTS, "resistance_scan"), 100)

    def test_extract_exactly_one_flank_each_side(self):
        sequence, _, _, loci = fixture()
        locus = loci[0]
        self.assertEqual((locus["start"], locus["end"]), (10, 249))
        self.assertEqual(locus["sequence"], sequence[9:249])
        self.assertEqual([f["role"] for f in locus["features"]], ["upstream", "core", "downstream"])
        self.assertEqual(locus["context_status"], "complete")

    def test_reverse_orientation_flanks(self):
        sequence, _, _, loci = fixture(strand="-")
        locus = loci[0]
        self.assertEqual(locus["sequence"], revcomp(sequence[9:249]))
        self.assertEqual([f["id"] for f in locus["features"]], ["a_2", "a_1", "a_0"])
        self.assertEqual([f["role"] for f in locus["features"]], ["upstream", "core", "downstream"])

    def test_nearby_duplicate_determinants_are_separate_loci(self):
        sequence, features, hit, _ = fixture()
        other = {**hit, "start": 220, "end": 249}
        loci = extract_loci("a", "g", {"c": sequence}, features, [hit, other], {}, 10000)
        self.assertEqual(len(loci), 2)
        self.assertEqual(len({l["locus_id"] for l in loci}), 2)

    def test_fragmented_contigs_not_stitched_and_missing_flank_unknown(self):
        sequence, features, hit, _ = fixture()
        first = {**hit, "symbol": "vanA"}
        second = {**hit, "contig": "c2", "symbol": "vanX-A"}
        loci = extract_loci("a", "g", {"c": sequence, "c2": sequence}, [features[1]], [first, second], definitions(DEFAULTS), 10000)
        self.assertEqual(len(loci), 2)
        self.assertTrue(all(l["context_status"] == "incomplete" for l in loci))
        self.assertTrue(all(l["upstream_missing"] for l in loci))
        self.assertTrue(all(l["unobserved_expected_genes"] for l in loci))

    def test_defined_multigene_span_includes_internal_non_amr_cds(self):
        sequence, features, hit, _ = fixture()
        features[0]["gene"] = "vanR-A"
        features[2]["gene"] = "vanX-A"
        loci = extract_loci("a", "g", {"c": sequence}, features, [{**hit, "symbol": "vanA"}], definitions(DEFAULTS), 10000)
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]["observed_genes"], ["vanA", "vanR-A", "vanX-A"])
        self.assertEqual(len(loci[0]["features"]), 3)

    def test_amrfinder_new_and_old_headers_filter_non_amr(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "amr.tsv"
            for field in ("Element symbol", "Gene symbol"):
                write_tsv(path, [field, "Contig id", "Start", "Stop", "Strand", "Type", "Method"],
                    [["tet(M)", "c", 179, 80, "-", "AMR", "PARTIALX"], ["stress", "c", 200, 300, "+", "STRESS", "EXACTP"]])
                hits = parse_hits(path)
                self.assertEqual(len(hits), 1)
                self.assertEqual((hits[0]["start"], hits[0]["end"]), (80, 179))
            path.write_text("bad\theader\n")
            with self.assertRaises(ValueError):
                parse_hits(path)

    @unittest.skipUnless(HAS_EDLIB, "edlib required")
    def test_point_mutations_bin_by_gene_not_amino_acid_allele(self):
        sequence, features, hit, _ = fixture()
        loci = []
        for iso, symbol in (("a", "gyrA_S83L"), ("b", "gyrA_S83Y")):
            loci += extract_loci(iso, "g", {"c": sequence}, features,
                [{**hit, "symbol": symbol, "subtype": "POINT", "method": "POINTP"}], {}, 10000)
        assign_families(loci, {})
        self.assertTrue(all(l["operon_type"] == "gyrA" for l in loci))
        self.assertEqual(loci[0]["category_id"], loci[1]["category_id"])
        self.assertNotEqual(loci[0]["reported_amr_symbols"], loci[1]["reported_amr_symbols"])

    def test_gff_inconsistent_coordinates_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "a.gff"
            path.write_text("c\ttest\tCDS\t1\t101\t.\t+\t0\tID=x\n")
            with self.assertRaises(ValueError):
                features_from_gff(path, {"c": "A"*100})

    def test_pileup_quality_indels_deletions_and_read_boundaries(self):
        self.assertEqual(pileup_counts("^F.,$Aa+2tt.,-3ACG*"), (4, 7, 2))
        self.assertEqual(pileup_counts(".,><N#"), (2, 4, 0))
        with self.assertRaises(ValueError):
            pileup_counts(".+A")

    def test_snp_insertion_deletion_and_ambiguous_events(self):
        events = variant_events("ACG--TACN", "ATGCGT-CG")
        self.assertEqual([e["event"] for e in events], ["SNP", "insertion", "deletion", "uncertain"])
        self.assertEqual(events[1]["reference_position"], 3)
        self.assertEqual(events[1]["alt"], "CG")
        self.assertEqual(events[2]["ref"], "A")
        self.assertEqual(variant_events("-AC", "TAC")[0]["reference_position"], 0)
        self.assertEqual(variant_events("A"*51, "-"*51)[0]["scale"], "large")
        with self.assertRaises(ValueError):
            variant_events("AC", "A")

    @unittest.skipUnless(HAS_EDLIB, "edlib required")
    def test_identity_95_boundary_and_global_length_penalty(self):
        a = "A" * 100
        self.assertAlmostEqual(global_identity(a, "C"*5 + "A"*95, .95), .95)
        self.assertLess(global_identity(a, "C"*6 + "A"*94, .95), .95)
        self.assertEqual(global_identity("A"*50, a, .95), 0)
        self.assertLess(global_identity("N"*100, "N"*100), .95)
        sequences = {"a": a, "b": "C"*5+"A"*95, "c": "C"*10+"A"*90}
        groups = cluster_sequences(sequences, {"a": 3, "b": 2, "c": 1}, .95)
        self.assertEqual(groups, {"a": "a", "b": "a", "c": "c"})

    @unittest.skipUnless(HAS_EDLIB, "edlib required")
    def test_gene_content_independent_of_order_and_hypothetical_label(self):
        _, _, _, loci = fixture()
        second = copy.deepcopy(loci[0])
        second["isolate_id"] = "b"
        second["features"].reverse()
        assign_families([loci[0], second], {})
        self.assertEqual(loci[0]["category_id"], second["category_id"])
        self.assertNotEqual(loci[0]["gene_order"], second["gene_order"])
        self.assertNotEqual(loci[0]["features"][0]["family"], loci[0]["features"][2]["family"])

    def test_no_reads_is_explicit_assembly_only(self):
        _, _, _, loci = fixture()
        read_support(loci, Path("unused"), {}, Path("unused"), DEFAULTS)
        self.assertEqual(loci[0]["evidence"], "assembly_only")
        self.assertEqual(loci[0]["read_supported_fraction"], "NA")

    def test_read_discordance_is_not_hidden_by_breadth(self):
        from cleangene.fasta import write_fasta
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            assembly = out / "input.fa"
            write_fasta(assembly, [("c", "ACGT"), ("competing_copy", "ACGT")])
            _, _, _, loci = fixture()
            loci[0].update(start=1, end=4)
            process = SimpleNamespace(stdout=io.StringIO(
                "c\t1\tA\t5\t.....\tIIIII\n"
                "c\t2\tC\t5\t.....\tIIIII\n"
                "c\t3\tG\t5\tAAAAA\tIIIII\n"
                "c\t4\tT\t5\t.....\tIIIII\n"), wait=lambda: 0)
            cfg = {**DEFAULTS, "RESISTANCE_MIN_BREADTH": "0.75"}
            with patch("cleangene.evidence.map_reads") as mapping, patch("cleangene.resistance.run"), patch("cleangene.resistance.subprocess.Popen", return_value=process):
                read_support(loci, assembly, {"R1": "reads.fq"}, out, cfg)
            self.assertEqual(loci[0]["read_supported_fraction"], .75)
            self.assertEqual(loci[0]["discordant_positions"], 1)
            self.assertEqual(loci[0]["evidence"], "assembly_read_warning")
            self.assertTrue((out / "base_support.tsv.gz").is_file())
            self.assertEqual(mapping.call_args.args[0].name, "mapping_reference.fasta")

    def test_amr_protein_identity_does_not_label_overlapping_neighbor(self):
        sequence, features, hit, _ = fixture()
        features.append({**features[0], "id": "overlap", "start": 170, "end": 195})
        hit["protein"] = features[1]["id"]
        loci = extract_loci("a", "g", {"c": sequence}, features, [hit], {}, 10000)
        labelled = [f for f in loci[0]["features"] if f["amr_symbols"]]
        self.assertEqual([f["id"] for f in labelled], [features[1]["id"]])

    def test_missing_mafft_fails_without_fake_alignment(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("cleangene.resistance.run", side_effect=FileNotFoundError("mafft")):
                with self.assertRaises(FileNotFoundError):
                    align_sequences({"a": "ACGT", "b": "ACGGT"}, Path(temp)/"out.fasta", 1)

    def test_origin_evidence_and_invalid_settings(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"origins.tsv"
            write_tsv(path, ["isolate_id", "contig", "origin", "evidence"], [["a", "c", "plasmid", "closed curated replicon"]])
            self.assertEqual(origins(path)["a", "c"]["origin"], "plasmid")
            write_tsv(path, ["isolate_id", "contig", "origin", "evidence"], [["a", "c", "plasmid", ""]])
            with self.assertRaises(ValueError): origins(path)
        for changes in ({"RESISTANCE_IDENTITY": "95"}, {"RESISTANCE_TOP": "0"}, {"RESISTANCE_MIN_DEPTH": "-1"}):
            with self.assertRaises(ValueError): validate_config({**DEFAULTS, **changes})

    @unittest.skipUnless(HAS_EDLIB and HAS_PLOTS and shutil.which("mafft"), "alignment/plot environment required")
    def test_real_mafft_reporting_distinct_isolates_bins_and_zero_hits(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp)
            atomic_json(run / "provenance/resolved_config.json", {**DEFAULTS, "RESISTANCE_MERGE_CPUS": "1"})
            write_tsv(run / "provenance/manifest.tsv", ["isolate_id", "group_id"], [[i, "Efaecium"] for i in "abcd"])
            write_tsv(run / "state/isolate_tasks.tsv", ["isolate_id", "group_id"], [[i, "Efaecium"] for i in "abcd"])
            for iso in "abcd":
                loci = fixture(isolate=iso, mutation=iso == "b")[3] if iso in "ab" else []
                if iso == "a":
                    loci += [{**copy.deepcopy(loci[0]), "locus_id": "duplicate_copy"}]
                atomic_json(root_dir(run) / "isolates" / iso / "loci.json",
                    {"isolate_id": iso, "group_id": "Efaecium", "status": "excluded" if iso == "d" else "evaluated", "loci": loci})
                atomic_json(run / "state/resistance_scan" / f"{iso}.done.json", {"status": "complete"})
            merge(run)
            tables = root_dir(run) / "tables"
            prevalence = read_tsv(tables / "operon_prevalence.tsv")[0]
            self.assertEqual((prevalence["n_present"], prevalence["n_evaluated"], prevalence["n_copies"]), ("2", "3", "3"))
            self.assertEqual((prevalence["n_categories"], prevalence["n_bins"], prevalence["n_exact_variants"]), ("1", "1", "2"))
            timestamp = (tables / "loci.tsv").stat().st_mtime_ns
            merge(run)
            self.assertEqual(timestamp, (tables / "loci.tsv").stat().st_mtime_ns)
            (root_dir(run)/"figures/operon_prevalence.png").unlink()
            merge(run)
            self.assertTrue((root_dir(run)/"figures/operon_prevalence.png").is_file())
            self.assertNotEqual(timestamp, (tables / "loci.tsv").stat().st_mtime_ns)
            events = read_tsv(tables / "variant_events.tsv")
            self.assertEqual(len(events), 1)
            self.assertEqual(events[0]["event"], "SNP")
            self.assertTrue(list((root_dir(run)/"figures").rglob("alignment_sites.pdf")))
            self.assertTrue((run/"state/resistance_merge.done.json").is_file())
            matrix = read_tsv(tables/"isolate_presence_absence.tsv")
            self.assertEqual(list(matrix[-1].values())[-1], "NA")

    @unittest.skipUnless(HAS_PLOTS, "matplotlib required")
    def test_zero_hit_cohort_has_valid_empty_reports(self):
        with tempfile.TemporaryDirectory() as temp:
            run = Path(temp)
            atomic_json(run/"provenance/resolved_config.json", DEFAULTS)
            write_tsv(run/"provenance/manifest.tsv", ["isolate_id", "group_id"], [["a", "g"]])
            write_tsv(run/"state/isolate_tasks.tsv", ["isolate_id", "group_id"], [["a", "g"]])
            atomic_json(root_dir(run)/"isolates/a/loci.json", {"isolate_id": "a", "group_id": "g", "status": "evaluated", "loci": []})
            atomic_json(run/"state/resistance_scan/a.done.json", {"status": "complete"})
            merge(run)
            self.assertEqual(read_tsv(root_dir(run)/"tables/loci.tsv"), [])
            self.assertTrue((root_dir(run)/"figures/operon_prevalence.png").is_file())


if __name__ == "__main__":
    unittest.main()
