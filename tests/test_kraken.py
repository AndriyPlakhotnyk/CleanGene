import multiprocessing
import tempfile
import time
import unittest
from pathlib import Path

from cleangene.kraken import Kraken2DbError, Kraken2DbNotReady, managed_kraken2_db_path, normalize_kraken2_database_size, resolve_kraken2_db, validate_kraken2_db
from cleangene.util import atomic_json, load_json, write_tsv
from cleangene.workers import ensure_kraken2_db, kraken_db_setup


def _write_valid_db(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (path / name).write_text(name + "\n")


def _concurrent_resolve_worker(root: str, calls_dir: str, results_dir: str, errors_dir: str) -> None:
    cfg = {"KRAKEN2_DATABASE_ROOT": root, "KRAKEN2_DATABASE_SIZE": "standard-8", "KRAKEN2_AUTO_DOWNLOAD": "true", "KRAKEN2_CLEAN_BUILD_FILES": "true"}

    def runner(command):
        Path(calls_dir, f"call-{time.time_ns()}-{Path(command[1]).name}").write_text("\t".join(command))
        time.sleep(0.2)
        _write_valid_db(Path(command[1]))

    try:
        result = resolve_kraken2_db(cfg, allow_download=True, runner=runner)
        Path(results_dir, f"result-{time.time_ns()}").write_text(str(result.path))
    except Exception as error:
        Path(errors_dir, f"error-{time.time_ns()}").write_text(str(error))


class Kraken2ResolutionTests(unittest.TestCase):
    def test_existing_shared_db_is_selected_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"
            cfg = {"TAXONOMY_MODE": "kraken2", "KRAKEN2_DATABASE_ROOT": str(Path(d) / "dbs"), "KRAKEN2_DATABASE_SIZE": "8", "KRAKEN2_AUTO_DOWNLOAD": "true"}
            db, _, _ = managed_kraken2_db_path(cfg)
            _write_valid_db(db)
            atomic_json(run / "provenance" / "resolved_config.json", cfg)
            write_tsv(run / "provenance" / "manifest.tsv", ["isolate_id", "group_id", "grouping_source"], [["i1", "g", "manifest_group_id"]])
            kraken_db_setup(run)
            resolved = load_json(run / "provenance" / "resolved_config.json")
            marker = load_json(run / "state" / "kraken_db_setup.done.json")
            self.assertEqual(resolved["KRAKEN2_DB"], str(db))
            self.assertEqual(resolved["KRAKEN2_DATABASE_SIZE"], "standard-8")
            self.assertEqual(marker["source"], "shared_existing")

    def test_missing_shared_db_downloads_once_and_records_path(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = {"KRAKEN2_DATABASE_ROOT": str(Path(d) / "dbs"), "KRAKEN2_DATABASE_SIZE": "standard-16", "KRAKEN2_AUTO_DOWNLOAD": "true", "KRAKEN2_CLEAN_BUILD_FILES": "true"}
            calls = []

            def runner(command):
                calls.append(command)
                _write_valid_db(Path(command[1]))

            result = resolve_kraken2_db(cfg, allow_download=True, runner=runner)
            self.assertEqual(len(calls), 1)
            self.assertEqual(result.source, "auto_download")
            self.assertEqual(result.path.name, "kraken2_standard-16")
            validate_kraken2_db(result.path)

    def test_partial_and_empty_databases_are_invalid(self):
        with tempfile.TemporaryDirectory() as d:
            partial = Path(d) / "partial"; partial.mkdir(); (partial / "hash.k2d").write_text("hash\n")
            with self.assertRaisesRegex(Kraken2DbError, "opts.k2d"):
                validate_kraken2_db(partial)
            empty = Path(d) / "empty"; _write_valid_db(empty); (empty / "opts.k2d").write_text("")
            with self.assertRaisesRegex(Kraken2DbError, "empty"):
                validate_kraken2_db(empty)

    def test_explicit_override_valid_and_invalid(self):
        with tempfile.TemporaryDirectory() as d:
            custom = Path(d) / "custom"; _write_valid_db(custom)
            cfg = {"KRAKEN2_DB": str(custom), "KRAKEN2_DATABASE_ROOT": str(Path(d) / "shared"), "KRAKEN2_DATABASE_SIZE": "standard-8"}
            result = resolve_kraken2_db(cfg, allow_download=True, runner=lambda command: (_ for _ in ()).throw(AssertionError("download should not run")))
            self.assertEqual(result.path, custom.resolve())
            self.assertEqual(result.source, "explicit")
            cfg["KRAKEN2_DB"] = str(Path(d) / "missing")
            with self.assertRaisesRegex(Kraken2DbError, "Explicit KRAKEN2_DB is invalid"):
                resolve_kraken2_db(cfg, allow_download=True, runner=lambda command: (_ for _ in ()).throw(AssertionError("download should not run")))

    def test_size_selection_and_aliases(self):
        aliases = {"standard-8": "standard-8", "standard_8": "standard-8", "8": "standard-8", "standard-16": "standard-16", "standard_16": "standard-16", "16": "standard-16", "standard": "standard", "full": "standard", "full-standard": "standard"}
        for alias, canonical in aliases.items():
            self.assertEqual(normalize_kraken2_database_size(alias), canonical)
        with tempfile.TemporaryDirectory() as d:
            for value, dirname in (("standard-8", "kraken2_standard-8"), ("standard-16", "kraken2_standard-16"), ("standard", "kraken2_standard")):
                path, size, _ = managed_kraken2_db_path({"KRAKEN2_DATABASE_ROOT": d, "KRAKEN2_DATABASE_SIZE": value})
                self.assertEqual(path.name, dirname)
                self.assertEqual(size, value)
        with self.assertRaisesRegex(Kraken2DbError, "must be one of"):
            normalize_kraken2_database_size("mini")

    def test_concurrent_missing_db_downloads_once(self):
        with tempfile.TemporaryDirectory() as d:
            calls = Path(d) / "calls"; results_dir = Path(d) / "results"; errors_dir = Path(d) / "errors"
            calls.mkdir(); results_dir.mkdir(); errors_dir.mkdir()
            processes = [multiprocessing.Process(target=_concurrent_resolve_worker, args=(str(Path(d) / "dbs"), str(calls), str(results_dir), str(errors_dir))) for _ in range(2)]
            for process in processes: process.start()
            for process in processes: process.join(10)
            for process in processes: self.assertEqual(process.exitcode, 0)
            self.assertEqual(list(errors_dir.iterdir()), [])
            self.assertEqual(len(list(calls.iterdir())), 1)
            results = [path.read_text() for path in results_dir.iterdir()]
            self.assertEqual(len(set(results)), 1)
            validate_kraken2_db(Path(results[0]))

    def test_auto_download_false_uses_existing_but_fails_missing(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = {"KRAKEN2_DATABASE_ROOT": str(Path(d) / "dbs"), "KRAKEN2_DATABASE_SIZE": "standard", "KRAKEN2_AUTO_DOWNLOAD": "false"}
            db, _, _ = managed_kraken2_db_path(cfg)
            _write_valid_db(db)
            self.assertEqual(resolve_kraken2_db(cfg).path, db)
            for child in db.iterdir(): child.unlink()
            db.rmdir()
            with self.assertRaisesRegex(Kraken2DbNotReady, "KRAKEN2_AUTO_DOWNLOAD=false"):
                resolve_kraken2_db(cfg)

    def test_kraken_not_required_does_not_validate_or_download(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"; run.mkdir()
            cfg = {"TAXONOMY_MODE": "off", "KRAKEN2_DB": str(Path(d) / "missing")}
            rows = [{"isolate_id": "i1", "group_id": "g", "grouping_source": "manifest_group_id"}]
            self.assertEqual(ensure_kraken2_db(run, cfg, rows), "")

    def test_resume_reuses_valid_recorded_db(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"; db = Path(d) / "recorded"; _write_valid_db(db)
            cfg = {"TAXONOMY_MODE": "kraken2", "KRAKEN2_DB": str(db), "KRAKEN2_DATABASE_ROOT": str(Path(d) / "shared")}
            rows = [{"isolate_id": "i1", "group_id": "g", "grouping_source": "manifest_group_id"}]
            self.assertEqual(ensure_kraken2_db(run, cfg, rows), str(db.resolve()))


if __name__ == "__main__":
    unittest.main()
