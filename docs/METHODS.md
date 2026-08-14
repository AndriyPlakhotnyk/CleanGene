# CleanGene methods

CleanGene groups isolates by `group_id`. Paired reads can be screened by Kraken2 before assembly; when an expected organism is supplied, the isolate is excluded when the summed species-rank percentage outside that expected species is strictly greater than 5%, matching the ImaGene default rule.

Assemblies supplied in the manifest are reused. Otherwise Shovill is run with its own defaults. All retained assemblies are annotated by Prokka and submitted together to Panaroo using `--clean-mode strict`.

Panaroo calls are normalized to binary cluster x isolate calls. The standalone CleanGene default validates all clusters. Selected cluster representatives are recovered from `pan_genome_reference.fa`, falling back to a sequence-bearing member in `gene_data.csv` when required.

For every retained isolate, paired reads are aligned jointly to the selected cluster representatives with BWA-MEM. Secondary, supplementary, duplicate and low-MAPQ alignments are removed. SAMtools coverage measures breadth/depth, BCFtools constructs a haploid consensus with low-depth positions masked, and minimap2 aligns the consensus back to its cluster reference. Presence is confirmed only when breadth >=0.90, mean depth >=5, and consensus identity >=0.95. Reads below MAPQ 20 are excluded and BCFtools mpileup uses base quality >=30.

The final binary matrix changes selected Panaroo calls according to this read-supported evidence. A separate long evidence table retains initial call, final call, validation state, breadth, depth, identity and mapped-read count for every cluster/isolate combination.
