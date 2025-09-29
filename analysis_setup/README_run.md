# Quick start

```bash
cd analysis_setup

# (Recommended) Use mamba/conda env with dorado + samtools available in $PATH
# mamba activate nanopore

# Dry-run
snakemake -n

# Real run (8 threads rule-wise; adjust with --cores)
snakemake --cores 8
