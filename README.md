# mIHoP
Accurate detection of pathogens from low biomass samples

## SRA finder (new)

A small helper script to search NCBI SRA for human whole-genome runs (WGS) meeting basic criteria (read length, platform, number of reads). It uses NCBI EFetch `rettype=runinfo` to fetch run metadata and filters client-side.

Basic example:

```bash
python scripts/sra_finder.py \
  --min-read-length 100 \
  --min-reads 1000000 \
  --platform Illumina \
  --max-results 2000 \
  --output sra_runs.csv
```

Notes:
- Set `NCBI_EMAIL` or pass `--email` to identify yourself to NCBI.
- You can pass multiple platforms as a comma-separated list (e.g., `"Illumina,Oxford Nanopore"`).
- The script uses `avg_read_len = bases/spots` when available to estimate average read length.
