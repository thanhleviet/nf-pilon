A simple nextflow script to run iteratively pilon for polishing long read contigs with Illumina reads

## Requirements

- conda and mamba or docker for managing bioinformatics software
- This pipeline requires nextflow edge version to use mamba.
  Installing the edge version:
  ```bash
  export NXF_EDGE=1
  nextflow self-update
  ```
## Usage

Prepare a sample_sheet.csv like below:

```
sample_id,sr1,sr2,contigs
sample_a,/path/to/short_read_1_sample_a,/path/to/short_read_2_sample_a,/path/to/contigs_sample_a
sample_b,/path/to/short_read_1_sample_b,/path/to/short_read_2_sample_b,/path/to/contigs_sample_b
```

### Run
```
nextflow run thanhleviet/nf-pilon \
--fofn /path/to/sample_sheet.csv \
--iterations 10 \
-profile conda
```