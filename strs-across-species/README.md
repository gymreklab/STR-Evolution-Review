# Analysis of STR patterns across species

A filtered list of repeats identified in each species can be found here: https://drive.google.com/drive/folders/1bED0S_BY_BNWKyxaN3JFYaqoxzAfxxLB?usp=sharing

This folder contains a snakemake workflow used to generate these STR sets. The workflow takes in fasta files for each genome and performs the following steps:

1. Index fasta files
2. Split fastas by chromosomes
3. Run Tandem Repeat Finder (TRF)
4. Filter the TRF results
5. Perform Statistical Analysis on STRs identified in each species
6. Visualize Statistical Analysis

Currently, it assumes genome fasta files are in:

```
$GENOMESDIR/$species/$species.fa
```

Note, the directory name is the same as the short species name to make automated changing of paths easier. The output filtered STRs are written to `$TRFDIR`. The Statistical Analysis as well as figures generated are stored in `$YJDIR`. 

The paths for `$GENOMESDIR`, `$TRFDIR`, and `$YJDIR` can be changed in the snakefile.

To run, just type `snakemake`

To add more species, update the list of `SPECIES` in the snakefile, assuming their genome fasta files are available in the relevant location.
