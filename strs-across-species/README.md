This snakemake workflow (in progress) does the following:

1. Index fasta files
2. Split fastas by chromosomes
3. Call script to run Tandem Repeat Finder
4. Filter Tandem Repeats founded
5. Perform Statistical Analysis on STRs identified in each species
6. Visualize Statistical Analysis


Currently, it assumes genome fasta files are in:

```
$GENOMESDIR/$species/$species.fa
```

e.g. `/storage/mgymrek/TReeofLife/genomes/bosTau7/bosTau7.fa`. Note, the directory name is the same as the short species name to make automated changing of paths easier. This could be replaced in the snakefile with a path to a different directory containing genomes.

The output filtered STRs are written to `$TRFDIR`. The path to this can be changed in the snakefile.

The Statistical Analysis as well as figure generated are stored in `$YJDIR`. The path to this can also be changed in the snakefile.

To run, just type `snakemake`

To add more species, update the list of `SPECIES` in the snakefile, assuming their genome fasta files are available in the relevant location.
