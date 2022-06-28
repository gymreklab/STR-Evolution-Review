#!/bin/bash

GENOMESDIR=$1
SPECIES=$2
OUTFILE=$3

echo "Running TRF on $SPECIES"

# trf config
matchscore=2
mismatchscore=5
indelscore=17
maxperiod=20 # Largest repeat unit
pm=80
pi=10
minscore=24 # Require at least 12 bp perfect matching
maxlen=1000

chroms=$(cat ${GENOMESDIR}/${SPECIES}/bychrom/${SPECIES}.chroms)

outdir=$(dirname $OUTFILE)
tmpdir=$(mktemp -d -t trf-XXXXXXXXXX --tmpdir=$outdir)

echo "writing to $tmpdir"

cd $tmpdir
# run trf one chromosome 
for chrom in $chroms
do
    trf ${GENOMESDIR}/${SPECIES}/bychrom/${chrom}.fa \
	${matchscore} ${mismatchscore} ${indelscore} \
	${pm} ${pi} \
	${minscore} \
	${maxperiod} \
	-d -h -l 1
done

ls $tmpdir # todo remove this

# concatenate results
for chrom in $chroms
do
    cat ${tmpdir}/${chrom}.fa.${matchscore}.${mismatchscore}.${indelscore}.${pm}.${pi}.${minscore}.${maxperiod}.dat | \
	awk -F' ' '(NF==15)' | \
	awk -v "chrom=$chrom" -F' ' '{print chrom,$1,$2,$3,$14,$15}' OFS='\t' | \
	awk -v "maxlen=${maxlen}" '(($3-$2) <= maxlen)'
done > ${OUTFILE}

rm -rf $tmpdir

exit 0
