#!/bin/bash

GENOMESDIR=$1
SPECIES=$2

mkdir -p ${GENOMESDIR}/${SPECIES}/bychrom

# Get chromosomes
rm -f ${GENOMESDIR}/${SPECIES}/bychrom/${SPECIES}.chroms

for line in $(cat ${GENOMESDIR}/${SPECIES}/${SPECIES}.fa.fai | sed 's/\t/,/g')
do
    echo $chrom
    chrom=$(echo $line | cut -f 1 -d',')
    samtools faidx ${GENOMESDIR}/${SPECIES}/${SPECIES}.fa ${chrom} > ${GENOMESDIR}/${SPECIES}/bychrom/${chrom}.fa
    samtools faidx ${GENOMESDIR}/${SPECIES}/bychrom/${chrom}.fa
    echo $chrom >> ${GENOMESDIR}/${SPECIES}/bychrom/${SPECIES}.chroms
done
