#!/usr/bin/env bash
set -ex

#need to have Mitty installed and running

FASTA=$1
SAMPLEVCF=$2
SAMPLENAME=$3
REGION_BED=Ch22.bed
FILTVCF=NA19471-filt.vcf.gz
READMODEL=1kg-pcr-free.pkl
COVERAGE=$4
READ_GEN_SEED=7
FASTQ_PREFIX=$5
READ_CORRUPT_SEED=8

mitty -v4 filter-variants \
  ${SAMPLEVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  - \
  2> vcf-filter.log | bgzip -c > ${FILTVCF}

tabix -p vcf ${FILTVCF}

mitty -v4 generate-reads \
  ${FASTA} \
  ${FILTVCF} \
  ${SAMPLENAME} \
  ${REGION_BED} \
  ${READMODEL} \
  ${COVERAGE} \
  ${READ_GEN_SEED} \
  >(gzip > ${FASTQ_PREFIX}.fq.gz) \
   ${FASTQ_PREFIX}-lq.txt \
   --unpair \
   --threads 2

mitty -v4 corrupt-reads \
  ${READMODEL} \
  ${FASTQ_PREFIX}.fq.gz >(gzip > ${FASTQ_PREFIX}-corrupt.fq.gz) \
  ${FASTQ_PREFIX}-lq.txt \
  ${FASTQ_PREFIX}-corrupt-lq.txt \
  ${READ_CORRUPT_SEED} \
  --threads 2
