#!/bin/bash

# usage: set a working directory on the commmand line - inside the working directory should be a directory called 'fastq' containing raw fastq reads with 'name'.raw.fastq.gz - 'name' is also set on the command line - also set trim quality score on command line

##### Set variables #####
time(

# command-line variables
WD=$1         # Working directory with adapter clipped fastq
READ1=$2	  # A comma separated list of read 1 files
READ2=$3      # A comma separated list of read 2 files
MINCONTIG=$4  # Minimum contig length 
PRESET=$5     # meta, meta-sensitive, meta-large, bulk

# software versions
MEGAHIT=1.0.2

# Assemble libraries
mkdir -p ${WD}/megahit_output
OUT=${WD}/megahit_output

MHITEXEC=~/marmic2017/software/megahit/megahit-1.0.2/megahit

# Run MEGAHIT
${MHITEXEC} -1 ${READ1} -2 ${READ2} -o ${OUT} --presets ${PRESET} --min-contig-len ${MINCONTIG}

# Collect stats
~/marmic2017/software/bbmap/stats.sh in=${OUT}/final.contigs.fa > ${OUT}/megahit.stats





