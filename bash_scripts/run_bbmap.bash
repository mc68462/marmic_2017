# Read mapping
# For downstream use in metabat and anvio

WD=$1         # full path of parent directory 
NAME=$2       # sample name / prefix to individual NAME fastq files
READDIR=$3    # full path of directory with quality trimmed reads to be mapped to assembly - fastq
ASSEMDIR=$4   # full path of directory with assembly contigs -fasta
THREADS=$5    # number of threads for samtools-rocksort

# software versions
BBMAP=35.14
SAMTOOLS=1.2
ROCKSORT=0.2
PICCARD=1.119

# Make directories
mkdir -p ${WD}/bbmap_output
MAPDIR=${WD}/bbmap_output

cd ${MAPDIR}

# run bbmap (global alignments)
BBEXEC=~/marmic2017/software/bbmap/bbmap-35.14/bbmap.sh

${BBEXEC} ref=${ASSEMDIR}/final.contigs.fa
${BBEXEC} in1=${READDIR}/${NAME}.R1.clean.fastq.gz in2=${READDIR}/${NAME}.R2.clean.fastq.gz out=${MAPDIR}/${NAME}.sam fast=t idfilter=0.97 idtag=t covstats=${MAPDIR}/${NAME}.cov.stats

SAMEXEC=~/marmic2017/software/samtools/samtools-1.2/bin/samtools
# fast hacker version of samtools called samtools-rocksort (http://devblog.dnanexus.com/tag/samtools-rocksort/)
ROCEXEC=~/marmic2017/software/samtools-rocksort/samtools-rocksort-0.2/bin/samtools

${SAMEXEC} faidx ${ASSEMDIR}/final.contigs.fa 

# convert paired sam to sorted bam and sort using rocksort

${SAMEXEC} view -@ ${THREADS} -q 10 -F 4 -bt ${ASSEMDIR}/final.contigs.fa.fai ${MAPDIR}/${NAME}.sam | ${ROCEXEC} rocksort -@ ${THREADS} -m 3 - ${NAME}.sorted

rmdir ${MAPDIR}/${NAME}.sorted.rocksort
rm ${MAPDIR}/${NAME}.sam 

# mark duplicates with piccard tools

MARKDUPEXEC=~/marmic2017/picard-tools-1.119/MarkDuplicates.jar

JAVA_OPT="-Xms2g -Xmx32g -XX:ParallelGCThreads=4 -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled"
java ${JAVA_OPT} \
     -jar ${MARKDUPEXEC} \
     INPUT=${MAPDIR}/${NAME}.sorted.bam \
     OUTPUT=${MAPDIR}/${NAME}.sorted.markdup.bam \
     METRICS_FILE=${MAPDIR}/${NAME}.sorted.markdup.metrics \
     AS=TRUE \
     VALIDATION_STRINGENCY=LENIENT \
     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
     REMOVE_DUPLICATES=TRUE

# resort, index and collect stats on resulted markdup bam file

${ROCEXEC} rocksort -@ ${THREADS} -m 3 ${MAPDIR}/${NAME}.sorted.markdup.bam ${MAPDIR}/${NAME}.sorted.markdup.sorted
${SAMEXEC} index ${MAPDIR}/${NAME}.sorted.markdup.sorted.bam
${SAMEXEC} flagstat ${MAPDIR}/${NAME}.sorted.markdup.sorted.bam > ${MAPDIR}/${NAME}.sorted.markdup.sorted.flagstat

rmdir ${MAPDIR}/${NAME}.sorted.markdup.sorted.rocksort
