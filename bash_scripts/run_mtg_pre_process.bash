#!/bin/bash

# usage: set a working directory on the commmand line - inside the working directory should be a directory called 'fastq' containing raw fastq reads with 'name'.raw.fastq.gz - 'name' is also set on the command line - also set trim quality score on command line

##### Set variables #####

time(
# command-line variables
NAME=$1     # sample name (ie 20110404)
WD=$2		# working directory starting from root
RAWDIR=$3   # raw read directory (ie fastq_files)
TRIMQS=$4   # quality score threshold (ie 10, 20, 30)

# hard-coded variables
MINLEN=100

# software versions
FASTQCVER=0.11.2
BBDUKVER=35.14

FASTQ=${RAWDIR}/${NAME}.raw.fastq.gz

##### Split R1 and R2 #####

SPLITEXE=~/marmic2017/software/bbmap/bbmap-35.14/reformat.sh
R1=${RAWDIR}/${NAME}.R1.raw.fastq.gz
R2=${RAWDIR}/${NAME}.R2.raw.fastq.gz

echo "Splitting reads 1 and 2 for sample:" ${NAME}

${SPLITEXE} in=${FASTQ} out1=${R1} out2=${R2}


##### Trim adapters and quality filter using bbduk #####

mkdir -p ${WD}/clean_fastq
TRIMWD=${WD}/clean_fastq

TRIMEXE=~/marmic2017/software/bbmap/bbmap-35.14/bbduk.sh
ADAPTERS=~/marmic2017/software/bbmap/bbmap-35.14/resources/truseq.fa.gz
R1CLEAN=${TRIMWD}/${NAME}.R1.clean.fastq.gz
R2CLEAN=${TRIMWD}/${NAME}.R2.clean.fastq.gz

echo "Clipping and quality trimming reads for sample:" ${NAME}

${TRIMEXE} in=${FASTQ} out1=${R1CLEAN} out2=${R2CLEAN} ref=${ADAPTERS} ktrim=r k=28 mink=12 hdist=1 tbo=t tpe=t qtrim=rl trimq=${TRIMQS} minlength=${MINLEN} 

##### Run FASTQC #####

mkdir -p ${WD}/fastqc_pre-clean
mkdir -p ${WD}/fastqc_post-clean
mkdir -p ${WD}/qc_reports
FASTQCWDPRE=${WD}/fastqc_pre-clean
FASTQCWDPOST=${WD}/fastqc_post-clean
REP=${WD}/qc_reports

FASTQCEXE=~/marmic2017/software/fastqc/fastqc-0.11/fastqc

echo "Running fastqc for sample:" ${NAME}

${FASTQCEXE} ${R1} ${R2} -o ${FASTQCWDPRE}

${FASTQCEXE} ${R1CLEAN} ${R2CLEAN} -o ${FASTQCWDPOST}

illu1REP=${REP}/${NAME}.R1
illu2REP=${REP}/${NAME}.R2

touch ${illu1REP}
touch ${illu2REP}

cd ${FASTQCWDPRE}

unzip ${NAME}.R1.raw_fastqc.zip
unzip ${NAME}.R2.raw_fastqc.zip

cd ${FASTQCWDPOST}

unzip ${NAME}.R1.clean_fastqc.zip
unzip ${NAME}.R2.clean_fastqc.zip

cd ${WD}

echo "Analyzing fastqc & collecting sequence stats for sample:" ${NAME}

#Analyze fastqc output
R1FQCPREA=$(grep Adapter ${FASTQCWDPRE}/${NAME}.R1.raw_fastqc/summary.txt | cut -f 1)
R2FQCPREA=$(grep Adapter ${FASTQCWDPRE}/${NAME}.R2.raw_fastqc/summary.txt | cut -f 1)
R1FQCPRED=$(grep Duplication ${FASTQCWDPRE}/${NAME}.R1.raw_fastqc/summary.txt | cut -f 1)
R2FQCPRED=$(grep Duplication ${FASTQCWDPRE}/${NAME}.R2.raw_fastqc/summary.txt | cut -f 1)

R1FQCPOSTA=$(grep Adapter ${FASTQCWDPOST}/${NAME}.R1.clean_fastqc/summary.txt | cut -f 1)
R2FQCPOSTA=$(grep Adapter ${FASTQCWDPOST}/${NAME}.R2.clean_fastqc/summary.txt | cut -f 1)
R1FQCPOSTD=$(grep Duplication ${FASTQCWDPOST}/${NAME}.R1.clean_fastqc/summary.txt | cut -f 1)
R2FQCPOSTD=$(grep Duplication ${FASTQCWDPOST}/${NAME}.R2.clean_fastq/summary.txt | cut -f 1)

#Count sequence statistics
L1PRE=$(zcat ${RAWDIR}/${NAME}.R1.raw.fastq.gz | awk 'BEGIN{L=0;i=0; max=0; min=0;}{if(NR%4==2){if (i == 0){max = length($0); min = length($0)}else{if (length($0) >= max){max = length($0)}; if (length($0) <= min){min = length($0)} }; L=L+length($0); i++}}END{print i"\t"int(L/i)"\t"max"\t"min}')
L2PRE=$(zcat ${RAWDIR}/${NAME}.R2.raw.fastq.gz | awk 'BEGIN{L=0;i=0; max=0; min=0;}{if(NR%4==2){if (i == 0){max = length($0); min = length($0)}else{if (length($0) >= max){max = length($0)}; if (length($0) <= min){min = length($0)} }; L=L+length($0); i++}}END{print i"\t"int(L/i)"\t"max"\t"min}')

L1POST=$(zcat ${TRIMWD}/${NAME}.R1.clean.fastq.gz | awk 'BEGIN{L=0;i=0; max=0; min=0;}{if(NR%4==2){if (i == 0){max = length($0); min = length($0)}else{if (length($0) >= max){max = length($0)}; if (length($0) <= min){min = length($0)} }; L=L+length($0); i++}}END{print i"\t"int(L/i)"\t"max"\t"min}')
L2POST=$(zcat  ${TRIMWD}/${NAME}.R2.clean.fastq.gz | awk 'BEGIN{L=0;i=0; max=0; min=0;}{if(NR%4==2){if (i == 0){max = length($0); min = length($0)}else{if (length($0) >= max){max = length($0)}; if (length($0) <= min){min = length($0)} }; L=L+length($0); i++}}END{print i"\t"int(L/i)"\t"max"\t"min}')

PRER1=$(echo ${L1PRE} | cut -f1 -d ' ')
POSTR1=$(echo ${L1POST} | cut -f1 -d ' ')
PCTRETR1=$(bc <<< "scale=2;${POSTR1}/${PRER1}")

PRER2=$(echo ${L2PRE} | cut -f1 -d ' ')
POSTR2=$(echo ${L2POST} | cut -f1 -d ' ')
PCTRETR2=$(bc <<< "scale=2;${POSTR2}/${PRER2}")

DATE=$(date)
#Cleanup, organize files and generate report
printf "COGITO Sample: ${NAME}\n" >> ${illu1REP}
printf "Illumina read pair: 1\n" >> ${illu1REP}
printf "Analysis date: ${DATE}\n" >> ${illu1REP}
printf "Number of reads (post/pre): $(echo ${L1POST} | cut -f1 -d ' ')/$(echo ${L1PRE} | cut -f1 -d ' ')\n" >> ${illu1REP}
printf "Percent reads retained: $(echo ${PCTRETR1})\n" >> ${illu1REP}
printf "Mean read length (post/pre): $(echo ${L1POST} | cut -f2 -d ' ')/$(echo ${L1PRE} | cut -f2 -d ' ')\n" >> ${illu1REP}
printf "Maximum read length (post/pre): $(echo ${L1POST} | cut -f3 -d ' ')/$(echo ${L1PRE} | cut -f3 -d ' ')\n" >> ${illu1REP}
printf "Minimum read length (post/pre): $(echo ${L1POST} | cut -f4 -d ' ')/$(echo ${L1PRE} | cut -f4 -d ' ')\n" >> ${illu1REP}
printf "Adapter FASTQC control (post/pre): ${R1FQCPOSTA}/${R1FQCPREA}\n" >> ${illu1REP}
printf "Duplication FASTQC control (post/pre): ${R1FQCPOSTD}/${R1FQCPRED}\n" >> ${illu1REP}
printf "BBDUK version: ${BBDUKVER}\n" >> ${illu1REP}
printf "FASTQC version: ${FASTQCVER}\n" >> ${illu1REP}

printf "COGITO Sample: ${NAME}\n" >> ${illu2REP}
printf "Illumina read pair: 2\n" >> ${illu2REP}
printf "Analysis date: ${DATE}\n" >> ${illu2REP}
printf "Number of reads (post/pre): $(echo ${L2POST} | cut -f1 -d ' ')/$(echo  ${L2PRE} | cut -f1 -d ' ')\n" >> ${illu2REP}
printf "Percent reads retained: $(echo ${PCTRETR2})\n" >> ${illu2REP}
printf "Mean read length (post/pre): $(echo ${L2POST} | cut -f2 -d ' ')/$(echo  ${L2PRE} | cut -f2 -d ' ')\n" >> ${illu2REP}
printf "Maximum read length (post/pre): $(echo ${L2POST} | cut -f3 -d ' ')/$(echo  ${L2PRE} | cut -f3 -d ' ')\n" >> ${illu2REP}
printf "Minimum read length (post/pre): $(echo ${L2POST} | cut -f4 -d ' ')/$(echo  ${L2PRE} | cut -f4 -d ' ')\n" >> ${illu2REP}
printf "Adapter FASTQC control (post/pre): ${R2FQCPOSTA}/${R2FQCPREA}\n" >> ${illu2REP}
printf "Duplication FASTQC control (post/pre): ${R2FQCPOSTD}/${R2FQCPRED}\n" >> ${illu2REP}
printf "BBDUK version: ${BBDUKVER}\n" >> ${illu2REP}
printf "FASTQC version: ${FASTQCVER}\n" >> ${illu2REP}

# remove unzipped fastqc directories
    
rm -r ${FASTQCWDPRE}/${NAME}.raw.R1_fastqc
rm -r ${FASTQCWDPRE}/${NAME}.raw.R2_fastqc
rm -r ${FASTQCWDPOST}/${NAME}.clean.R1_fastqc
rm -r ${FASTQCWDPOST}/${NAME}.clean.R2_fastqc

echo "Pre-processing analysis complete for sample:" ${NAME}
)
													  








