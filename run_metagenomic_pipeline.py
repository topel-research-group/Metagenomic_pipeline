#$ -cwd
#$ -q Annotation-1
#$ -pe mpich 40
#$ -S /bin/bash

# Modules used
module load Megahit/v1.1.3

# Variables
MEM=232000000000

# Megahit assembly
for SAMPLE in {106..150}; do

	FILE1=../../00_data/Trimmed/P8352_${SAMPLE}/*R1.Pair.fastq.gz
	FILE2=../../00_data/Trimmed/P8352_${SAMPLE}/*R2.Pair.fastq.gz

	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o P8352_${SAMPLE} --out-prefix P8352_${SAMPLE} 2> P8352_${SAMPLE}_assemble.log

	done
