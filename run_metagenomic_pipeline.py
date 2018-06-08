#$ -cwd
#$ -q sandbox
#$ -pe mpich 16
#$ -S /bin/bash

# Modules used
module load Megahit/v1.1.3

# Variables
MEM=232000000000

# Megahit assembly
for DIR in `find ./* -type d`; do

	FILE1=${DIR}/*R1*.fastq.gz
	FILE2=${DIR}/*R2*.fastq.gz

	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o ${DIR}/Megahit --out-prefix ${DIR} 2> ${DIR}/${DIR}_assemble.log

done


