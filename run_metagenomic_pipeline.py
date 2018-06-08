#$ -cwd
#$ -q sandbox
#$ -pe mpich 16
#$ -S /bin/bash

# Modules used
module load Megahit/v1.1.3
module load MetaBAT/v0.32.4
module load Bowtie2/v2.2.7
module load samtools/v1.3.1

# Variables
MEM=232000000000

for DIR in `find ./* -type d`; do

	DIR=${DIR#./}

	FILE1=${DIR}/*R1*.fastq.gz
	FILE2=${DIR}/*R2*.fastq.gz

	# Megahit assembly
	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o ${DIR}/Megahit --out-prefix ${DIR} 2> ${DIR}/${DIR}_assemble.log
	
	# Bowtie2 mapping
	REF=${DIR}/${DIR}.contigs.fa

	bowtie2 -x $REF --no-unal --very-sensitive -p $NSLOTS --mm \
		-1 $FILE1 -2 $FILE2 -S ${DIR}/${DIR}.sam 2> ${DIR}/${DIR}_bowtie2.log

	samtools view -Sb ${DIR}/${DIR}.sam > ${DIR}/${DIR}.bam
		rm ${DIR}/${DIR}.sam
	samtools sort ${DIR}/${DIR}.bam -o ${DIR}/${DIR}_sorted.bam
		rm ${DIR}/${DIR}.bam

	# MetaBat binning
	jgi_summarize_bam_contig_depths --outputDepth ${DIR}/${DIR}_depth.txt ${DIR}/${DIR}_sorted.bam
	
	metabat -i $REF -a ${DIR}/${DIR}_depth.txt -o ${DIR}/bin --sensitive -t $NSLOTS --saveCls --unbinned --keep

done


