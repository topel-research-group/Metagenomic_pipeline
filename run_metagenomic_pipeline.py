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

	FILE1=${DIR}/*R1*.fastq.gz
	FILE2=${DIR}/*R2*.fastq.gz

	# Megahit assembly
	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o ${DIR}/Megahit --out-prefix ${DIR} 2> ${DIR}/${DIR}_assemble.log
	
	# Bowtie2 mapping
	DB=all_contigs_10000
	FILE1=../00_data/Trimmed/P8352_${SAMPLE}/*R1.Pair.fastq.gz
	FILE2=../00_data/Trimmed/P8352_${SAMPLE}/*R2.Pair.fastq.gz

	bowtie2 -x $DB --no-unal --very-sensitive -p $NSLOTS --mm \
		-1 $FILE1 -2 $FILE2 -S P8352_${SAMPLE}_${DB}.sam 2> P8352_${SAMPLE}_${DB}.log

	samtools view -Sb P8352_${SAMPLE}_${DB}.sam > P8352_${SAMPLE}_${DB}.bam
		rm P8352_${SAMPLE}_${DB}.sam
	samtools sort P8352_${SAMPLE}_${DB}.bam -o P8352_${SAMPLE}_${DB}_sorted.bam
		rm P8352_${SAMPLE}_${DB}.bam

	# MetaBat binning
	REF=${DIR}.contigs.fa
	jgi_summarize_bam_contig_depths --outputDepth depth.txt *_sorted.bam
	metabat -i $REF -a depth.txt -o bin --sensitive -t 40 --saveCls --unbinned --keep

done


