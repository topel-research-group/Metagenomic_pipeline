#$ -cwd
#$ -q sandbox
#$ -pe mpich 16
#$ -S /bin/bash

# Modules used
module load Megahit/v1.1.3
module load MetaBAT/v0.32.4
module load Bowtie2/v2.2.7
module load samtools/v1.3.1
module load CheckM/v1.0.7

# Variables
MEM=232000000000

for DIR in `find ./* -maxdepth 0 -type d`; do
	# Make sure this sample has been analysed
	if [ -f ${DIR}/pipeline.done ]; then
		break
	else

	DIR=${DIR#./}

	FILE1=${DIR}/*R1*.fastq.gz
	FILE2=${DIR}/*R2*.fastq.gz

	# Megahit assembly
#	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o ${DIR}/Megahit --out-prefix ${DIR} 2> ${DIR}/${DIR}_assemble.log
	
	# Bowtie2 mapping
#	mkdir ${DIR}/Bowtie2

	REF=${DIR}/Megahit/${DIR}.contigs.fa
#	DB=${DIR}/Bowtie2/${DIR}.contigs
#	bowtie2-build --threads $NSLOTS $REF $DB 2> ${DIR}/${DIR}_bowtie2.log

#	bowtie2 -x $DB --no-unal --very-sensitive -p $NSLOTS --mm \
#		-1 $FILE1 -2 $FILE2 -S ${DIR}/Bowtie2/${DIR}.sam 2> ${DIR}/${DIR}_bowtie2.log

#	samtools view -Sb ${DIR}/Bowtie2/${DIR}.sam > ${DIR}/Bowtie2/${DIR}.bam
#		rm ${DIR}/Bowtie2/${DIR}.sam
#	samtools sort ${DIR}/Bowtie2/${DIR}.bam -o ${DIR}/Bowtie2/${DIR}_sorted.bam
#		rm ${DIR}/Bowtie2/${DIR}.bam

	# MetaBat binning
#	mkdir ${DIR}/Metabat
#	jgi_summarize_bam_contig_depths --outputDepth ${DIR}/Bowtie2/${DIR}_depth.txt \
#					${DIR}/Bowtie2/${DIR}_sorted.bam 2> ${DIR}/${DIR}_metabat.log
	
#	metabat -i $REF -a ${DIR}/Bowtie2/${DIR}_depth.txt -o ${DIR}/Metabat/bin \
#		--sensitive -t $NSLOTS --saveCls --unbinned --keep 2> ${DIR}/${DIR}_metabat.log

	checkm lineage_wf -t $NSLOTS -x fa ${DIR}/Metabat/bin ${DIR}/CheckM
	
	# Mark the sample as finished
	touch ${DIR}/pipeline.done
	fi 
done


