#$ -cwd
#$ -q Annotation-3
#$ -pe mpich 40
#$ -S /bin/bash

# Modules used
module load Megahit/v1.1.3
module load MetaBAT/v0.32.4
module load Bowtie2/v2.2.7
module load samtools/v1.3.1
module load CheckM/v1.0.7

# Variables
MEM=232000000000

# Print some informative error meassages
err() {
	echo "$1 exited unexpectedly";
	break;
}

# Function for checking the exit code of a child process
checkExit() {
if [ "$1" == "0" ]; then
	echo "[Done] $2 `date`";
	else
	err "[Error] $2 returned non-0 exit code $1";
	fi
}


for DIR in `find ./* -maxdepth 0 -type d`; do
	# Make sure this sample has not been analysed already
	if [ -f ${DIR}/pipeline.done ]; then
		continue
	else

	DIR=${DIR#./}

	FILE1=${DIR}/*R1*.fastq.gz
	FILE2=${DIR}/*R2*.fastq.gz

	# Megahit assembly
	megahit -m $MEM -1 $FILE1 -2 $FILE2 -t $NSLOTS -o ${DIR}/Megahit --out-prefix ${DIR} 2> ${DIR}/${DIR}_assemble.log
		checkExit $? "megahit"

	# Bowtie2 mapping
	mkdir ${DIR}/Bowtie2

	REF=${DIR}/Megahit/${DIR}.contigs.fa
	DB=${DIR}/Bowtie2/${DIR}.contigs
	bowtie2-build --threads $NSLOTS $REF $DB 2> ${DIR}/${DIR}_bowtie2.log
		checkExit $? "bowtie2-build"

	bowtie2 -x $DB --no-unal --very-sensitive -p $NSLOTS --mm \
		-1 $FILE1 -2 $FILE2 -S ${DIR}/Bowtie2/${DIR}.sam 2> ${DIR}/${DIR}_bowtie2.log
		checkExit $? "bowtie2"

	samtools view -Sb ${DIR}/Bowtie2/${DIR}.sam > ${DIR}/Bowtie2/${DIR}.bam
		rm ${DIR}/Bowtie2/${DIR}.sam
		checkExit $? "samtools"
	samtools sort ${DIR}/Bowtie2/${DIR}.bam -o ${DIR}/Bowtie2/${DIR}_sorted.bam
		rm ${DIR}/Bowtie2/${DIR}.bam
		checkExit $? "samtools"

	# MetaBat binning
	mkdir ${DIR}/Metabat
	jgi_summarize_bam_contig_depths --outputDepth ${DIR}/Bowtie2/${DIR}_depth.txt \
					${DIR}/Bowtie2/${DIR}_sorted.bam 2> ${DIR}/${DIR}_metabat.log
		checkExit $? "jgi_summarize_bam_contig_depths"
	
	metabat -i $REF -a ${DIR}/Bowtie2/${DIR}_depth.txt -o ${DIR}/Metabat/bin \
		--sensitive -t $NSLOTS --saveCls --unbinned --keep 2> ${DIR}/${DIR}_metabat.log
		checkExit $? "metabat"

	checkm lineage_wf -t $NSLOTS -x fa ${DIR}/Metabat ${DIR}/CheckM
		checkExit $? "checkm"

	# Mark the sample as finished
	touch ${DIR}/pipeline.done
	fi 
done


