#! /bin/bash
#PBS -l nodes=1:ppn=4:centos6,cput=24:00:00,walltime=48:00:00
#PBS -N assignment_1
#PBS -d /export/home/biostuds/2503736r/NGS/ASSIGNMENT
#PBS -m abe
#PBS -M 2503736r@student.gla.ac.uk
#PBS -q bioinf-stud
#
#MY DIRECTORIES (main, data location)
main_dir='/export/home/biostuds/2503736r/NGS/ASSIGNMENT'
HPC_data='/export/projects/polyomics/buzz/biostuds'
#
# RESOURCE FILES
reference_sequence='/export/projects/polyomics/Genome/Mus_musculus/mm10/genome/chr2.fa'
hs2index='/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2'
illumina_adapter='/export/projects/polyomics/biostuds/data/illumina_adapter.fa'
gtf_reference='/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf'
#
# MAKE SUBDIRS UNLESS THEY EXIST
fastqc_dir='/export/home/biostuds/2503736r/NGS/ASSIGNMENT/fastqc_dir'
hisat_dir='/export/home/biostuds/2503736r/NGS/ASSIGNMENT/hisat_dir'
stringtie_dir='/export/home/biostuds/2503736r/NGS/ASSIGNMENT/stringtie_dir'
mkdir -p ${fastqc_dir}
mkdir -p ${hisat_dir} 
mkdir -p ${stringtie_dir}
#
#loop going through all samples pertaining to the three groups and executing the cascade commands
#for each of them
for sample in s1.c2 s2.c2 s3.c2 s4.c2 s5.c2 s6.c2 s7.c2 s8.c2 s9.c2 s10.c2 s11.c2 s12.c2
do
	#VARIABLE DECLARATION
	raw_fastq="${sample}.fq"
	link_file="${sample}_link.fq"
	trim_adapt="${sample}.ta.fq"
	trim_quality="${sample}.tq.fq"
	sam_file="${sample}.sam"
	bam_file="${sample}.bam"
	bam_sorted="${sample}_sorted"
	gtf_file="${sample}.gtf"
	#
	#CREATION OF SOFT LINKS OF FILES 
	ln -s ${HPC_data}/${raw_fastq} ${main_dir}/${link_file}
	#
	#assess raw data quality with fastq software 
	fastqc ${link_file} --extract --outdir ${fastqc_dir}
	#Trim reads based on the QC assessment to improve read mapping
	scythe -o ${trim_adapt} -a ${illumina_adapter} -q sanger ${link_file}
	sickle se -f ${trim_adapt} -t sanger -o ${trim_quality} -q 10 -l 49 
	#
	#Alignment of the reads of RNA-seq to the reference genome (offset 33, Illumina 1.9 (superior to 1.8))
	hisat2 -p 4 --phred33 --rna-strandness R -x ${hs2index} -U ${trim_quality} -S ${sam_file} 
	#Sequence alignment mapping (output for transcript assembly)
	samtools view -b -S -o ${bam_file} ${sam_file}
	samtools sort ${bam_file} ${bam_sorted}
	#
	#removing bam, sam, ta and tq files
	rm ${sam_file} ${bam_file} ${trim_adapt} ${trim_quality}
	#
	#path to sample-specific subdirectory for stringtie results
	stringtie ${bam_sorted}.bam -p 4 --fr -b ${stringtie_dir}/${sample} -e -G -${gtf} -o ${stringtie_dir}/${sample}/${gtf_file}
	echo "${sample} /export/home/biostuds/2503736r/tx2dm/data/assignment/stringtie_dir/${sample}/${gtf_file}" >> sample_list.txt

done
#Preparation command to obtain the two differential expression tables (gene and transcript)
python2.7 /export/projects/polyomics/App/prepDE.py -i ${main_dir}/sample_list.txt



