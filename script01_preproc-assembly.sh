#! /bin/bash

# -------------------------------------------------- #
# PRE-PROCCESSING AND ASSEMBLY OF SPONGE METAGENOMES #
# -------------------------------------------------- #

# NOTES:
## Input files "tmp-names-truseq.txt" and "tmp-names-nextera.txt" contain the file names corresponding to the samples/metagenomes with such adapters.
## Input file "tmp-names-all.txt" contains a list of the file names corresponding to all the samples/metagenomes.
## nohup command was used for each new execution, making a backup of each resulting "nohup.out" file to avoid overwriting. Execution: $ nohup ./{script} & 


#--- Quality checking of paired-end reads using FastQC and MultiQC:

fastqc *.fastq.gz -o quality --noextract -t 70
fastqc *.fastq -o quality -t 70
multiqc .


#--- Adapter trimming using Skewer:

# For Truseq adapters:
while IFS= read -r line
do
	skewer ${line}_R1.fastq.gz ${line}_R2.fastq.gz -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -r 0.03 -l 50 -k 16 -t 40 -z -o ${line}_truseq
done < tmp-names-truseq.txt

# For Nextera adapters:
while IFS= read -r line
do
	skewer ${line}_R1.fastq.gz ${line}_R2.fastq.gz -x CTGTCTCTTATACACATCT -y CTGTCTCTTATACACATCT -r 0.05 -l 50 -k 9 -t 40 -z -o ${line}_nextera
done < tmp-names-nextera.txt


#--- Quality filtering of paired-end reads using BBDuk:

while IFS= read -r line
do
	sh ../bbmap/bbduk.sh in=${line}_trimmed-pair1.fastq.gz in2=${line}_trimmed-pair2.fastq.gz out=${line}_q28-R1.fastq.gz out2=${line}_q28-R2.fastq.gz outm=${line}_failq28-R1.fastq.gz outm2=${line}_failq28-R2.fastq.gz maq=28 t=40
done < tmp-names-all.txt


#--- Assembly of metagenomes using MEGAHIT:

# Assemble paired-end reads to contigs without merging pairs:
gunzip -k *.fastq.gz

while IFS= read -r line
do
 	megahit -1 ${line}_q28-R1.fastq -2 ${line}_q28-R2.fastq -o ./megahit_${line} -t 35 # also works with ".gz" files
done < tmp-names-all.txt

# Verify that all assemblies finished well:
grep -c "DONE" nohup.out


#--- Quality assessment of metagenome assemblies using QUAST:

# Retrieve the contigs' files and change their names:
while IFS= read -r line
do
	mv megahit_${line}/final.contigs.fa ${line}_MH.fasta
done < tmp-names-all.txt

# Run assessment (QUAST doesn't use references, is faster than MetaQUAST):
while IFS= read -r line
do
	quast.py ${line}*.fasta -o ./quast-${line} -t 30
done < tmp-names-all.txt

# Retrieve the assessment reports and change their names:
while IFS= read -r line
do
	mv quast-${line}/report.html report-${line}.html
	mv quast-${line}/report.tsv report-${line}.tsv
done < tmp-names-all.txt

# Verify that execution finished well:
grep -c "Error" nohup.out
grep -c "Finished" nohup.out


