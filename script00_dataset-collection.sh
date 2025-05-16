#! /bin/bash

# ---------------------------------------- #
# DATASET COLLECTION OF SPONGE METAGENOMES #
# ---------------------------------------- #

# NOTE:
## Full metadata available in "Supplementary Table S1".


#--- Downloading info on sponge metagenomes from tropical, temperate, and Antarctic environments, linked to 17 published studies from the NCBI web page:

# ACCESSIONS:
: '
##	_Reference study_		_Accession_
##	Engelberts et al., 2020	PRJNA555144
##	Glasl et al., 2020		PRJNA594068
##	Knobloch et al., 2020		PRJNA508373
##	Rust et al., 2020		PRJNA603662
##	Karimi et al., 2017		PRJEB11585
##	Gao et al., 2017		PRJNA383957
##	Bayer et al., 2020		PRJNA613976
##	Rubin-Blum et al., 2019 (1)	PRJNA475438
##	Rubin-Blum et al., 2019 (2)	PRJNA475442
##	Botté et al., 2019		PRJNA488959
##	Robbins et al., 2020		PRJNA602572
##	Storey et al., 2020		PRJNA515312
##	Horn et al., 2016		PRJNA318959
##	Li et al., 2014		SRA052801
##	Taylor et al., 2020		PRJNA589708
##	Slaby et al., 2017		PRJNA366444
##	Moreno-Pino et al., 2024	PRJNA874040
'

# DOWNLOADS:
## SRA accessions lists, input files "SRR_Acc_List-$study.txt" 
## SRA run metadata, input files "SraRunTable-$study.csv"


#--- Extracting metagenome run IDs from SRA accession lists to each sponge species:

# For accession lists from a single sponge species with several replicates (random selection of a replicate):
array1=( Engelberts2020 Glasl2020 Knobloch2020 Rust2020 Karimi2017 Gao2017 Bayer2020 Rubin-Blum2019-1 Rubin-Blum2019-2 )

for study in "${array1[@]}"
do
	a=$( rl -c1 SRR_Acc_List-$study.txt  )
	echo -e "$a\t$study" >> selected-IDs.txt
done

# For accession lists from several sponge species with several replicates (random selection of a replicate):
b=$( awk "FNR>=1 && FNR<=3" SRR_Acc_List-Botté2019.txt | rl -c1 )
c=$( awk "FNR>=4 && FNR<=6" SRR_Acc_List-Botté2019.txt | rl -c1 )
d=$( awk "FNR>=7 && FNR<=9" SRR_Acc_List-Botté2019.txt | rl -c1 )
echo -e "$b\tBotté2019" >> selected-IDs.txt
echo -e "$c\tBotté2019" >> selected-IDs.txt
echo -e "$d\tBotté2019" >> selected-IDs.txt

e=$( awk "FNR>=1 && FNR<=4" SRR_Acc_List-Robbins2020.txt | rl -c1 )
f=$( awk "FNR>=5 && FNR<=8" SRR_Acc_List-Robbins2020.txt | rl -c1 )
g=$( awk "FNR>=9 && FNR<=12" SRR_Acc_List-Robbins2020.txt | rl -c1 ) 
h=$( awk "FNR>=13 && FNR<=15" SRR_Acc_List-Robbins2020.txt | rl -c1 )
i=$( awk "FNR>=16 && FNR<=19" SRR_Acc_List-Robbins2020.txt | rl -c1 )
j=$( awk "FNR>=20 && FNR<=23" SRR_Acc_List-Robbins2020.txt | rl -c1 )
echo -e "$e\tRobbins2020" >> selected-IDs.txt
echo -e "$f\tRobbins2020" >> selected-IDs.txt
echo -e "$g\tRobbins2020" >> selected-IDs.txt
echo -e "$h\tRobbins2020" >> selected-IDs.txt
echo -e "$i\tRobbins2020" >> selected-IDs.txt
echo -e "$j\tRobbins2020" >> selected-IDs.txt

k=$( awk "FNR>=2 && FNR<=5" SRR_Acc_List-Storey2020.txt | rl -c1 )
echo -e "$k\tStorey2020" >> selected-IDs.txt

# For accession lists from several/single sponge species with a single replicate (non-random selection):
array2=( Horn2016 Li2014 Taylor2020 Moreno-Pino2024) # Except SRR25758459 from Moreno-Pino2024

for study in "${array2[@]}"
do
	l=$( cat SRR_Acc_List-$study.txt )
	echo -e "$l\t$study" >> selected-IDs.txt
done

# For accession lists from a single sponge species with several replicates (targeted selection of a replicate according to metadata "SraRunTable-$study.csv"):
m=$( head -n1 SRR_Acc_List-Slaby2017.txt )
echo -e "$m\tSlaby2017" >> selected-IDs.txt

n=$( tail -n1 SRR_Acc_List-Storey2020.txt )
echo -e "$n\tStorey2020" >> selected-IDs.txt


#--- Download metagenome datasets from SRA database:

while IFS= read -r line
do
	prefetch $line -X 40G -O SRA-raw
	cd SRA-raw/$line/
	fasterq-dump $line.sra -S -x -O . 2>> ../details-SRA-raw.txt
	echo "----------------------------------------------------" >> ../details-SRA-raw.txt
	cd ../../
done < selected-IDs.txt


