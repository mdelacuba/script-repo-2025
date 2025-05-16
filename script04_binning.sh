#! /bin/bash

# ----------------------------------------------------------- #
# BINNING OF ANTARCTIC METAGENOMES FROM SPONGES AND SEAWATER) #
# ----------------------------------------------------------- #

# NOTES:
## Input file "tmp-names-all.txt" contains a list of the file names corresponding to all the samples/metagenomes.
## nohup command was used for each new execution, making a backup of each resulting "nohup.out" file to avoid overwriting. Execution: $ nohup ./{script} &


#--- Binning and related processes using MetaWRAP:

# Uncompress high-quality paired-end reads:
gunzip -k *.gz

while IFS= read -r line
do
# Binning:
	metawrap binning -o out_${line}_binning -t 30 -m 50 -a ${line}_MH.fasta --metabat2 --maxbin2 reads_q28/${line}_q28_*.fastq
# Bin refinement:	
	metawrap bin_refinement -o bin_refinement_${line} -t 40 -c 70 -x 5 -A out_${line}_binning/maxbin2_bins/ -B out_${line}_binning/metabat2_bins 
# Bin quantification:	
	metawrap quant_bins -t 50 -b bin_refinement_${line}/metawrap_70_5_bins -o quant_bins_${line} -a ${line}_MH.fasta reads_q28/${line}_q28_*.fastq
# Bin reassembly:	
	metawrap reassemble_bins -o bin_reassembly_${line} -1 reads_q28/${line}_q28_1.fastq -2 reads_q28/${line}_q28_2.fastq -t 50 -m 80 -c 70 -x 5 -b bin_refinement_${line}/metawrap_70_5_bins

# Change names of quantication output files:
        cd quant_bins_${line}/
	rename "s/bin_abundance_table/${line}_bin_abundance/g" *.tab
        cd ../
        	
# Change names of reassembled bins (high-quality bins):
	cd bin_reassembly_${line}/reassembled_bins/
	rename "s/bin./${line}_bin./g" *.fa
	cd ../../
	
# Extract information of the high-quality bins:
	awk -v var="$line" '{print $0 "\t"var}' reassembled_bins_$line.stats > tmp_$line.stats
done < tmp-names-all.txt

# Join information of the high-quality bins:
echo -e "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tsample" > bin_stats_join.txt
tail -n +2 tmp_*.stats >> bin_stats_join.txt	

#--- Inspecting bin numbers:

# Initial bins number:
tail -n +2 bin_refinement_GM2034-*/*2_bins.stats | cut -f1 | grep -v -E "==>|unbinned" | grep "bin" | wc -l

# High-quality bins (MAGs) number:
tail -n +2 bin_refinement_GM2034-*/metawrap*.stats | cut -f1 | grep -v "==>" | grep "bin" | wc -l


