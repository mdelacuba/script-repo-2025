#! /bin/bash

# -------------------------------------------------- #
# GENE ABUNDANCE NORMALIZATION OF SPONGE METAGENOMES #
# -------------------------------------------------- #

# NOTES:
## Input file "tmp-names-all.txt" contains a list of the file names corresponding to all the samples/metagenomes.
## nohup command was used for each new execution, making a backup of each resulting "nohup.out" file to avoid overwriting. Execution: $ nohup ./{script} &
## Output file "table-cogs-abnorm.txt" was imported in R to make the barplot of "Figure 1B".
## Output file "tmp-all-rpkg.txt" was imported in R to construct an ordination matrix and make the NMDS plot of "Figure 1C" (see "script07_plots-stats.r"). 
## Output files "table-all-abnorm.txt" and "table-shannon.txt" were imported in R to make boxplots similar to "Figure 2" and "Supplementary Figure S3" but were not added to the main manuscript or the additional files (see "script07_plots-stats.r").


#--- Mapping reads against genes for calculating gene abundace using Bowtie 2 and BBMap:

while IFS= read -r line
do
        mkdir bw2-${line}
	../bowtie2-2.5.1-linux-x86_64/bowtie2-build ${line}-nucl.fasta bw2-${line}/${line} --large-index --threads 30 --seed 1000
	../bowtie2-2.5.1-linux-x86_64/bowtie2 -x bw2-${line}/${line} -1 ${line}_q28-R1.fastq -2 ${line}_q28-R2.fastq -S ${line}.sam --reorder --threads 30 --very-sensitive --seed 1000
	../bbmap/pileup.sh in=${line}.sam rpkm=${line}-rpkm.txt
done < tmp-names-all.txt


#--- Abundance normalization using MicrobeCensus:

# Calculate Average Genome Size (AGS) and Genome Equivalents (GE):
while IFS= read -r line
do
	run_microbe_census.py -t 30 ${line}_q28-R1.fastq,${line}_q28-R2.fastq ${line}-AGS.out
done < tmp-names-all.txt

# Extract GE:
grep "genome_equivalents:" *AGS.out | sed "s/-AGS.out:genome_equivalents://g" > GE-all.txt

# Calculate RPKG:
while IFS= read -r line
do
	paste -d "\t" ${line}-rpkm.txt <(yes '' | sed 4q && echo "GE" && yes "$(grep -w "${line}" GE-all.txt | cut -f2)" | head -$(grep -v "^#" ${line}-rpkm.txt | wc -l)) > ${line}-ge.txt
	paste -d "\t" ${line}-ge.txt <(yes '' | sed 4q && echo "RPKG" && grep -v "^#" ${line}-ge.txt | awk '{print $5/($2/1000)/$9}') > ${line}-rpkg.txt
	echo "$(head -n5 ${line}-rpkg.txt | tail -n1)" > ${line}-annot-rpkg.txt && grep -wFf list-annot-${line}.txt ${line}-rpkg.txt >> ${line}-annot-rpkg.txt
done < tmp-names-all.txt

# Make table joining abundance and annotation info:
while IFS= read -r line
do
	paste -d "\t" ${line}-annot-rpkg.txt <(grep -v "^##" $line.annotations) > $line-annot-fulltab.txt
done < tmp-names-all.txt


#--- Assessing COG categories abundance per metagenome:

# Extract COG abundances to a table:
array=( - A B C D E F G H I J K L M N O P Q R S T U V W X Y Z )

echo -e "ID\tCOG\tcount" > table-cogs-abnorm-tmp.txt

while IFS= read -r line
do
	for cog in "${array[@]}"
	do
		echo -e "$line\t$cog\t$(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F'\t' '{print $17"\t"$10}' | grep "$cog" | awk '{sum+=$2} END {print sum}')" >> table-cogs-abnorm-tmp.txt
	done
done < tmp-names-all.txt

# Rename unassigned COG annotations:
sed -i 's/\t-\t/\tUnassigned\t/g' table-cogs-abnorm-tmp.txt

# Isolate COG letters from table:
cut -f2 table-cogs-abnorm-tmp.txt > tmp-cog-list.txt

# Extract names of COG categories and groups from external table:
echo -e "cog\tdescription\tgroup" > tmp-categs.txt

while IFS= read -r line
do
	grep -w "$line" list-cog-categories.txt >> tmp-categs.txt
done < tmp-cog-list.txt

# Join columns in a new table:
paste -d "\t" table-cogs-abnorm-tmp.txt tmp-categs.txt > table-cogs-abnorm.txt

# Extract the belonging to habitat (Antarctic sponge, tropical sponge, temperate sponge, or Antarctic seawater):
cut -f1 table-cogs-abnorm.txt > tmp-ids.txt

while IFS= read -r line
do
	a=$( grep -w "$line" ../metadata.csv | cut -f2) # Paste manually to "table-cogs-abnorm.txt"
	echo -e "$a\t$line"
done < tmp-ids.txt

# Remove temporal files:
rm tmp-ids.txt
rm tmp-categs.txt
rm table-cogs-abnorm-tmp.txt


#--- Assessing the gene ortholog composition of metagenomes:

# Make a table of orthologs abundances by metagenome:
while IFS= read -r line
do
	paste -d "\t" <(yes "$line" | head -$(wc -l < $line-annot-fulltab.txt) | tail -n +2) <(awk '{print $12"\t"$10}' $line-annot-fulltab.txt | tail -n +2) > tmp-rpkg-$line.txt
done < tmp-names-all.txt

# Join all files in a single large table:
echo -e "Sample\tSeed_ortholog\tRPKG" >> tmp-all-rpkg.txt
cat tmp-rpkg-* >> tmp-all-rpkg.txt
mv tmp-rpkg-* tmp_rpkg_files/


#--- Assessing the gene abundance of selected functional groups related to cold adaptation, metabolism, resistance, and conjugation:

# Make table from RPKG calculation:
echo -e "ID\tcsp\tchaperone\tosmoprotectant\tantioxidant\tantifreeze\tfatdesat\thsp\trepair\tlipid\tenerg\tamino\tcarbo\tresis\tices" > table-all-abnorm.txt

while IFS= read -r line
do
	paste -d "\t" <(echo -e "$line\t$(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f10 | awk '{sum+=$1} END {print sum}')") <(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|ABC transporter|MsrPQ" | grep -i "chaperone" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "osmoprotectant" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "antifreeze" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}' | grep "I" | awk '{sum+=$2} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "C" | awk '{sum+=$2} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "E" | awk '{sum+=$2} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "G" | awk '{sum+=$2} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -iE "resistance|antibiotic|drug|antimicrobial|heavy metal|heavy-metal" | grep -Ev "bacteriophage resistance|PhoQ|Antitoxin|nucleosome-like|oxidative stress|low-pH|cold resistance|resistance to hypoosmotic shock|stress resistance|macrophage|phage-resistance|extreme acid resistance|Toll-Interleukin 1-resistance|Serum resistance|TraT complement|Ultraviolet light resistance|Ultra-violet resistance" | cut -f10 | awk '{sum+=$1} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -Evi "eukaryota|CRISPR|phage" | grep -iE "integrase|relaxase|relaxosome|conjugative|Type IV secretion|Type-IV|Type IV secretory|T4SS|tyrosine recombinase|ice cds17|TrbI-like|TrbI family|TraG family|unidirectional conjugation|Tra gene|TraG-like|TraB family|TraC-like" | cut -f10 | awk '{sum+=$1} END {print sum}') >> table-all-abnorm.txt
done < tmp-names-all.txt
	
# Make table from Shannon diversity index calculation:
echo -e "ID\tcsp\tchaperone\tosmoprotectant\tantioxidant\tantifreeze\tfatdesat\thsp\trepair\tlipid\tenerg\tamino\tcarbo\tresis\tices" > table-shannon.txt

while IFS= read -r line
do
	paste -d "\t" <(echo -e "$line\t$(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}')") <(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|ABC transporter|MsrPQ" | grep -i "chaperone" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "osmoprotectant" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -i "antifreeze" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "I" | cut -f2 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "C" | cut -f2 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "E" | cut -f2 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | awk -F '\t' '{print $17"\t"$10}'  | grep "G" | cut -f2 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -v -i "eukaryota" | grep -iE "resistance|antibiotic|drug|antimicrobial|heavy metal|heavy-metal" | grep -Ev "bacteriophage resistance|PhoQ|Antitoxin|nucleosome-like|oxidative stress|low-pH|cold resistance|resistance to hypoosmotic shock|stress resistance|macrophage|phage-resistance|extreme acid resistance|Toll-Interleukin 1-resistance|Serum resistance|TraT complement|Ultraviolet light resistance|Ultra-violet resistance" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') <(grep -v "^#" $line-annot-fulltab.txt | grep -Evi "eukaryota|CRISPR|phage" | grep -iE "integrase|relaxase|relaxosome|conjugative|Type IV secretion|Type-IV|Type IV secretory|T4SS|tyrosine recombinase|ice cds17|TrbI-like|TrbI family|TraG family|unidirectional conjugation|Tra gene|TraG-like|TraB family|TraC-like" | cut -f10 | awk '{sum += $1; value[NR] = $1} END {for (i = 1; i <= NR; i++) {a = value[i]/sum; printf("%0.4f\n", a)}}' | awk '{a = $1*log($1); printf("%0.4f\n", a)}' | sed "s/nan/0/g" | awk '{sum+=-($1)} END {print sum}') >> table-shannon.txt
done < tmp-names-all.txt


