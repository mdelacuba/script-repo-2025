#! /bin/bash

# ----------------------------------------------------------- #
# FUNCTIONAL ANNOTATION AND FILTERING FROM SPONGE METAGENOMES #
# ----------------------------------------------------------- #

# NOTES:
## Input file "tmp-names-all.txt" contains a list of the file names corresponding to all the samples/metagenomes.
## nohup command was used for each new execution, making a backup of each resulting "nohup.out" file to avoid overwriting. Execution: $ nohup ./{script} & 
## Output file "table-all-func.txt" was used to manually convert sequence numbers in percentages and make the boxplots of "Figure 2" and "Supplementary Figure S3" (see "script07_plots-stats.r").
## Output files "list-orthologs-*.txt" and "list-csps-*.txt" were imported in R to make the Venn diagrams of "Supplementary Figure S1" and "Supplementary Figure S4", respectively (see "script07_plots-stats.r")
## Output files "all-subs-coldfuncs.txt" and "annot-subs.txt" were imported in R to make the heatmap of "Supplementary Figure S2" (see "script07_plots-stats.r").


#--- Gene prediction using Prodigal:

while IFS= read -r line
do
	prodigal -i ${line}_MH.fasta -o ${line}-cds.gbk -a ${line}-prot.faa -p meta # amino acid sequences
	prodigal -i ${line}_MH.fasta -d ${line}-nucl.fasta -p meta # nucleotide sequences
	prodigal -i ${line}_MH.fasta -f gff -o ${line}-cds.gff -p meta # save results in gff format
done < tmp-names-all.txt

# Verify that execution finished well:
grep -i "warning" nohop.out
grep -i "error" nohop.out


#--- Functional Annotation using eggNOG database and eggNOG-mapper:

while IFS= read -r line
do
	python ../eggnog-mapper-2.1.11/emapper.py -i ${line}-prot.faa -m diamond -o ${line}_annot_eggn --cpu 35 # rename annotation files as "${line}.annotations"
done < tmp-names-all.txt


#--- Assessing the gene percentage of selected functional groups:

# For cold adaptation functions:
echo -e "ID\tcsp\tchaperone\tosmoprotectant\tantioxidant\tantifreeze\tfatdesat\thsp\trepair" > table-all-func.txt

while IFS= read -r line
do
	a=$( echo -e "$line\t$(grep -v "^#" $line.annotations | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f8 | wc -l)" ) # Cold shock proteins (Csp)
	b=$( grep -v "^#" $line.annotations | grep -E -v -i "eukaryota|ABC transporter|MsrPQ" | grep -i "chaperone" | cut -f8 | wc -l ) # Chaperones
	c=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | grep -i "osmoprotectant" | cut -f8 | wc -l ) # Osmoprotectant-related proteins
	d=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f8 | wc -l ) # Antioxidant-related proteins
	e=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | grep -i "antifreeze" | cut -f8 | wc -l ) # Antifreeze proteins
	f=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" | cut -f8 | wc -l ) # Fatty acid desaturases
	g=$( grep -v "^#" $line.annotations | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" | cut -f8 | wc -l ) # Heat shock proteins (Hsp)
	h=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" | cut -f8 | wc -l ) # Nucleotide repair proteins
	echo -e "$line\t$a\t$b\t$c\t$d\t$e\t$f\t$g\t$h" >> table-all-func.txt
done < tmp-names-all.txt

# For metabolism functions:
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f7 | grep "I" | sort | uniq > list-lipid-cogs.txt # Lipid transport and metabolism (LTM)
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f7 | grep "C" | sort | uniq > list-energ-cogs.txt # Energy production and conversion (EPC)
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f7 | grep "E" | sort | uniq > list-amino-cogs.txt # Amino acid transport and metabolism (ATM)
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f7 | grep "G" | sort | uniq > list-carbo-cogs.txt # Carbon transport and metabolism (CTM)

echo -e "ID\tlipid\t\energ\tamino\tcarbo"

for categ in lipid energ amino carbo
do
	echo -e "\t$categ"
	for file in *.annotations
	do
		a=$( grep -v "^#" $file | grep -v -i "eukaryota" | cut -f 1,2,7 | grep -wf list-$categ-cogs.txt | wc -l )
		echo -e "$file\t$a" # Paste manually to "table-all-func.txt"
	done
done

# For resistance and conjugation functions:
echo -e "ID\tresis\t\ices"

for file in *.annotations
do
	a=$( grep -Evi "^#|eukaryota" $file | grep -iE "resistance|antibiotic|drug|antimicrobial|heavy metal|heavy-metal" | grep -Ev "bacteriophage resistance|PhoQ|Antitoxin|nucleosome-like|oxidative stress|low-pH|cold resistance|resistance to hypoosmotic shock|stress resistance|macrophage|phage-resistance|extreme acid resistance|Toll-Interleukin 1-resistance|Serum resistance|TraT complement|Ultraviolet light resistance|Ultra-violet resistance" | cut -f8 | wc -l) # Metal and antibiotic resistance (MAR)
	b=$( grep -Evi "^#|eukaryota|CRISPR|phage" $file | grep -iE "integrase|relaxase|relaxosome|conjugative|Type IV secretion|Type-IV|Type IV secretory|T4SS|tyrosine recombinase|ice cds17|TrbI-like|TrbI family|TraG family|unidirectional conjugation|Tra gene|TraG-like|TraB family|TraC-like" | cut -f8 | wc -l ) # Machinary of integrative and conjugative elements (ICE)
	echo -e "$file\t$a\t$b" # Paste manually to "table-all-func.txt"
done

# Count total genes annotated:
while IFS= read -r line
do
	a=$( grep -v "^#" $line.annotations | grep -v -i "eukaryota" | wc -l ) # Without Eukarya
	echo -e "$a\t$line" # Add manually to table-all-func.txt
done < tmp-names-all.txt


#--- Assessing the shared and exclusive gene orthologs among sponge metagenomes of different environments:

# Separate files manually into per-environment folders:
mkdir Antarctic Temperate Tropical

#  Search all orthologs IDs per each environment:
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f2 | sort | uniq > list-orthologs-anta.txt # Into "Antarctic" folder
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f2 | sort | uniq > list-orthologs-temp.txt # Into "Temperate" folder
grep -v "^#" *.annotations | grep -v -i "eukaryota" | cut -f2 | sort | uniq > list-orthologs-trop.txt # Into "Tropical" folder

# Search Cold shock protein (Csp) orthologs IDs per each environment:
grep -v "^#" *.annotations | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f2 | sort | uniq > list-csps-anta.txt # Into "Antarctic" folder
grep -v "^#" *.annotations | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f2 | sort | uniq > list-csps-temp.txt # Into "Temperate" folder
grep -v "^#" *.annotations | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f2 | sort | uniq > list-csps-trop.txt # Into "Tropical" folder


#--- Assessing the presence/absence of genes for cold adaptation from subsampling of the quality-filtered reads:

# Subsample metagenomic reads to the sample with the lowest sequence number using Seqtk:
for file in *.fastq; do seqtk sample -s 1000 ${file} 266595 > ${file}_subs; done # SRA ID SRR7783610

# Follow the assembly, gene prediction, and functional annotation as before for non-subsampled metagenomes...

# ...Then, make tables of orthologs by metagenome:
for file in *.annotations
do
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | wc -l)) <(grep -v "^#" $file | grep -E -v -i "eukaryota|toxin" | grep -i "cold" | cut -f2) > subs-csps-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -E -v -i "eukaryota|ABC transporter|MsrPQ" | grep -i "chaperone" | wc -l)) <(grep -v "^#" $file | grep -E -v -i "eukaryota|ABC transporter|MsrPQ" | grep -i "chaperone" | cut -f2) > subs-chap-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "osmoprotectant" | wc -l)) <(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "osmoprotectant" | cut -f2) > subs-osmo-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "antioxidant" | wc -l)) <(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f2) > subs-anox-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "antifreeze" | wc -l)) <(grep -v "^#" $file | grep -v -i "eukaryota" | grep -i "antifreeze" | cut -f2) > subs-anfr-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" | wc -l)) <(grep -v "^#" $file | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" | cut -f2) > subs-fatd-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" | wc -l)) <(grep -v "^#" $file | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" | cut -f2) > subs-hsps-$file.txt
	paste -d "\t" <(yes "$file 1" | head -$(grep -v "^#" $file | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" | wc -l)) <(grep -v "^#" $file | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" | cut -f2) > subs-repa-$file.txt
done

# Join all files in a single large table:
sed -i "s/ /\t/g" subs-*.txt
echo -e "Sample\tFreq\tSeed_ortholog" >> all-subs-coldfuncs.txt
cat subs-* >> all-subs-coldfuncs.txt
sed -i "s/_subs.annotations//g" all-subs-coldfuncs.txt

# Reorder orthologs according to each function for cold adaptation:
while IFS= read -r line
do
	a=$( grep -w -m1 "$line" subs-*.txt | head -n1 | cut -d "-" -f2 )
	echo -e "$line\t$a" >> tmp.txt
done < tmp-subs-orthologs.txt # List of orthologs exported from R

sort -k 2,2 tmp.txt > annot-subs.txt
rm tmp.txt

# Replace short names with complete ones:
sed -i "s/anfr/Antifreeze/g" annot-subs.txt
sed -i "s/anox/Antioxidant/g" annot-subs.txt
sed -i "s/chap/Chaperone related/g" annot-subs.txt
sed -i "s/csps/Cold-shock related/g" annot-subs.txt
sed -i "s/fatd/Fatty acid desaturase/g" annot-subs.txt
sed -i "s/hsps/Heat-shock related/g" annot-subs.txt
sed -i "s/osmo/Osmoprotectant related/g" annot-subs.txt
sed -i "s/repa/Nucleotide repair/g" annot-subs.txt


