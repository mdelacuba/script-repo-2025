#! /bin/bash

# --------------------------------------------------------- #
# HGT DETECTION IN MAGs FROM ANTARCTIC SPONGES AND SEAWATER #
# --------------------------------------------------------- #

# NOTES:
## HT: Horizontally transferred
## Csp: Cold shock proteins
## Hsp: Heat shock proteins
## CAF: Cold adaptacion functions
## Input file "list-bins.txt" contains a list of the file names corresponding to all MAGs.
## Input file "list-bins-hgts.txt" contains a list of the file names corresponding to MAGs that present HT genes.
## Input files "list-recipient-hab-onlySWtoSP.txt" and "list-recipient-tax-onlySWtoSP.txt" contain a manually selected list of the recipient MAGs from sponges presenting HT genes from donor MAGs from seawater.
## Input file "list_hgts_all.txt" contains a list of the file names corresponding to MAGs involved in HGT.
## nohup command was used for each new execution, making a backup of each resulting "nohup.out" file to avoid overwriting. Execution: $ nohup ./{script} &
## Output file "tab-count-hgts.txt" was used to calculate percentages of HT genes and manually construct the customized table "tab-corr.txt". This new output file was imported in R to make the correlation plots of "Supplementary Figure S5" (see "script07_plots-stats.r").
## Output file "df-count-hgts.txt" was manually modified to leave only % HT genes, resulting in the dataframe file "df-count-hgts.txt". This new output file was imported in R to make the boxplots of "Figure 4" (see "script07_plots-stats.r").
## Output file "tab-count-donors.txt" was adapted and improved using the metadata of MAGs, resulting in the file "plotting_donors_recipients.csv". This new output file was imported in R to make the dotplot of "Figure 5" see "script07_plots-stats.r".
## Output files "tab-count-donors.txt" was also used to calculate gene numbers of undetermined taxonomic distances of HGTs, resulting in the file "count_donors_NA.csv". This new output file was imported in R to make the barplot of "Figure 5" (see "script07_plots-stats.r").
## Output files "recipient-hgt-hab-onlySWtoSP.annotations and "recipient-hgt-tax-onlySWtoSP.annotations" were manually merged and improved, resulting in the file "recipient-hgt-onlySWtoSP.annotations". This new output file was imported in R to make the dotplot of "Supplementary Figure S6" (see "script07_plots-stats.r"). 
## Output files "list-csps-*.txt" and "list-anox-*.txt" were imported in R to make Venn diagrams of "Figure 6" (see "script07_plots-stats.r").
## Output files "table-csp-all.txt" and "table-anox-all.txt" were imported in R to make alluvial plots of "Figure 6" (see "script07_plots-stats.r").


#--- HGT detection between recipient MAGs (from sponges or seawater) and donor reference genomes using RefSeq database and HGTector2:

# Detect putative HT genes:
while IFS= read -r line
do
	hgtector search -i input_proteins -o out_hgtector_search -m diamond -d ../hgtector_db_dir/diamond/db -t ../hgtector_db_dir/taxdump -p 35
	hgtector analyze -i out_hgtector_search/${line}_prot.tsv -o out_hgtector_analyze_${line} -t ../hgtector_db_dir/taxdump --donor-name
done < list-bins.txt

# Extract putative HT genes (generate single-line fasta proteins):
while IFS= read -r line
do
	cut -f1 tabs-hgts-hgtector/${line}_hgtector.txt > tmp-${line}-hgts.txt
	grep -wA1 -f tmp-${line}-hgts.txt proteins_singleline/${line}_prot.fasta >> ${line}_hgts.fasta
	head -n5 ${line}_func.annotations > ${line}_hgts.annotations && grep -w -f tmp-${line}-hgts.txt ${line}_func.annotations >> ${line}_hgts.annotations
done < list-bins-hgts.txt

sed -i.bak '/--/d' *hgts.fasta


#--- Functional annotation of HT genes using eggNOG database and eggNOG-mapper (from HGTector2 outputs):

while IFS= read -r line
do
	python ../../eggnog-mapper-2.1.11/emapper.py -i ${line}_hgts.fasta -m diamond -o ${line}_hgts --go_evidence all --cpu 30 --output_dir ../annotated-hgtector
done < list-bins-hgts.txt


#--- Filtering the annotations of selected functional groups presenting HT genes (from HGTector2 outputs):

# For cold adaptation functions (CAF):
while IFS= read -r line
do
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -i "cold" >> tmp-${line}-csps.annotations # Csp
	grep -v "^#" ${line}_hgts.annotations | grep -E -v -i "eukaryota|ABC transporter|MsrPQ|MeaI" | grep -i "chaperone" >> tmp-${line}-chap.annotations @ Chaperones
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -i "osmoprotectant" >> tmp-${line}-osmo.annotations # Osmorpotectants
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -i "antioxidant" >> tmp-${line}-anox.annotations # Antioxidants
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -i "antifreeze" >> tmp-${line}-anfr.annotations # Antifreeze proteins
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -E -i "acid desaturase|Phytoene desaturase|acyl-CoA desaturase|Sphingolipid Delta4-desaturase|squalene-associated FAD-dependent desaturase" >> tmp-${line}-fatd.annotations # Fatty acid desaturases
	grep -v "^#" ${line}_hgts.annotations | grep -E -v -i "eukaryota|chaperone|S1P" | grep -E -i "heat shock" >> tmp-${line}-hsps.annotations # Hsp
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | grep -E -i "DNA repair|RNA repair|repair protein|repair response|repair endonuclease|repair exonuclease|mismatch repair|photorepair protein|base-excision repair|DNA alkylation repair|damage repair|double-strand break repair|repair of stalled replication forks|repairing DNA|DNA damage lesion repair|base excision repair|double-strand break (DSB) repair|repair of damaged DNA|DNA-repair|UvrABC repair|repair of mismatches in DNA|DNA damage recognition|patch repair" >> tmp-${line}-repa.annotations # Nucleotide repair proteins
	
	head -n5 ${line}_hgts.annotations | tail -n1 > tmp-${line}_CAF.annotations && cat tmp-${line}-*.annotations >> tmp-${line}_CAF.annotations
	cut -f1 tmp-${line}_CAF.annotations | tail -n +2 > tmp-${line}.txt
	paste -d "\t" <(echo -e "#protein_ID\tsilhouette_score\tpotential_donor" && grep -wf tmp-${line}.txt tabs-hgts-hgtector/sw/${line}_hgtector.txt) tmp-${line}_CAF.annotations >> ${line}_CAF.annotations
done < list-bins-hgts.txt

rm tmp*.txt tmp-*CAF.annotations
mkdir tmp-cafs && mv tmp-*-*.annotations tmp-cafs/

while IFS= read -r line; do a=$( grep -w "$line" table-adapt-count-MAGs.txt | awk -F "\t" '{sum+=$3} END {print sum}' ); echo -e "$line\t$a"; done < list_hgts_all.txt

for function in "Cold shock proteins" "Chaperones" "Osmoprotectant-related proteins" "Antioxidants" "Antifreeze proteins" "Fatty acid desaturases" "Heat shock proteins" "Nucleotide repair proteins"
do
	echo -e "\n$function"
	grep -w "$function" table-adapt-count-MAGs.txt | grep -wf list_hgts_all.txt | cut -f3
done

# For lipid transport and metabolism (LTM, COG letter "I"):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_lipid.annotations && grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "I"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CI"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EI"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "IM"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "IQ"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "IT"' >> ${line}_lipid.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "IU"' >> ${line}_lipid.annotations
done < list-bins-hgts.txt

# For energy production and conversion (EPC, COG letter "C"):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_energ.annotations && grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "C"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CE"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CEH"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CG"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CH"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CI"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CJQ"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CK"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CM"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CO"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CP"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CQ"' >> ${line}_energ.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CT"' >> ${line}_energ.annotations
done < list-bins-hgts.txt

# For amino acid transport and metabolism" (ATM, COG letter "E"):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_amino.annotations && grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "E"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CE"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CEH"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EF"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EFGP"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EG"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EGP"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EH"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EHJ"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EI"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EJ"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EK"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EM"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EP"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EQ"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "ET"' >> ${line}_amino.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EU"' >> ${line}_amino.annotations
done < list-bins-hgts.txt

# For carbohydrate transport and metabolism" (CTM, COG letter "G"):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_carbo.annotations && grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "G"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "CG"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EFGP"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EG"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "EGP"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "FG"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GH"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GK"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GM"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GN"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GOQ"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GP"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GPU"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GT"' >> ${line}_carbo.annotations
	grep -v "^#" ${line}_hgts.annotations | grep -v -i "eukaryota" | awk -F "\t" '$7 == "GV"' >> ${line}_carbo.annotations
done < list-bins-hgts.txt

# For metal and antibiotic resistance (MAR):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_resis.annotations && grep -Evi "^#|eukaryota" ${line}_hgts.annotations | grep -iE "resistance|antibiotic|drug|antimicrobial|heavy metal|heavy-metal" | grep -Ev "bacteriophage resistance|PhoQ|Antitoxin|nucleosome-like|oxidative stress|low-pH|cold resistance|resistance to hypoosmotic shock|stress resistance|macrophage|phage-resistance|extreme acid resistance|Toll-Interleukin 1-resistance|Serum resistance|TraT complement|Ultraviolet light resistance|Ultra-violet resistance" >> ${line}_resis.annotations
done < list-bins-hgts.txt

# For machinary of integrative and conjugative elements (ICE):
while IFS= read -r line
do
	head -n5 ${line}_hgts.annotations | tail -n1  > ${line}_ices.annotations && grep -Evi "^#|eukaryota|CRISPR|phage" ${line}_hgts.annotations | grep -iE "integrase|relaxase|relaxosome|conjugative|Type IV secretion|Type-IV|Type IV secretory|T4SS|tyrosine recombinase|ice cds17|TrbI-like|TrbI family|TraG family|unidirectional conjugation|Tra gene|TraG-like|TraB family|TraC-like" >> ${line}_ices.annotations
done < list-bins-hgts.txt


#--- Counting HT genes:

# Count total HT genes annotated:
grep -Evic "^#|eukaryota" *_hgts.annotations | sed "s/_hgts.annotations:/\t/g" # Without Eukarya
grep -vc "^#" *_hgts.annotations | sed "s/_hgts.annotations:/\t/g" # With Eukarya

# Count HT genes annotated of the filtered functional groups and manually make the table "tab-count-hgts.txt" of HGT counts:
echo -e "\tCAF_hgts" & grep -vc "^#" *_CAF.annotations | sed "s/_CAF.annotations:/\t/g" # For cold adaptation
echo -e "\tlipid_hgts" & grep -vc "^#" *_lipid.annotations | sed "s/_lipid.annotations:/\t/g" # For LTM
echo -e "\tenerg_hgts" & grep -vc "^#" *_energ.annotations | sed "s/_energ.annotations:/\t/g" # For EPC
echo -e "\tamino_hgts" & grep -vc "^#" *_amino.annotations | sed "s/_amino.annotations:/\t/g" # For ATM
echo -e "\tcarbo_hgts" & grep -vc "^#" *_carbo.annotations | sed "s/_carbo.annotations:/\t/g" # For CTM
echo -e "\tresis_hgts" & grep -vc "^#" *_resis.annotations | sed "s/_resis.annotations:/\t/g" # For MAR
echo -e "\tices_hgts" & grep -vc "^#" *_ices.annotations | sed "s/_ices.annotations:/\t/g" # For ICE


# Modify table of HGT counts to include extra metadata:
echo -e "mag\tclass\tvariable\tpercentage" > df-count-hgts.txt
while IFS= read -r line
do
	 paste -d "\t" <(yes "$line" | head -n22) <(yes $(grep -w "$line" metadata-bins.tsv | cut -f4) | head -n22) <(cut -f29-35,37,62-68,70-76 tab-count-hgts.txt | head -n1 | tr '\t' '\n' | sed "s/%//g") <(grep -w "$line" tab-count-hgts.txt | cut -f29-35,37,62-68,70-76 | tr '\t' '\n')  >> df-count-hgts.txt
done < list-bins.txt


#--- Assessing the potential donor taxa (from HGTector2 outputs):

# Count total potential donor taxa:
cut -f3 *CAF.annotations | grep -vw "potential_donor" | sort | uniq -c | sort -nr | awk '{print$2"\t"$1}' > tab-count-donors.txt

# Make table of donor taxa by function and MAG:
while IFS= read -r line
do
	yes "${line}" | head -$(cut -f3 ${line}_CAF.annotations | grep -vw "potential_donor" | wc -l) >> tmp-mag.txt
	cut -f3 ${line}_CAF.annotations | grep -vw "potential_donor" >> tmp-tax.txt
	
	for function in csps chap osmo anox anfr fatd hsps repa
	do
		yes "${function}" | head -$(cut -f1 tmp-cafs/tmp-${line}-${function}.annotations | wc -l) >> tmp-fun.txt
	done
done < list-bins-hgts.txt
paste -d "\t" tmp-*.txt > tab-detail-hgts-sw.txt
rm tmp-*.txt

# Search taxonomy of potential donor taxa in Silva and GTDB databases:
while IFS= read -r line
do
	grep -P -m1 "$line" bac120_taxonomy.tsv | cut -f2 | sed s"/;/\t/g" | sed s"/p__//g" | sed s"/c__//g" | sed s"/o__//g" | sed s"/f__//g" | sed s"/g__//g" | sed s"/d__//g" | cut -f1-6 # GTDB database
	grep -wm1 "$line" tax_slv_ssu_138.2.txt | tr ";" "\t" | cut -f1-6 # Silva database
done < list-tmp.txt


#--- HGT detection between recipient MAGs from sponges and donor MAGs from seawater using MetaCHIP:

# Detect putative HT genes between groups (sponges and seawater):
MetaCHIP PI -i input_bins/ -o out_metachip -g hab_list -p hab -x fa -t 5
MetaCHIP BP -o out_metachip -g hab_list -p hab -t 6 -pfr -force

# Detect putative HT genes between taxa (including the taxonomic levels of phylum, class, order, family, and genus):
MetaCHIP PI -i input_bins/ -o out_metachip_tax -taxon tax_list -p tax -r pcofg -x fa -t 6
MetaCHIP BP -o out_metachip_tax -r pcofg -p tax -t 6 -pfr


#--- Functional annotation of putative HT genes (from MetaCHIP outputs) using eggNOG database and eggNOG-mapper:

# For recipient MAGs:
## Between-groups analysis input:
python ../eggnog-mapper-2.1.11/emapper.py -i hab_x_detected_HGTs_recipient_genes.faa -m diamond -o recipient-hgt-hab --go_evidence all --cpu 30 --output_dir ./annotated-hgt-genes
## Between-taxa analysis input:
python ../../eggnog-mapper-2.1.11/emapper.py -i tax_pcofg_detected_HGTs_recipient_genes.faa -m diamond -o recipient-hgt-tax --go_evidence all --cpu 10 --output_dir ../annotated-hgt-genes

# For donors MAGs (these annotations coincide with those for recipient MAGs):
## Between-groups analysis input:
python ../eggnog-mapper-2.1.11/emapper.py -i hab_x_detected_HGTs_donor_genes.faa -m diamond -o donor-hgt-hab --go_evidence all --cpu 30 --output_dir ./annotated-hgt-genes
## Between-taxa analysis input:
python ../../eggnog-mapper-2.1.11/emapper.py -i tax_pcofg_detected_HGTs_donor_genes.faa -m diamond -o donor-hgt-tax --go_evidence all --cpu 10 --output_dir ../annotated-hgt-genes


#--- Extracting annotations with the HGT direction MAGs_seawater --> MAGs_sponge (from MetaCHIP outputs):

# Extract annotations:
head -n5 recipient-hgt.annotations | tail -n1 > recipient-hgt-hab-onlySWtoSP.annotations
head -n5 recipient-hgt-tax.annotations | tail -n1 > recipient-hgt-tax-onlySWtoSP.annotations

while IFS= read -r line; do grep -w "$line" recipient-hgt-hab.annotations ; done < list-recipient-hab-onlySWtoSP.txt >> recipient-hgt-hab-onlySWtoSP.annotations
while IFS= read -r line; do grep -w "$line" recipient-hgt-tax.annotations ; done < list-recipient-tax-onlySWtoSP.txt >> recipient-hgt-tax-onlySWtoSP.annotations

# Extract taxonomy of MAGs with the selected annotations:
while IFS= read -r line; do grep -w "$line" gtdbtk.bac120.summary_sponge.tsv | cut -f7 ; done < list-recipient-hab-onlySWtoSP.txt
while IFS= read -r line; do grep -w "$line" gtdbtk.bac120.summary_sponge.tsv | cut -f7 ; done < list-recipient-tax-onlySWtoSP.txt


#--- Assessing the shared and exclusive Csp and antioxidant orthologs between MAGs from Antarctic sponges and seawater:

# Separate files manually into per-habitat folders:
mkdir Antarctic-sponge Antarctic-seawater

# Search Csp and antioxidant orthologs IDs per each habitat:
## For Csp:
grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "cold" | cut -f2 | sort | uniq > list-csps-sponge.txt # Into "Antarctic-sponge" folder
grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "cold" | cut -f2 | sort | uniq > list-csps-seawater.txt # Into "Antarctic-seawater" folder

## For antioxidants:
grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f2 | sort | uniq > ../list-anox-sponge.txt # Into "Antarctic-sponge" folder
grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f2 | sort | uniq > ../list-anox-seawater.txt # Into "Antarctic-seawater" folder

# Isolate table of Csp and antioxidants:
while IFS= read -r line
do
	grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "cold" | cut -f1,2,8,9 >> table-csp-all.txt
	grep -v "^#" *_func.annotations | grep -v -i "eukaryota" | grep -i "antioxidant" | cut -f1,2,8,9 >> table-anox-all.txt
done < list-bins.txt


#--- Searching HT genes for Csp and antioxidants:

# Check Csp family members and antioxidant genes:
cut -f8 table-csp-all.txt | sort | uniq -c
cut -f8 table-anox-all.txt | sort | uniq -c

# Check HT genes of the different Csp-family members:
grep -i "cold" *hgts.annotations | grep -v -i "eukaryota" | cut -f9 | sort | uniq -c
grep -i "antioxidant" *hgts.annotations | grep -v -i "eukaryota" | cut -f9 | sort | uniq -c

grep -w "cspA" *hgts.annotations | grep -v -i "eukaryota" | cut -f1,9 | sort | uniq -c
grep -w "cspB" *hgts.annotations | grep -v -i "eukaryota" | cut -f1,9 | sort | uniq -c
grep -w "cspC" *hgts.annotations | grep -v -i "eukaryota" | cut -f1,9 | sort | uniq -c
# ... and so on for the rest of the Csp family members and antioxidant genes. HT numbers were added to files "table-*all.txt".


