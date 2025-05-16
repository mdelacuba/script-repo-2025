# ----------------------- #
# PLOTTING AND STATISTICS #
# ----------------------- #

#--- Loading packages ----
library("ggplot2")
library("Biostrings")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("devtools")
library("kableExtra")
library("vegan")
library("ape")
library("tidyverse")
library(wesanderson)
library(scales)
library("data.table")
library(psadd)
library(wrapr)
library(grid)
library(ggpubr)
library(ggtext)
library(png)
library(pairwiseAdonis)
library(RVAideMemoire)
library("picante")
library("coin")
library("factoextra")
library("dendextend")
library("BiodiversityR")
library("maps")
library("clustsig")
library("ggords")
library("pvclust")
library("meowR")
library(jcolors)
library(ggVennDiagram)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(corrplot)
library(RColorBrewer)
library("RColorBrewer")
library(ranacapa)
library(fossil)
library(ggalluvial)
library(grid)
library(plyr)
library(treemapify)


#--- Loading tables and dataframes ----

# NOTES:
## Input file "list-ids-species.txt" contains a list of metagenome names such as "ID_sponge-genus".
## Objects "df.rpkg" and  "df.shan" can be replaced as input data to the boxplots in "Comparing gene numbers of selected functional groups".

# Metadata from metagenomes:
mdata.g2 <- read.table(file = "./metadata.csv", header = T, 
                       sep = "\t", fill = T)

# For COG abundance:
df.cogs <- read.table(file = "./table-cogs-abnorm.txt", 
                      header = TRUE, sep = "\t")

# For gene ordination:
all.rpkg <- read.table("./tmp-all-rpkg.txt", header = T, 
                       sep = "\t", dec = ".")

# For gene content of selected functional groups:
df.all <- read.table(file = "./table-all-func.txt", 
                     header = TRUE, sep = "\t")
df.rpkg <- read.table(file = "./table-all-abnorm.txt", header = TRUE, sep = "\t")
df.shan <- read.table("./table-shannon.txt", header = T, sep = "\t", dec = ".")

# For shared and exclusive orthologs:
## From metagenomes:
orto.anta <- read.table(file = "./Antarctic/list-orthologs-anta.txt", 
                        header = FALSE, sep = "\t")
orto.temp <- read.table(file = "./Temperate/list-orthologs-temp.txt", 
                        header = FALSE, sep = "\t")
orto.trop <- read.table(file = "./Tropical/list-orthologs-trop.txt", 
                        header = FALSE, sep = "\t")
csp.anta <- read.table(file = "./Antarctic/list-csps-anta.txt", 
                       header = FALSE, sep = "\t")
csp.temp <- read.table(file = "./Temperate/list-csps-temp.txt", 
                       header = FALSE, sep = "\t")
csp.trop <- read.table(file = "./Tropical/list-csps-trop.txt", 
                       header = FALSE, sep = "\t")
## From MAGs:
csp.spon <- read.delim(file = "./list-csps-sponge.txt", 
                       header = FALSE, sep = "\t")
csp.seaw <- read.delim(file = "./list-csps-seawater.txt", 
                       header = FALSE, sep = "\t")
anox.spon <- read.delim(file = "./list-anox-sponge.txt", 
                        header = FALSE, sep = "\t")
anox.seaw <- read.delim(file = "./list-anox-seawater.txt", 
                        header = FALSE, sep = "\t")

# For presence/absence of gene orthologs for cold adaptation:
all.subs.cold <- read.table("./subsampling/all-subs-coldfuncs.txt", header = T, 
                            sep = "\t", dec = ".")
annot.subs <- read.table(file = "./subsampling/annot-subs.txt", header = FALSE, 
                         sep = "\t", fill = TRUE)
ids.new <- read.table(file = "./subsampling/list-ids-genera.txt", header = F, 
                      sep = "\t", fill = TRUE)

# Metadata from MAGs:
tab.bins <- read.table(file = "./metadata-bins.tsv", 
                       header = TRUE, sep = "\t")

# For correlation between gene content and HGT genes:
tab.corr <- read.table(file = "./tab-corr.txt", 
                       header = TRUE, sep = "\t")

# For HGT numbers of selected functional groups from HGTector 2 outputs:
tab.hgt.h <- read.table(file = "./tab-count-hgts.txt", 
                        header = TRUE, sep = "\t")
df.hgt.h <- read.table(file = "./df-count-hgts.txt", 
                        header = TRUE, sep = "\t")
                        
# For HGT recipients and donors (HGTector 2 outputs):
hgt.taxa <- read.table(file = "./plotting_donors_recipients.csv", 
                       header = TRUE, sep = "\t")
hgt.taxa.na <- read.table(file = "./count_donors_NA.csv", 
                          header = TRUE, sep = "\t")
                          
# For HGT detection from MetaCHIP outputs:
df.hgt.m <- read.table(file = "./recipient-hgt-onlySWtoSP.annotations", 
                       header = TRUE, sep = "\t")                        

# For distribution and proportions of Csp and antioxidants:
df.csps.all <- read.table(file = "./table-csp-all.txt", 
                          header = TRUE, sep = "\t", na.strings = "-")
df.anox.all <- read.table(file = "./table-anox-all.txt", 
                          header = TRUE, sep = "\t", na.strings = "-")                       


#--- Loading color palettes ----

colors.env.s <- c("deepskyblue2", "#2166ac", "#DD8B71", "#b2182b")

colors <- c("coral3", "darkblue", "darkorchid4", "darkgoldenrod3", 
            "seagreen4", "tomato4", "deepskyblue2", "turquoise", "lightpink", 
            "darkcyan", "orchid1",  "olivedrab4", "red1", "skyblue3", "maroon4", 
            "darkorange", "gold2", "khaki", "plum4", "maroon2", "lightgreen", 
            "tan", "grey50", "mediumpurple", "sienna4", "aquamarine", "limegreen", 
            "yellow", "thistle", "gray81")


colors.env <- c('Antarctic sponge' = "#2166ac", 
                'Temperate sponge' = "#DD8B71", 'Tropical sponge'= "#b2182b",
                'Antarctic seawater' = "deepskyblue2")

colors.env2 <- list(Habitat = 
                     c('Antarctic sponge' = "#2166ac",
                       'Temperate sponge' = "#DD8B71", 
                       'Tropical sponge'= "#b2182b", 
                       'Antarctic seawater' = "deepskyblue2"))

colors.func <- list(Function= c('Antifreeze'= "darkcyan",
                                'Antioxidant'= "darkorchid4",
                                'Chaperone related'= "darkgoldenrod3",
                                'Cold-shock related'= "darkblue",
                                'Fatty acid desaturase'= "seagreen4",
                                'Heat-shock related'= "maroon4",
                                'Osmoprotectant related'= "darkorange",
                                'Nucleotide repair'= "skyblue3"))

color.class <- c('Nitrospinia' = "#674ea7",
                 'Alphaproteobacteria' = "#d9b7ce",
                 'Bacteroidia' = "#eda137",
                 'Bdellovibrionia' = "#ef7605",
                 'Gammaproteobacteria' = "#a6558b",
                 'Actinomycetia'= "#49BEAA",
                 'Spirochaetia' = "#456990",
                 'Acidimicrobiia' = "#6DAEDB",
                 'Nitrososphaeria' = "#b8b893")

colors.func2 <- c('Antioxidant'= "darkorchid4",
                  'Chaperone'= "darkgoldenrod3",
                  'Cold-shock protein'= "darkblue",
                  'Fatty acid desaturase'= "seagreen4",
                  'Heat-shock protein'= "maroon4",
                  'Osmoprotectant'= "darkorange",
                  'Nucleotide repair'= "skyblue3")

colors.2 <- c("darkorange", "deepskyblue2")

colors.3 <- c("deepskyblue2", "#ffc100", "#ff9a00", "#ff4d00", "#ff0000")


#--- Mapping of sampling sites (Figure 1A) ----

# Include world data:
world <- map_data("world") %>% filter(!long > 180)

# Include Marine Ecoregions Of the World (MEOW) data:
data(regions)
data(provinces)
data(realms)
realms.df
provinces.df
regions.df

# Plot world map including geographic coordinates of habitats/environments:
plot.map <- ggplot() +
  geom_map(data = world, map = world, aes(map_id = region),
           color = "white", fill = "lightgray", linewidth = 0.1) +
  expand_limits(x = world$long, y = world$lat) +
  geom_point(data = mdata.g2, aes(x = longitude, y = latitude, color = Environment),
             size = 2, show.legend = TRUE, shape = 0) +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_light() +
  scale_x_continuous(breaks = seq(-150, 150, 50)) +
  scale_y_continuous(breaks = seq(-80, 80, 20)) +
  scale_color_manual(values = colors.env.s, labels = c('Antarctic seawater (n=4)',
                                                       'Antarctic sponge (n=17)',
                                                       'Temperate sponge (n=13)',
                                                       'Tropical sponge (n=15)')) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 13))

# Export plot as image:
png("../graphics/map-sampling.png", units = "in", width = 8, height = 3.5, res = 300)
plot.map
dev.off()


#--- Assessing COG categories abundance (Figure 1B) ----

# Reorder according to abundance:
df.cogs$description <- reorder(df.cogs$description, df.cogs$count)
df.cogs$description <- factor(df.cogs$description, levels = rev(levels(df.cogs$description)))

# Plot barplot of COGs:
plot.cogs <- ggplot(data = df.cogs, aes(x = ID, y = count, fill = description)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  theme_classic() +
  facet_wrap(~Environment, nrow = 1, scale = "free_x") +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 15), legend.text = element_text(size = 14),
        legend.key.size = unit(0.7, "cm"),
        strip.text.x = element_text(size = 17)) +
  #theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = colors) +
  ylab("RPKG") +
  labs(x = NULL, fill = "COG category") +
  theme(legend.position = "right")

# Export plot as image:
png("../graphics/cogs-abnor.png", units = "in", width = 16, height = 8, res = 300)
plot.cogs.tmp
dev.off()


#--- Ordination of gene orthologs (Figure 1C) ----

# Create a matrix of abundances:
orthotab <- create.matrix(all.rpkg, tax.name = "Seed_ortholog",
                          locality = "Sample",
                          abund.col = "RPKG",
                          abund = TRUE)

# Convert to dataframe and transpose:
orthodf <- as.data.frame(orthotab)
orthodf <- t(orthodf)

# Make ordination:
ord.ortho <- metaMDS(orthodf, distance = 'bray', k=3, try = 10,
                     trymax = 100, autotransform = F, pc = TRUE)

# First plot needed:
plot.ord <- ordiplot(ord.ortho, choices=c(1,2,3))

# Attach and include metadata:
attach(mdata.g2)
sites.long <- sites.long(plot.ord, env.data = mdata.g2)

# Plot NMDS:
plot.nmds <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data = sites.long.all, 
             aes(x = axis1, y = axis2, colour = Environment), 
             shape = 15, size = 4) +
  scale_color_manual(values = colors.env) +
  theme_classic() +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_blank()) +
  annotation_custom(textGrob(label = paste0("stress = ", 
                                            round(ord.ortho.all$stress, 3)), 
                             x = 0.75, y = 0.95, hjust = 0))

# Export plot as image:
png("../graphics/nmds-ortho.png", units = "in", width = 7.5, height = 5, res = 300)
plot.nmds
dev.off()


# Perform PERMANOVA:
set.seed(1000)
ortho.dist <- vegdist(orthodf, method = "bray")

ortho.perm <- adonis2(ortho.dist ~ Habitat,
                      data = mdata.g2,
                      permutations = 999)

set.seed(1000)
ortho.perw <- pairwise.adonis(ortho.dist,
                              factors = mdata.g2$Habitat,
                              p.adjust.m = "bonferroni",
                              perm = 999)

set.seed(1000)
ortho.disp <- betadisper(ortho.dist, mdata.g2$Habitat)
ortho.homog <- permutest(ortho.disp, permutations = 999)

# Perform ANOSIM:
set.seed(1000)
ortho.anosim <- anosim(x = ortho.dist, grouping = mdata.g2$Habitat,
                       permutations = 999)


#--- Comparing gene content of selected functional groups (Figure 2 and Supplementary Figure S3) ----

# For cold adaptation functions:
## Cold shock proteins (Csp): 
plot.csp <- ggplot(data = df.all, aes(x = Environment, y = csp_per, fill = Environment)) +
  geom_boxplot(alpha = 0.9) + #AQUÍ , width = 0.5
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) + #width = 0.2
  labs(x = NULL, y = "% Cold-shock related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "right", text = element_text(size = 21),
        axis.text.x = element_blank(),#element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Temperate sponge"),
                       c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))
## Chaperones:
plot.chap <- ggplot(data = df.all, aes(x = Environment, y = chaperone_per, fill = Environment)) +
  geom_boxplot(alpha = 0.9) +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Chaperone-related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),#element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Temperate sponge"),
                       c("Antarctic sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))
## Osmoprotectants:
plot.osmo <- ggplot(data = df.all, aes(x = Environment, y = osmo_per, fill = Environment)) +
  geom_boxplot(alpha = 0.9) +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Osmoprotectant-related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),#element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       c("Antarctic sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Antarctic seawater"),
                       #c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))
## Antioxidants:
plot.aox <- ggplot(data = df.all, aes(x = Environment, y = antiox_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Antioxidant-related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge")))
                        #c("Antarctic seawater", "Tropical sponge")))
# Heat shock proteins (Hsp):
plot.hsp <- ggplot(data = df.all, aes(x = Environment, y = hsp_per, fill = Environment)) +
  geom_boxplot(alpha = 0.9) +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Heat-shock related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),#element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge")))
                        #c("Antarctic seawater", "Tropical sponge")))
## Fatty acid desaturases:
plot.fad <- ggplot(data = df.all, aes(x = Environment, y = fatdesat_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Fatty acid desaturases", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23))
# stat_compare_means(method = "t.test",
#                    label = c("p.signif"), label.x.npc = "center",
#                    label.y.npc = "top", hide.ns = F, size = 6.5,
#                    comparisons = list(c("Temperate sponge", "Tropical sponge"),
#                                       c("Antarctic sponge", "Temperate sponge"),
#                                       c("Antarctic sponge", "Tropical sponge"),
#                                       c("Antarctic sponge", "Antarctic seawater"),
#                                       c("Antarctic seawater", "Temperate sponge"),
#                                       c("Antarctic seawater", "Tropical sponge")))
## Nucleotide repair proteins:
plot.rep <- ggplot(data = df.all, aes(x = Environment, y = repair_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% Nucleotide-repair related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23))
# stat_compare_means(method = "t.test",
#                    label = c("p.signif"), label.x.npc = "center",
#                    label.y.npc = "top", hide.ns = F, size = 6.5,
#                    comparisons = list(c("Temperate sponge", "Tropical sponge"),
#                                       c("Antarctic sponge", "Temperate sponge"),
#                                       c("Antarctic sponge", "Tropical sponge"),
#                                       c("Antarctic sponge", "Antarctic seawater"),
#                                       c("Antarctic seawater", "Temperate sponge"),
#                                       c("Antarctic seawater", "Tropical sponge")))

# For Lipid transport and metabolism (LTM):
plot.lip <- ggplot(data = df.all, aes(x = Environment, y = lipid_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% LTM", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))

# For energy production and conversion (EPC):
plot.ene <- ggplot(data = df.all, aes(x = Environment, y = energ_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% EPC", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23))# +
# stat_compare_means(method = "t.test",
#                    label = c("p.signif"), label.x.npc = "center",
#                    label.y.npc = "top", hide.ns = F, size = 6.5,
#                    comparisons = list(c("Temperate sponge", "Tropical sponge"),
#                      c("Antarctic sponge", "Temperate sponge"),
#                      c("Antarctic sponge", "Tropical sponge"),
#                      c("Antarctic sponge", "Antarctic seawater"),
#                      c("Antarctic seawater", "Temperate sponge"),
#                      c("Antarctic seawater", "Tropical sponge")))

# For amino acid transport and metabolism (ATM):
plot.ami <- ggplot(data = df.all, aes(x = Environment, y = amino_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% ATM", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))

# For carbon transport and metabolism (CTM):
plot.car <- ggplot(data = df.all, aes(x = Environment, y = carbo_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% CTM", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge")))
                        #c("Antarctic seawater", "Tropical sponge")))

# For metal and antibiotic resistance (MAR):
plot.res <- ggplot(data = df.all, aes(x = Environment, y = resis_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% MAR", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater")))
                        #c("Antarctic seawater", "Temperate sponge"), 
                        #c("Antarctic seawater", "Tropical sponge")))

# For machinary of integrative and conjugative elements (ICE):
plot.ice <- ggplot(data = df.all, aes(x = Environment, y = ices_per, fill = Environment)) +
  geom_boxplot() +
  scale_fill_manual(values = colors.env) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% ICE-related", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", text = element_text(size = 21),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 23)) +
  stat_compare_means(method = "t.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = F, size = 6.5,
                     comparisons = list(#c("Temperate sponge", "Tropical sponge"),
                       #c("Antarctic sponge", "Temperate sponge"),
                       #c("Antarctic sponge", "Tropical sponge"),
                       c("Antarctic sponge", "Antarctic seawater"),
                       c("Antarctic seawater", "Temperate sponge"),
                       c("Antarctic seawater", "Tropical sponge")))

# Export plots as images:
png("../graphics/box-csp.png", units = "in", width = 8, height = 6, res = 300)
plot.csp
dev.off()

png("../graphics/box-chaperone.png", units = "in", width = 8, height = 6, res = 300)
plot.chap
dev.off()

png("../graphics/box-osmop.png", units = "in", width = 8, height = 6, res = 300)
plot.osmo
dev.off()

png("../graphics/box-antiox.png", units = "in", width = 8, height = 6, res = 300)
plot.aox
dev.off()

png("../graphics/box-hsp.png", units = "in", width = 8, height = 6, res = 300)
plot.hsp
dev.off()

png("../graphics/box-fatdesat.png", units = "in", width = 8, height = 6, res = 300)
plot.fad
dev.off()

png("../graphics/box-repair.png", units = "in", width = 8, height = 6, res = 300)
plot.rep
dev.off()

png("../graphics/box-lipid.png", units = "in", width = 8, height = 6, res = 300)
plot.lip
dev.off()

png("../graphics/box-energ.png", units = "in", width = 8, height = 6, res = 300)
plot.ene
dev.off()

png("../graphics/box-amino.png", units = "in", width = 8, height = 6, res = 300)
plot.ami
dev.off()

png("../graphics/box-carbo.png", units = "in", width = 8, height = 6, res = 300)
plot.car
dev.off()

png("../graphics/box-resis.png", units = "in", width = 8, height = 6, res = 300)
plot.res
dev.off()

png("../graphics/box-ices.png", units = "in", width = 8, height = 6, res = 300)
plot.ice
dev.off()


#--- Shared and exclusive orthologs (Supplementary Figures S1 and S4 and Figure 6) ----

# Create lists from metagenomes:
## For all genes:
list.venn <- list('Antarctic sponge' = as.vector(as.matrix(orto.anta)), 
                  'Temperate sponge' = as.vector(as.matrix(orto.temp)),
                  'Tropical sponge' = as.vector(as.matrix(orto.trop)))
## For genes encoding cold shock protein (Csp):
list.venn.csp <- list('Antarctic sponge' = as.vector(as.matrix(csp.anta)), 
                      'Temperate sponge' = as.vector(as.matrix(csp.temp)),
                      'Tropical sponge' = as.vector(as.matrix(csp.trop)))

# Plot Venn diagrams from metagenomes:
## For all genes:
plot.ortho <- ggVennDiagram(list.venn, 
                            category.names = c("Antarctic", "Temperate", "Tropical"),
                            color = "black", lwd = 0.8, lty = 1,
                            label = "percent", label_alpha = 0, label_size = 5,
                            set_size = 5) + 
  scale_fill_gradient(low = "#fcf2d6", high = "#d8ae2d") +
  scale_color_manual(values = colors.env.s) +
  labs(fill = "# orthologs") +
  theme(text = element_text(size = 15)) +
  scale_x_continuous(expand = expansion(mult = .1))
## For genes encoding cold shock protein (Csp):
plot.ortho.csp <- ggVennDiagram(list.venn.csp,
                                category.names = c("Antarctic", "Temperate", "Tropical"),
                                color = "black", lwd = 0.8, lty = 1,
                                label = "percent", label_alpha = 0, label_size = 5,
                                set_size = 5) + 
  scale_fill_gradient(low = "#e0f1f3", high = "#64b9c6") +
  scale_color_manual(values = colors.env.s) +
  labs(fill = "# Csp orthologs") +
  theme(text = element_text(size = 15)) +
  scale_x_continuous(expand = expansion(mult = .1))

# Create lists from MAGs:
## For genes encoding Csp:
list.venn.csp <- list('Antarctic sponge' = as.vector(as.matrix(csp.spon)), 
                      'Antarctic seawater' = as.vector(as.matrix(csp.seaw)))
## For genes encoding antioxidants:
list.venn.anox <- list('Antarctic sponge' = as.vector(as.matrix(anox.spon)), 
                       'Antarctic seawater' = as.vector(as.matrix(anox.seaw)))

# Plot Venn diagrams from MAGs:
## For genes encoding Csp:
plot.ortho.csp <- ggVennDiagram(list.venn.csp,
                                # category.names = c("Antarctic sponge", 
                                #                    "Antarctic seawater"),
                                color = "black", lwd = 0.8, lty = 1,
                                label = "both", label_alpha = 0, label_size = 6,
                                set_size = 5.5) +
  scale_fill_gradient(low = "#e0f1f3", high = "#64b9c6") +
  scale_color_manual(values = colors.2) +
  labs(fill = "# Orthologs") +
  theme(text = element_text(size = 18)) +
  coord_flip() +
  scale_x_continuous(expand = expansion(mult = .1))
## For genes encoding antioxidants:
plot.ortho.anox <- ggVennDiagram(list.venn.anox,
                                 # category.names = c("Antarctic sponge", 
                                 #                    "Antarctic seawater"),
                                 color = "black", lwd = 0.8, lty = 1,
                                 label = "both", label_alpha = 0, label_size = 6,
                                 set_size = 5.5) +
  scale_fill_gradient(low = "#e0f1f3", high = "#64b9c6") +
  scale_color_manual(values = colors.tmp) +
  labs(fill = "# Orthologs") +
  theme(text = element_text(size = 18)) +
  coord_flip() +
  scale_x_continuous(expand = expansion(mult = .1))

# Export plots as images:
png("../graphics/venn-ortho-metag.png", units = "in", width = 7, height = 5, res = 300)
plot.ortho
dev.off()

png("../graphics/venn-csp-metag.png", units = "in", width = 7, height = 5, res = 300)
plot.ortho.csp
dev.off()

png("./graphics/venn-csp-mag.png", units = "in", width = 7.5, height = 5, res = 300)
plot.ortho.csp
dev.off()

png("./graphics/venn-anox-mag.png", units = "in", width = 7.5, height = 5, res = 300)
plot.ortho.anox
dev.off()


#--- Presence/absence of gene orthologs for cold adaptation (Supplementary Figure S2) ----

# Create a matrix of presence/absence:
matrix.subs <- create.matrix(all.subs.cold, tax.name = "Seed_ortholog",
                             locality = "Sample",
                             abund.col = "Freq",
                             abund = FALSE)

# Prepare annotation format:
## For metadata:
rownames(mdata.g2) <- mdata.g2$ID_dataset
annot.data <- data.frame(mdata.g2$Environment)
rownames(annot.data) <- rownames(mdata.g2)
colnames(annot.data) <- "Habitat"
## For functions:
rownames(annot.subs) <- annot.subs$V1
annot.gene <- data.frame(annot.subs$V2)
rownames(annot.gene) <- rownames(annot.subs)
colnames(annot.gene) <- "Function"

# Order orthologs in matrix according to annotation order:
matrix.subs.sort <- matrix.subs[,rownames(annot.data)]
matrix.subs.sort <- matrix.subs[rownames(annot.gene), ]

# Change names of metagenomes:
rownames(annot.data) <- ids.new$V1
colnames(matrix.subs.sort) <- ids.new$V1

# Plot heatmap of presence/absence:
set.seed(1000)
plot.heatmap.large <- Heatmap(matrix.subs.sort, show_row_names = F,
                              cluster_rows = F, cluster_columns = T,
                              col = c("whitesmoke", "#000"),
                              column_names_side = "top",
                              column_dend_side = "top",
                              column_dend_height = unit(6.5, "cm"),
                              # row_dend_width = unit(4, "cm"),
                              # clustering_distance_rows = "euclidean",
                              # clustering_method_rows = "ward.D",
                              clustering_distance_columns = "euclidean",
                              clustering_method_columns = "ward.D",
                              # column_dend_reorder = TRUE,
                              # row_dend_reorder = TRUE,
                              # row_split = matrix.cold.sort[, 1],
                              # column_split = 3,
                              # row_split = 3,
                              # column_km = 3,
                              # row_km = 3,
                              # cluster_column_slices = F,
                              top_annotation = HeatmapAnnotation(df = annot.data,
                                                                 col = colors.env2),
                              left_annotation = rowAnnotation(df = annot.gene,
                                                              col = colors.func),
                              use_raster = TRUE, raster_resize_mat = max)

# Export figure:
png("../graphics/heatmap-hier-all.png", units = "in", width = 10, height = 20, res = 300)
plot.heatmap.large
dev.off()


#--- Abundance of MAG classes (Figure 3) ----

# Reorder habitat and class data:
tab.bins$habitat <- factor(tab.bins$habitat, levels = c("Antarctic sponge", 
                                                        "Antarctic seawater"))
tab.bins$class <- reorder(tab.bins$class, tab.bins$abundance)
tab.bins$class <- factor(tab.bins$class, levels = rev(levels(tab.bins$class)))

# Plot barplot of classes:
class.plot <- ggplot(tab.bins, aes(x = habitat, y = abundance, fill = class)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  labs(x = NULL, y = "relative abundance") +
  theme(legend.position = "none",
        text = element_text(size = 16), legend.text = element_text(size = 16),
        legend.key.size = unit(0.7, "cm"),
        axis.text.x = element_text(size = 16)) +
  scale_fill_manual(values = color.class)

# Export plot as image:
png("./graphics/classes-mags.png", units = "in", width = 5, height = 5, res = 300)
class.plot
dev.off()


#--- Correlation between % gene content and % HGT genes (Supplementary Figure S5) ----

# Fix names of the functional groups:
tab.corr$genes <- gsub("caf", "Cold adaptation", tab.corr$genes)
tab.corr$genes <- gsub("amino", "ATM", tab.corr$genes)
tab.corr$genes <- gsub("energ", "EPC", tab.corr$genes)
tab.corr$genes <- gsub("lipid", "LTM", tab.corr$genes)
tab.corr$genes <- gsub("carbo", "CTM", tab.corr$genes)
tab.corr$genes <- gsub("resis", "MAR", tab.corr$genes)
tab.corr$genes <- gsub("ices", "ICE", tab.corr$genes)
tab.corr$genes <- factor(tab.corr$genes, levels = c("ATM", "EPC", 
                                                    "LTM", "CTM", 
                                                    "Cold adaptation",
                                                    "MAR", "ICE"))
# Separate tables:
# tab.corr.sp <- tab.corr[df2.hgt.h$habitat == "Sponge",]
# tab.corr.sw <- tab.corr[df2.hgt.h$habitat == "Seawater",]

# Linear regressions:
ggplot(tab.corr, aes(genes_per, hgts_per)) +
  geom_point() +
  geom_smooth(method='lm')

# Plot correlations:
## Measure based on dividing by the total HGTs number:
plot.corr <- ggscatter(tab.corr, x = "genes_per", y = "hgts_per", 
                       color = "genes", palette = "Dark2") +#, add = "reg.line"
  stat_cor(label.x = 3, label.y = 30) +
  # stat_regline_equation(label.x = 3, label.y = 27) +
  labs(x = "% genes 'X'", y = "% genes 'X' HT per total genes HT") +
  theme(legend.title = element_blank(), legend.position = "none") +
## Measure based on dividing by the gene content: 
  ggscatter(tab.corr, x = "genes_per", y = "hgts.genes_per", 
            color = "genes", palette = "Dark2") + #add = "reg.line",
  stat_cor(label.x = 3, label.y = 100) +
  # stat_regline_equation(label.x = 3, label.y = 90) +
  labs(x = "% genes 'X'", y = "% genes 'X' HT per total genes 'X'") +
  theme(legend.title = element_blank(), legend.position = "right")

# Export plot as image:
png("./graphics/corr-spsw.png", units = "in", width = 11.5, height = 5, res = 300)
plot.corr
dev.off()


#--- Comparing HGT numbers of selected functional groups  (Figure 4) ----

# Customize tables:
tab.hgt.h$habitat <- factor(tab.hgt.h$habitat, levels = c("Sponge", "Seawater"))
tab2.hgt.h <- rbind(tab.hgt.h[tab.hgt.h$class == "Gammaproteobacteria",],
                    tab.hgt.h[tab.hgt.h$class == "Alphaproteobacteria",],
                    tab.hgt.h[tab.hgt.h$class == "Bacteroidia",])
tab2.hgt.h$habitat <- factor(tab2.hgt.h$habitat, levels = c("Sponge", "Seawater"))

# Customize dataframes:
df.hgt.h$variable <- gsub("_hgts/genes", "", df.hgt.h$variable)
df.hgt.h$variable <- gsub("caf", "Cold adaptation", df.hgt.h$variable)
df.hgt.h$variable <- gsub("amino", "ATM", df.hgt.h$variable)
df.hgt.h$variable <- gsub("energ", "EPC", df.hgt.h$variable)
df.hgt.h$variable <- gsub("lipid", "LTM", df.hgt.h$variable)
df.hgt.h$variable <- gsub("carbo", "CTM", df.hgt.h$variable)
df.hgt.h$variable <- gsub("resis", "MAR", df.hgt.h$variable)
df.hgt.h$variable <- gsub("ices", "ICE", df.hgt.h$variable)
df.hgt.h$habitat <- factor(df.hgt.h$habitat, levels = c("Sponge", "Seawater"))
df.hgt.sp <- df.hgt.h[df.hgt.h$habitat == "Sponge",]
df.hgt.h <- rbind(df.hgt.h[df.hgt.h$class == "Gammaproteobacteria",],
                   df.hgt.h[df.hgt.h$class == "Alphaproteobacteria",],
                   df.hgt.h[df.hgt.h$class == "Bacteroidia",])

# Plot boxplot of cold adaptation functions:
caf.plot <- ggplot(data = tab2.hgt.h, 
                   aes(x = habitat, 
                       y = X.caf, color = habitat)) +
  # scale_color_manual(values =  "grey22") +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) + #width = 0.2
  labs(x = NULL, y = "% genes for cold adaptation", fill = NULL, color = NULL) +
  theme_classic() +
  #with facet_wrap: nrow = 2, strip.position = "right", y cambiar ~ al inicio
  theme(legend.position = "none", axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) +
  # scale_y_reverse() +
  # coord_flip() +
  geom_boxplot(alpha = 0.4) + #width = 0.5
  scale_color_manual(values = colors.2) +
  stat_compare_means(method = "wilcox.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = FALSE, size = 6)

# Plot boxplot of total HGTs:
tot.hgt.plot <- ggplot(data = tab2.hgt.h, 
                       aes(x = habitat, 
                           y = X.total_hgts, color = habitat)) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% horizontally-transferred genes", fill = NULL, color = NULL) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 16),
        axis.title = element_text(size = 18)) +
  geom_boxplot(alpha = 0.4) +
  scale_color_manual(values = colors.2) +
  stat_compare_means(method = "wilcox.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = FALSE, size = 6)

# Plot boxplots of HGTs of functional groups:
spsw.plot <- ggplot(data = df.hgt.h, 
                    aes(x = habitat, 
                        y = percentage, color = habitat)) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% horizontally transferred genes", fill = NULL, color = NULL) +
  theme_classic() +
  facet_wrap(~factor(variable, c("EPC", "LTM", "ATM","Cold adaptation",
                                 "CTM", "MAR","ICE")), 
             nrow = 1, scale = "free_x") +
  geom_boxplot(alpha = 0.4) +
  theme(legend.position = "bottom", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), text = element_text(size = 16),
        axis.title = element_text(size = 18), strip.text = element_text(size = 13.5)) +
  scale_color_manual(values = colors.2) +
  stat_compare_means(method = "wilcox.test",
                     label = c("p.signif"), label.x.npc = "center",
                     label.y.npc = "top", hide.ns = T, size = 6)

# Plot boxplots of HGTs of functional groups only for MAGs from Antarctic sponges:
sp.hgts.plot <- ggplot(data = df.hgt.sp, 
                       aes(x = fct_reorder(variable, percentage, 
                                           .fun = median, .desc = T), 
                           y = percentage)) +
  geom_jitter(alpha = 0.3, width = 0.1, size = 2) +
  labs(x = NULL, y = "% genes horizontally transferred", fill = NULL, color = NULL) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), strip.text = element_text(size = 13.5)) +
  scale_color_manual(values = colors.2) +
  stat_summary(fun = "mean", geom = "point", 
               shape = 17, size = 3,  color = "blue", fill = "blue") +
  geom_boxplot(alpha = 0.4)

# Export plots as images:
png("./graphics/boxplot-spsw-caf.png", units = "in", width = 6, height = 4, res = 300)
caf.plot
dev.off()

png("./graphics/boxplot-spsw-tothgt.png", units = "in", width = 6, height = 4.1, res = 300)
tot.hgt.plot
dev.off()

png("./graphics/boxplot-spsw-hgts.png", units = "in", width = 13, height = 5.5, res = 300)
spsw.plot
dev.off()

png("./graphics/boxplot-sp-hgts.png", units = "in", width = 8, height = 5, res = 300)
sp.hgts.plot
dev.off()


#--- HGT recipients MAGs and reference donors for cold adaptation functions (Figure 5) ----

# Improve HGT distance info:
hgt.taxa$hgt_type <- str_replace(hgt.taxa$hgt_type, 
                                 "cross-domain", "Cross-Domain (7 donors)")
hgt.taxa$hgt_type <- str_replace(hgt.taxa$hgt_type, 
                                 "cross-phylum", "Cross-Phylum (30 donors)")
hgt.taxa$hgt_type <- str_replace(hgt.taxa$hgt_type, 
                                 "cross-class", "Cross-Class (3 donors)")
hgt.taxa$hgt_type <- str_replace(hgt.taxa$hgt_type, 
                                 "cross-order", "Cross-Order (6 donors)")
hgt.taxa$hgt_type <- str_replace(hgt.taxa$hgt_type, 
                                 "cross-family", "Cross-Family (2 donors)")

# Reorder facets (HGT distances):
hgt.taxa$hgt_type <- factor(hgt.taxa$hgt_type, levels=c("Cross-Domain (7 donors)",
                                                        "Cross-Phylum (30 donors)",
                                                        "Cross-Class (3 donors)",
                                                        "Cross-Order (6 donors)",
                                                        "Cross-Family (2 donors)"))

# Reorder cold adaptation functions:
hgt.taxa$function. <- factor(hgt.taxa$function.,
                             levels=c("Nucleotide repair",
                                      "Chaperone",
                                      "Heat-shock protein",
                                      "Antioxidant",
                                      "Osmoprotectant",
                                      "Cold-shock protein",
                                      "Fatty acid desaturase"))

# Plot dotplot with HGT distances for cold adaptation functions:
donors.plot <- ggplot(hgt.taxa, aes(x = hgt_type, 
                                    y = fct_reorder(recipient_name, order))) +
  geom_point(aes(fill = function.,), color = "black",
             shape = 21, stroke = 0.1, size = 5,
             position = position_dodge2(width = 1)) +
  facet_wrap(~hgt_type,
             nrow = 1,
             scales = "free_x",
             strip.position = "top") +
  labs(x = "Genes HT",
       y = "Recipient MAGs", fill = "Function for cold adaptation") +
  
  scale_fill_manual(values = colors.func2) +
  scale_x_discrete(breaks = NULL) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 150)),
        axis.title.x = element_text(margin = margin(t = 30, b = 30)),
        legend.position = "bottom",
        axis.text.y = element_text(face = 'italic'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"))

# Reorder functions as above for genes with undetermined HGT distances:
hgt.taxa.na$function. <- factor(hgt.taxa.na$function.,
                                levels = rev(c("Nucleotide repair",
                                               "Chaperone",
                                               "Heat-shock protein",
                                               "Antioxidant",
                                               "Osmoprotectant",
                                               "Cold-shock protein",
                                               "Fatty acid desaturase")))

# Plot barplot of genes with undetermined HGT distances:
donors.na.plot <- ggplot(hgt.taxa.na, aes(x = function., y = na_count, fill = function.)) +
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = colors.func2) +
  labs(x = NULL,
       y = "# HT genes with uncertain donor/distance", 
       fill = "Function for cold adaptation") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.1,
                                          linetype = 2),
        panel.grid.minor.x = element_line(size = 0.1,
                                          linetype = 2)) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140))

# Export plots as images:
png("./graphics/hgt-donors.png", units = "in", width = 13, height = 6, res = 300)
donors.plot
dev.off()

png("./graphics/hgt-donors-na.png", units = "in", width = 4, height = 2, res = 300)
donors.na.plot
dev.off()


#--- HGT detection between MAGs from MetaCHIP outputs (Figure S6) ----

# Fix details:
df.hgt.m$sponge_name <- str_replace(df.hgt.m$sponge_name, "\\*", "")

# Plot dotplot of HGTs between MAGs:
mc.hgt.plot <- ggplot(df.hgt.m, 
                      aes(x = sponge_tax_short, 
                          y = reorder(KEGG_ko_name, 
                                      desc(COG_category_name)))) + 
  geom_point(aes(fill = COG_category_name), 
             shape = 21, size = 5) +
  labs( x= NULL, y = NULL, fill = "COG Category") + 
  theme_bw() +
  facet_wrap(~sponge_name2,
             nrow = 1,
             scales = "free_x",
             strip.position = "top") +
  theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 6.8, face = "bold.italic")) +  
  scale_fill_brewer(palette = "Set3")

# Export Figure:
png("./graphics/hgt-metachip.png", units = "in", width = 19, height = 10, res = 300)
mc.hgt.plot
dev.off()


#--- Distribution, proportions, and acquisition of Csp and antioxidants (Figure 6) ----

# Omit NAs
df.csps.all <- na.omit(df.csps.all)
df.anox.all <- na.omit(df.anox.all)

# Include genes without standardized nomenclature:
# df.csps.all[is.na(df.csps.all)] <- "Non-ID genes"
# df.anox.all[is.na(df.anox.all)] <- "Non-ID genes"

# Add fequency to table:
df.csps.all$Freq <- rep(x = 1, 96) #175 with unassigned
df.anox.all$Freq <- rep(x = 1, 41) #186 with unassigned

# Change mistaken class name:
df.csps.all$taxa <- gsub("Actinomycetia", "Actinomycetes", df.csps.all$taxa)
df.anox.all$taxa <- gsub("Actinomycetia", "Actinomycetes", df.anox.all$taxa)

# Plot alluvial distributions:
##For Csp:
alluv.csps.all <- ggplot(data = df.csps.all, aes(axis1 = sample, axis2 = taxa, 
                                                 axis3 = gene, #axis4 = origin,
                                                 y = Freq)) +
  geom_alluvium(aes(fill = sample), 
                curve_type = "xspline", width = 0.60, linewidth = 0.3) + #0.3
  geom_stratum(alpha = .05, width = 0.60, linewidth = 0.6, #0.2
               color = "black", fill = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4, color = "black") +
  theme_void() +
  scale_fill_manual(values = colors.3) +
  theme(text = element_text(size = 17), legend.position = "none") +
  labs(fill = NULL) +
  geom_label(inherit.aes = FALSE,
             data = data.frame(x = c(1, 2, 3), y = c(101.5, 101.5, 101.5), # 4, 101.5
                               label = c('Habitat', 'Class', 'Csp gene')), #, 'Acquisition'
             aes(label = label, x = x, y = y))
## For antioxidants:
alluv.anox.all <- ggplot(data = df.anox.all, aes(axis1 = sample, axis2 = taxa, 
                                                 axis3 = gene, #axis4 = origin, 
                                                 y = Freq)) +
  geom_alluvium(aes(fill = sample), 
                curve_type = "xspline", width = 0.60, linewidth = 0.3) + #0.3
  geom_stratum(alpha = .05, width = 0.60, linewidth = 0.6, #0.2
               color = "black", fill = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4, color = "black") +
  theme_void() +
  scale_fill_manual(values = colors.3) +
  theme(text = element_text(size = 17), legend.position = "none") +
  labs(fill = NULL) +
  geom_label(inherit.aes = FALSE,
             data = data.frame(x = c(1, 2, 3), y = c(42, 42, 42), # 4, 42
                               label = c('Habitat', 
                                         'Class', 
                                         'Antioxidant gene')), #, 'Acquisition'
             aes(label = label, x = x, y = y))

# Export plots as images:
png("./graphics/alluvial-csps.png", units = "in", width = 11, height = 7, res = 300)
alluv.csps.all
dev.off()

png("./graphics/alluvial-anox.png", units = "in", width = 11, height = 7, res = 300)
alluv.anox.all
dev.off()


