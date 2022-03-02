## Covid analysis of taxonomic assignation data:

#-- Let's load the packages needed for the analysis --#
library("phyloseq")
library("ggplot2")
library("edgeR")
library("DESeq2")
library("pheatmap")
library("readr")
library("tidyr")
library("purrr")
library("ape")
library("plyr")
library("dplyr")
library(RColorBrewer)
library(ggpubr)
library(vegan)
library("picante")

setwd("D:/documentos/sideprojects/covid/from-reads")
#-- First, we need to import our data with phyloseq command --#                                                    --#

covid <- import_biom("non-human/covid.biom")
#-- Now, we will get trim our data for the analysis: --#
covid@tax_table@.Data <- substring(covid@tax_table@.Data, 4)
colnames(covid@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#-- Also, we need to load the metadata from our file --#
meta <- read.csv("metadata-covid3.csv", row.names = 1)
colnames(covid@otu_table@.Data) <- row.names(meta)
#-- Let's tell R that the `meta` object will be part of a phyloseq object, and merge them --#
meta <- sample_data(meta)
covid <- merge_phyloseq(covid,meta)

#-- Removing non bacterial OTUs --#
sample_sums(covid)

summary(covid@tax_table@.Data[,1] == "Bacteria")
summary(covid@tax_table@.Data[,5] != "mitocondria")
summary(covid@tax_table@.Data[,3] != "Chloroplast")

covid <- subset_taxa(covid, Kingdom == "Bacteria")

covid <- subset_taxa(covid, Family != "mitocondria" & Class != "Chloroplast")


#-- Analyzing the deepness of the sequenciation --#

deepn <- data.frame(
  samples=as.character(map(strsplit(colnames(covid@otu_table@.Data), "_"),1)),
  reads= sample_sums(covid))
ggplot(data = deepn, aes(y= reads, x = samples))+
  geom_bar(stat="identity", fill="violetred4") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept = 4242198, col= "cyan3", size = 1.5, alpha = 0.5) +
  xlab("Samples") + ylab("Number of reads") +
  ggtitle("Depth of covid libraries") +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))
#¡# By this result, we see that the sample SS09 is underepresented so we will
#¡# leave it out of the analysis
#-- Redefining the objets to comply with this result: --#
meta <- meta[-13,]
covid.otu <- covid@otu_table@.Data
covid.otu <- as.data.frame(covid.otu)
covid.otu <- select(covid.otu, -13)
colnames(covid.otu) <- row.names(meta)
covid.otu <- otu_table(covid.otu, taxa_are_rows = TRUE)
covid.tax <- covid@tax_table@.Data
covid.tax <- tax_table(covid.tax)
covid <- merge_phyloseq(covid.otu, covid.tax,meta)

#-- Defining the normalization method --#
edgeRnorm = function(physeq, ...) {
  require("edgeR")
  require("phyloseq")
  # physeq = simlist[['55000_1e-04']] z0 = simlisttmm[['55000_1e-04']] physeq
  # = simlist[['1000_0.2']] z0 = simlisttmm[['1000_0.2']] Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # See if adding a single observation, 1, everywhere (so not zeros) prevents
  # errors without needing to borrow and modify calcNormFactors (and its
  # dependent functions) It did. This fixed all problems.  Can the 1 be
  # reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  # Don't need the following additional steps, which are also built-in to some
  # of the downstream distance methods. z1 = estimateCommonDisp(z) z2 =
  # estimateTagwiseDisp(z1)
  return(z)
}
z<- edgeRnorm(covid, method = "TMM")
z1<-edgeRnorm(covid, method = "RLE")
z2<-edgeRnorm(covid, method = "upperquartile")

#-- Defining the new phyloseq object --#
covid.otu2 <- otu_table(z@.Data[[1]], taxa_are_rows = TRUE)
#-- Merging all the objects in the new normalized phyloseq object --#
covid2r <- merge_phyloseq(covid.otu2, covid.tax,meta)

#-- Let's transform the read counts in relative abundances --#
covid2 <- transform_sample_counts(physeq = covid2r, function(x) x*100/sum(x))

#-- In order to use vegan for the multivariate analysis, we will extract our data from the phyloseq objects --#
d.covid <- t(covid2@otu_table@.Data)

#-- Obtainment of beta diversity with NMDS --#
covid.ord <- ordinate(physeq = covid2r,method = "PCoA", distance = "bray")
#-- Let's plot the beta diversity analysis --#
plot_ordination(covid2r, covid.ord, color = "Covid") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Covid),
               type = "t", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  theme_bw()+
  geom_point(size=4, alpha=0.6) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "PCoA - Bray-Curtis")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=14))
#-- By Family clusterization --#
plot_ordination(covid2, covid.ord, color = "Family") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Family),
               type = "norm", linetype = 5,size = 2) +
  theme_bw()+
  geom_point(size=4, alpha=0.8) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis") 
#-- By Sympthoms clusterization --#
plot_ordination(covid2, covid.ord, color = "Symptoms") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Symptoms),
               type = "norm", linetype = 5,size = 2) +
  scale_fill_manual(values=c("#35978f", "#762a83", "#1b7837"))+
  theme_bw()+
  geom_point(size=5, alpha=0.8) +
  scale_color_manual(values=c("#35978f", "#762a83", "#1b7837")) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis") 

plot_ordination(covid2, covid.ord, color = "Sex") +
  stat_ellipse(geom = "polygon", alpha= 0.2,aes(fill=Sex),
               type = "norm", linetype = 5,size = 2) +
  theme_bw()+
  geom_point(size=5, alpha=0.8) +
  theme_bw()+
  labs(title = "NMDS - Bray-Curtis")

#-- Getting the data needed for the multivariate analysis --#
meta.covid <- data.frame(Age = as.factor(covid2@sam_data@.Data[[5]]),
                         Sex = as.factor(covid2@sam_data@.Data[[2]]),
                         Pacient = as.factor(covid2@sam_data@.Data[[3]]),
                         Symt = as.factor(covid2@sam_data@.Data[[6]]),
                         Diag = as.factor(covid2@sam_data@.Data[[7]]),
                         Covid = as.factor(covid2@sam_data@.Data[[10]]),
                         Family = as.factor(covid@sam_data@.Data[[8]]))
#-- Multivariate analysis --#
adonis(d.covid ~ Diag * Symt, data = meta.covid, permutations = 999, strata = meta.covid$Covid)
adonis(d.covid ~ Diag , data = meta.covid, permutations = 999)
adonis(d.covid ~ Symt , data = meta.covid, permutations = 999)
adonis(d.covid ~ Age , data = meta.covid, permutations = 999)
adonis(d.covid ~ Covid , data = meta.covid, permutations = 999)
adonis(d.covid ~ Family, data = meta.covid, permutations = 999)
## A combined analysis of multivariance ##
adonis(d.covid ~ Family*Symt, data = meta.covid, permutations = 999)

dis <- vegdist(d.covid, method = "bray")
mod <- betadisper(d = dis, group = meta.covid$Symt)
anova(mod)
plot(mod)

#N According to the MAGs annotation from Professor Nelly Selem, we have some taxa that made
#N the majority of the read counts, we will extract this Genera of microbes, and do the
#N same analysis as with the entire dataset

#-- Assembling the data for the pheatmap --##
covid.gene <- tax_glom(covid2r,taxrank = rank_names(covid2)[6])
covid.gene.top <- filter_taxa(covid.gene, function(x) mean(x) > 5000, TRUE)

covid.frame <- as.data.frame(covid.gene.top@otu_table@.Data)
rownames(covid.frame) <- covid.gene.top@tax_table@.Data[,6]

top.log <- log10(covid.frame)
top.log[top.log=="Inf"] <- 0
top.log[top.log=="-Inf"] <- 0

breaksList = seq(2, 7, by = .5)

my.col <- data.frame(Symptoms=meta$Symptoms, row.names = rownames(meta))
my.col$Family <- meta$Family
my.col$PCR <- meta$Covid
my.col$Same <- meta$Same

my_color <- list(Symptoms = c(`Yes`= "deeppink4", `No`= "aquamarine4"),
                 Family = c(`Fam1`="deepskyblue3" , `Fam2`="indianred4", 
                            `Fam3`="purple4", `Na`="#bababa"),
                 PCR = c(`Positive`="deeppink3", `P-Negative`="turquoise4", `Negative`="turquoise"),
                 Same = c(`0`="#bababa", `1`="brown4" , `2`="#377eb8", 
                          `3`="darkolivegreen4", `4`="#984ea3", `5`="darkgoldenrod3"))

pheatmap(top.log,
         color = (c("#ffffff","#fff7fb","#ece2f0","#d0d1e6","#a6bddb","#67a9cf","#3690c0",
                    "#02818a","#016c59","#014636")), 
        cluster_cols = TRUE, 
         cutree_cols = 3, cutree_rows = 5, border_color ="#000000",
         annotation_col = my.col, annotation_colors = my_color,)
