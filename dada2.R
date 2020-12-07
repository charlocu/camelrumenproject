set.seed(100)
setwd("/media/ankit/camel-rumen/BMJ-single-replicate-renamed")
ggplot2::theme_set(ggplot2::theme_grey())


library(dada2); packageVersion("dada2")

#####RUN-1#####
path <- "run-1/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE,  trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

errF <- learnErrors(filtFs, multithread = 30)
errR <- learnErrors(filtRs, multithread = 30)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread = 30)
dadaRs <- dada(filtRs, err=errR, multithread = 30)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:429]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution in run-1", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2, "run1.rds")

##############Just to check run-wise

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats-run1.txt", sep = "\t", quote = FALSE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_run-1.RData")


#####RUN-2#####
path <- "run-2/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE,  trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

errF <- learnErrors(filtFs, multithread = 30)
errR <- learnErrors(filtRs, multithread = 30)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread = 30)
dadaRs <- dada(filtRs, err=errR, multithread = 30)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:429]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution in run-2", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2, "run2.rds")

##############Just to check run-wise

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats-run2.txt", sep = "\t", quote = FALSE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_run-2.RData")


#####RUN-3#####
path <- "run-3/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE,  trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

errF <- learnErrors(filtFs, multithread = 30)
errR <- learnErrors(filtRs, multithread = 30)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread = 30)
dadaRs <- dada(filtRs, err=errR, multithread = 30)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:429]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution in run-3", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2, "run3.rds")


##############Just to check run-wise

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats-run3.txt", sep = "\t", quote = FALSE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_run-3.RData")

#####RUN-4#####
path <- "run-4/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, rm.phix=TRUE,  trimLeft = c(17,21), compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

errF <- learnErrors(filtFs, multithread = 30)
errR <- learnErrors(filtRs, multithread = 30)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread = 30)
dadaRs <- dada(filtRs, err=errR, multithread = 30)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:429]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), main = "Histogram of merged read length distribution in run-4", xlab = "Length", ylab = "Number of reads", col = "skyblue1")

saveRDS(seqtab2, "run4.rds")

##############Just to check run-wise

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(x = track, file = "reads-stats-run4.txt", sep = "\t", quote = FALSE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_run-4.RData")
#starting with fresh environment for merging and taxonomy assignment
rm(list=ls())


#####MERGE####
# Merge multiple runs 
st1 <- readRDS("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/run1.rds")
st2 <- readRDS("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/run2.rds")
st3 <- readRDS("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/run3.rds")
st4 <- readRDS("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/run4.rds")
st.all <- mergeSequenceTables(st1, st2, st3, st4)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_merged.RData")

# Remove chimeras
seqtab.merged <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

saveRDS(seqtab.merged, "/media/ankit/camel-rumen/BMJ-single-replicate-renamed/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_merged.RData")

# To check final sequences analysed and compare it with run-wise
rowSums(seqtab.merged)

# Assign taxonomy
taxGTDB89 <- assignTaxonomy(seqtab.merged, "/media/ankit/camel-rumen/DADA2-db/GTDB_bac120_arc122_ssu_r89.fa", multithread=TRUE)
# Write to disk
saveRDS(taxGTDB89, "/media/ankit/camel-rumen/BMJ-single-replicate-renamed/taxGTDB89.rds") 
save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_merged.RData")


#ADD SPECIES
speciesGTDB89 <- addSpecies(taxGTDB89, "/media/ankit/camel-rumen/DADA2-db/GTDB_dada2_assignment_species.fa")
saveRDS(speciesGTDB89, "/media/ankit/camel-rumen/BMJ-single-replicate-renamed/speciesGTDB89.rds")
save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_merged.RData")


#Generating a phyloseq object
library(phyloseq)

#importing metadata file
metadata <- read.delim("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/metadata.txt")
rownames(metadata) = metadata$Names

#loading in the phyloseq
ps <- phyloseq(otu_table(seqtab.merged, taxa_are_rows=FALSE), sample_data(metadata), tax_table(speciesGTDB89))


#renaming to ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(object = ps, file = "ps.RDS")
save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/dada2_merged.RData")

#clear the environment before starting downstream process
rm(list=ls())

library(ggplot2)
library(microbiome)
library(vegan)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggConvexHull)
library(grid)
library(data.table)
library(randomcoloR)
library(scales)
library(tidyr)
library(UpSetR)
library(Vennerable)
library(viridis)
library(ComplexHeatmap)

ps <- readRDS("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/ps.RDS")
sort(sample_sums(ps))

#removing ASVs with toal count less than 6
ps1 <- prune_taxa(taxa_sums(ps) > 5, ps)
ps1
#removing taxa observed in less than 6 samples
ps2 <- filter_taxa(ps1, function(x){sum(x > 0) > 5}, prune = TRUE)
ps2

saveRDS(object = ps2, file = "ps2.RDS")

#making defaults for saving images in A4 size
width = 8.27
height = 11.69

#####ALPHA DIVERSITY####
#measuring indices
qalpha <- plot_richness(ps2, measures = c("Observed", "Shannon"))
dfalpha <- as.data.frame(qalpha$data)
dfalpha$distance <- paste(dfalpha$variable,"-",dfalpha$Fraction," fraction", sep = "")
comparisonsList <- list(c("Collection-1","Collection-2"),c("Collection-1","Collection-3"), c("Collection-1","Collection-4"), c("Collection-1","Collection-5"), c("Collection-2","Collection-3"), c("Collection-3","Collection-4"), c("Collection-4","Collection-5"))
#plotting Observed ASVs
palphaObserved <- ggplot(dfalpha[dfalpha$variable == "Observed",], aes(Collection, value, color=Breed))  + geom_point(alpha = 0) + facet_grid(Fraction~Feed, scales="free")  +  geom_boxplot(position=position_dodge(1), alpha = 0.2) + geom_jitter(position=position_dodge(1)) + ggtitle("Alpha diversity measures") + theme(plot.title = element_text(hjust=0.5)) + stat_compare_means(label = "p.format", size =1.8, label.y = 2100, method = "wilcox.test") + stat_compare_means(label = "p.signif", comparisons = comparisonsList, size = 2.8, tip.length = 0.04,label.y = c(2400, 2750,3100,3450,2400,2400,2400), method = "wilcox.test") + expand_limits(y = 400)  + theme(legend.position="top") + ylab("Observed ASVs") + xlab("") + stat_compare_means( data = dfalpha[dfalpha$variable == "Observed",], mapping = aes(Collection, value) , inherit.aes = FALSE, label.x = 1.6, label.y = 3850, size =3, method = "kruskal.test") + theme(axis.title.x = element_blank() , axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(plot.margin = unit(c(5,5,1,5), "pt"))
#plotting Shannon
palphaShannon <- ggplot(dfalpha[dfalpha$variable == "Shannon",], aes(Collection, value, color=Breed))  + geom_point(alpha = 0) + facet_grid(Fraction~Feed, scales="free")  +  geom_boxplot(position=position_dodge(1), alpha = 0.2) + geom_jitter(position=position_dodge(1)) + stat_compare_means(label = "p.format" , size = 1.8, label.y = 7, method = "wilcox.test") + stat_compare_means(label = "p.signif", comparisons = comparisonsList, size = 2.8, tip.length = 0.04, label.y = c(7.2,7.5,7.8,8.1,7.2,7.2,7.2), method = "wilcox.test") + theme(legend.position="none") + ylab("Shannon Index") + xlab("Collection") + stat_compare_means( data = dfalpha[dfalpha$variable == "Shannon",], mapping = aes(Collection, value) , inherit.aes = FALSE, label.x = 1.6, label.y = 8.5, size =3, method = "kruskal.test") + theme(strip.background.x = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) + theme(plot.margin = unit(c(1,5,5,21), "pt"))
#combining the plots
palpha <- ggarrange(palphaObserved, palphaShannon, ncol =1, align = "v")
palpha
ggsave(filename = "alpha-plot.pdf", plot = palpha, device = "pdf", width = width, height = height, units = "in")
#comparing among feed, fraction, collection and breed and saving the results in a file
alpha_means_feed <- compare_means(formula = value ~ Feed, data = dfalpha, method = "kruskal.test", group.by = "variable", p.adjust.method = "BH" )
alpha_means_fraction <- compare_means(formula = value ~ Fraction, data = dfalpha, method = "kruskal.test", group.by = "variable", p.adjust.method = "BH" )
alpha_means_breed <- compare_means(formula = value ~ Breed, data = dfalpha, method = "kruskal.test", group.by = "variable", p.adjust.method = "BH" )
alpha_means_collection <- compare_means(formula = value ~ Collection, data = dfalpha, method = "kruskal.test", group.by = "variable", p.adjust.method = "BH" )
cat("Kruskal-wallis test between variables for Observed ASVs and Shannon Diversity\n", file = "alpha-diversity-tests.txt")
cat("\nAmong different Feeds\n", file = "alpha-diversity-tests.txt", append = TRUE)
capture.output(alpha_means_feed, file = "alpha-diversity-tests.txt", append = TRUE)
cat("\nBetween different Fractions\n", file = "alpha-diversity-tests.txt", append = TRUE)
capture.output(alpha_means_fraction, file = "alpha-diversity-tests.txt", append = TRUE)
cat("\nBetween different Breeds\n", file = "alpha-diversity-tests.txt", append = TRUE)
capture.output(alpha_means_breed, file = "alpha-diversity-tests.txt", append = TRUE)
cat("\nAmong different Collections\n", file = "alpha-diversity-tests.txt", append = TRUE)
capture.output(alpha_means_collection, file = "alpha-diversity-tests.txt", append = TRUE)
cat("\n\n", file = "alpha-diversity-tests.txt", append = TRUE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/phyloseq.RData")

#####BETA-DIVERSITY####
#calculating relative frequency
ps2.rel <- microbiome::transform(ps2, "compositional")
#subsetting according to different groups for easier downstream processing
Liquidps.rel <- subset_samples(ps2.rel, Fraction == "Liquid")
Solidps.rel <- subset_samples(ps2.rel, Fraction == "Solid")
Bajraps.rel <- subset_samples(ps2.rel, Feed == "Bajra")
Jowarps.rel <- subset_samples(ps2.rel, Feed == "Jowar")
Makaips.rel <- subset_samples(ps2.rel, Feed == "Makai")
#Exctracting metadata in dataframe
metadata <- data.frame(sample_data(ps2.rel))
metadata.liquid <- data.frame(sample_data(Liquidps.rel))
metadata.solid <- data.frame(sample_data(Solidps.rel))
metadata.bajra <- data.frame(sample_data(Bajraps.rel))
metadata.jowar <- data.frame(sample_data(Jowarps.rel))
metadata.makai <- data.frame(sample_data(Makaips.rel))


#Comparison among ALL samples
ps2.bray.dist <- phyloseq::distance(ps2.rel, method="bray")
permanova.bray <- vegan::adonis(ps2.bray.dist ~ Fraction + Feed + Collection + Breed, data = metadata)
permanova.bray
ps2.bray.ordination = ordinate(ps2.rel, method="NMDS", distance=ps2.bray.dist)
ps2.stress.label <- paste("NMDS Stress = ", round(ps2.bray.ordination$stress, 3))
#plot ordination
betaplot.ps2 <- phyloseq::plot_ordination(ps2.rel, ps2.bray.ordination, color= "Collection", shape = "Group.Fp") + theme(aspect.ratio=1) + labs(color = "Collection", shape = "Feed-Fraction") + geom_point(size=4) + scale_shape_manual(values=c(0,5,2,6,3,4)) + stat_ellipse(geom = "polygon", aes(fill = Fraction, x = NMDS1, y = NMDS2), inherit.aes = FALSE, alpha=0.1, show.legend = FALSE)  + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of all samples") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.ps2 <- betaplot.ps2 + annotation_custom(grob = textGrob(label = ps2.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = -0.45, xmax = -0.45, ymin = 0.27, ymax = 0.27)
betaplot.ps2
leg <- get_legend(betaplot.ps2)


#Comparison among LIQUID samples
liquid.bray.dist <- phyloseq::distance(Liquidps.rel, method="bray")
permanova.bray.liquid <- vegan::adonis(liquid.bray.dist ~ Feed + Collection + Breed, data = metadata.liquid)
permanova.bray.liquid
liquid.bray.ordination = ordinate(Liquidps.rel, method="NMDS", distance=liquid.bray.dist)
liquid.stress.label <- paste("NMDS Stress = ", round(liquid.bray.ordination$stress, 3))
#plot ordination
betaplot.liquid <- phyloseq::plot_ordination(Liquidps.rel, liquid.bray.ordination, color= "Collection", shape = "Group.Fp") + theme(aspect.ratio=1) + labs(color = "Collection", shape = "Feed-Fraction") + geom_point(size=4) + scale_shape_manual(values=c(0,2,3)) + geom_convexhull(inherit.aes=FALSE, alpha = 0.1, show.legend = FALSE, mapping = aes(x = NMDS1, y= NMDS2,fill = Group.Cp, color = Collection)) + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of Liquid fraction samples") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.liquid <- betaplot.liquid + annotation_custom(grob = textGrob(label = liquid.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = 0.12, xmax = 0.12, ymin = 0.3, ymax = 0.3)
betaplot.liquid

#Comparison among SOLID samples
solid.bray.dist <- phyloseq::distance(Solidps.rel, method="bray")
permanova.bray.solid <- vegan::adonis(solid.bray.dist ~ Feed + Collection + Breed, data = metadata.solid)
permanova.bray.solid
solid.bray.ordination = ordinate(Solidps.rel, method="NMDS", distance=solid.bray.dist)
solid.stress.label <- paste("NMDS Stress = ", round(solid.bray.ordination$stress, 3))
#plot ordination
betaplot.solid <- phyloseq::plot_ordination(Solidps.rel, solid.bray.ordination, color= "Collection", shape = "Group.Fp") + theme(aspect.ratio=1) + labs(color = "Collection", shape = "Feed-Fraction") + geom_point(size=4) + scale_shape_manual(values=c(5,6,4)) + geom_convexhull(inherit.aes=FALSE, alpha = 0.1, show.legend = FALSE, mapping = aes(x = NMDS1, y= NMDS2,fill = Group.Cp, color = Collection)) + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of solid fraction samples") + theme(plot.title = element_text(size = 10, hjust = 0.5)) + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.solid <- betaplot.solid + annotation_custom(grob = textGrob(label = solid.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = 0.12, xmax = 0.12, ymin = 0.22, ymax = 0.22)
betaplot.solid

#Comparison among BAJRA samples
bajra.bray.dist <- phyloseq::distance(Bajraps.rel, method="bray")
permanova.bray.bajra <- vegan::adonis(bajra.bray.dist ~ Fraction + Collection + Breed, data = metadata.bajra)
permanova.bray.bajra
bajra.bray.ordination = ordinate(Bajraps.rel, method="NMDS", distance=bajra.bray.dist)
bajra.stress.label <- paste("NMDS Stress = ", round(bajra.bray.ordination$stress, 3))
#plot ordination
betaplot.bajra <- phyloseq::plot_ordination(Bajraps.rel, bajra.bray.ordination, color = "Collection", shape = "Group.Fp") + theme(aspect.ratio=1)  + geom_point(size=4) + scale_shape_manual(values=c(0,5)) + geom_convexhull(inherit.aes=FALSE, alpha = 0.1, show.legend = FALSE, mapping = aes(x = NMDS1, y= NMDS2,fill = Group.Cp, color = Collection)) + scale_fill_manual(values = c("Red","Red", "chartreuse4","chartreuse4", "cyan4", "cyan4", "Purple", "Purple", "darkorange1", "darkorange1")) + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of Bajra fed samples") + theme(plot.title = element_text(size = 10, hjust = 0.5))+ labs(color = "Collection", shape = "Feed-Fraction") + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.bajra <- betaplot.bajra + annotation_custom(grob = textGrob(label = bajra.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = 0.2, xmax = 0.2, ymin = 0.2, ymax = 0.2)
betaplot.bajra

#Comparison among JOWAR samples
jowar.bray.dist <- phyloseq::distance(Jowarps.rel, method="bray")
permanova.bray.jowar <- vegan::adonis(jowar.bray.dist ~ Fraction +  Collection + Breed, data = metadata.jowar)
permanova.bray.jowar
jowar.bray.ordination = ordinate(Jowarps.rel, method="NMDS", distance=jowar.bray.dist)
jowar.stress.label <- paste("NMDS Stress = ", round(jowar.bray.ordination$stress, 3))
#plot ordination
betaplot.jowar <- phyloseq::plot_ordination(Jowarps.rel, jowar.bray.ordination, color = "Collection", shape = "Group.Fp") + theme(aspect.ratio=1)  + geom_point(size=4) + scale_shape_manual(values=c(2,6)) + geom_convexhull(inherit.aes=FALSE, alpha = 0.1, show.legend = FALSE, mapping = aes(x = NMDS1, y= NMDS2,fill = Group.Cp, color = Collection)) + scale_fill_manual(values = c("Red","Red", "chartreuse4","chartreuse4", "cyan4", "cyan4", "Purple", "Purple", "darkorange1", "darkorange1")) + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of Jowar fed samples") + theme(plot.title = element_text(size = 10, hjust = 0.5))+ labs(color = "Collection", shape = "Feed-Fraction") + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.jowar <- betaplot.jowar + annotation_custom(grob = textGrob(label = jowar.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = 0, xmax = 0, ymin = 0.15, ymax = 0.15)
betaplot.jowar


#Comparison among MAKAI samples
makai.bray.dist <- phyloseq::distance(Makaips.rel, method="bray")
permanova.bray.makai <- vegan::adonis(makai.bray.dist ~ Fraction +  Collection + Breed, data = metadata.makai)
permanova.bray.makai
makai.bray.ordination = ordinate(Makaips.rel, method="NMDS", distance=makai.bray.dist)
makai.stress.label <- paste("NMDS Stress = ", round(makai.bray.ordination$stress, 3))
#plot ordination
betaplot.makai <- phyloseq::plot_ordination(Makaips.rel, makai.bray.ordination, color = "Collection", shape = "Group.Fp") + theme(aspect.ratio=1)  + geom_point(size=4) + scale_shape_manual(values=c(3,4)) + geom_convexhull(inherit.aes=FALSE, alpha = 0.1, show.legend = FALSE, mapping = aes(x = NMDS1, y= NMDS2,fill = Group.Cp, color = Collection)) + scale_fill_manual(values = c("Red","Red", "chartreuse4","chartreuse4", "cyan4", "cyan4", "Purple", "Purple", "darkorange1", "darkorange1")) + ggtitle("NMDS plot using Bray-Curtis distance \n for relative abundance of Makai fed samples") + theme(plot.title = element_text(size = 10, hjust = 0.5))+ labs(color = "Collection", shape = "Feed-Fraction")  + scale_color_manual(values = c("Red", "chartreuse4", "cyan4",  "Purple", "darkorange1"))
#add NMDS stress on the graph....adjust the position based on the actual graph. 
betaplot.makai <- betaplot.makai + annotation_custom(grob = textGrob(label = makai.stress.label, hjust = 0, gp = gpar(fontsize = 10, fontface = "bold")), xmin = 0.12, xmax = 0.12, ymin = 0.3, ymax = 0.3)
betaplot.makai

beta1 <- ggarrange(betaplot.ps2, betaplot.liquid, betaplot.solid, ncol = 3, nrow = 1, align = "hv", legend = "none")
beta2 <- ggarrange(betaplot.bajra, betaplot.jowar, betaplot.makai, ncol = 3, nrow = 1, align = "hv", legend = "none")
beta3 <- ggarrange(ggarrange(beta1, beta2, nrow = 2, align = "hv", legend = "none", ncol = 1), as_ggplot(leg), ncol = 2, nrow = 1, widths = c(9, 1))
ggsave(filename = "beta-plot-all.pdf", plot = beta3, device = "pdf", width = height*1.5, height = width*1.5, units = "in", dpi = 300)

cat("PERMANOVA test among various samples on different parameters\n", file = "beta-diversity-permanova.txt")
cat("\nAmong ALL samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray, file = "beta-diversity-permanova.txt", append = TRUE)
cat("\n\nAmong LIQUID fraction samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray.liquid, file = "beta-diversity-permanova.txt", append = TRUE)
cat("\n\nAmong SOLID fraction samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray.solid, file = "beta-diversity-permanova.txt", append = TRUE)
cat("\n\nAmong BAJRA feed samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray.bajra, file = "beta-diversity-permanova.txt", append = TRUE)
cat("\n\nAmong JOWAR feed samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray.jowar, file = "beta-diversity-permanova.txt", append = TRUE)
cat("\n\nAmong MAKAI feed samples\n", file = "beta-diversity-permanova.txt", append = TRUE)
capture.output(permanova.bray.makai, file = "beta-diversity-permanova.txt", append = TRUE)

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/phyloseq.RData")


####COMPARING & PLOTTING TAXONOMY ####

########plotting taxonomy
#####PHYLUM
# define the levels to glom 
ps2.Phylum.rel <- tax_glom(physeq = ps2.rel, taxrank = "Phylum", NArm = FALSE)
#adding phylum names and changing names of NA to Unclassified taxa
taxa_names(ps2.Phylum.rel) <- tax_table(ps2.Phylum.rel)[, 2]
taxa_names(ps2.Phylum.rel)[is.na(taxa_names(ps2.Phylum.rel))] <- "Unclassified Phylum"
#create dataframe
ps2.Phylum.rel.df <- data.table(psmelt(ps2.Phylum.rel))
#group df by taxa and calculate mean rel. abundance
ps2.Phylum.rel.df[, avg := mean(Abundance, na.rm = TRUE), by = "OTU"]
# Change name of taxa less than 1%, storing the names in new column
ps2.Phylum.rel.df[(avg <= 0.01), OTU := "Less Abundant Phyla"]
# Creating df with summarized lesser abundant taxa abundance
ps2.Phylum.rel.df <- ps2.Phylum.rel.df[, sum(Abundance), by = list(OTU, Sample, Feed, Collection, Breed, Fraction, Group.Fc, Group.Fb, Group.Fp, Group.Cb, Group.Cp, Group.Bp, Group.Fcb, Group.Fcp, Group.Cbp, Group.all, BarPlot, avg)]
colnames(ps2.Phylum.rel.df)[19] <- "Abundance"
#get color codes, store in object and remember to not run that code again or different colors will generate everytime
colcodes.phylum <- distinctColorPalette(length(unique(ps2.Phylum.rel.df$OTU))+10)
#plotting stacked bar chart
Phylum.bar <- ggplot(data=ps2.Phylum.rel.df, aes(x=BarPlot, y=Abundance*100, fill=OTU)) + geom_bar(position="stack", stat="identity", color = "black", size = 0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Phyla")   + scale_fill_manual(values = colcodes.phylum) + labs(fill = "Phylum")  + coord_cartesian(expand = FALSE ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5))) + facet_grid(Fraction ~ Feed , scales="free") + theme(axis.text.x = element_text(angle = 90)) + geom_vline(xintercept=c(4.5,8.5,12.5,16.5), color= "red")
Phylum.bar

#####Genus
# # define the levels to glom 
ps2.Genus.rel <- tax_glom(physeq = ps2.rel, taxrank = "Genus", NArm = FALSE)
#melt to dataframe and convert factor to character
ps2.Genus.rel.df <- data.table(psmelt(ps2.Genus.rel))
ps2.Genus.rel.df$Kingdom <- as.character(ps2.Genus.rel.df$Kingdom)
ps2.Genus.rel.df$Phylum <- as.character(ps2.Genus.rel.df$Phylum)
ps2.Genus.rel.df$Class <- as.character(ps2.Genus.rel.df$Class)
ps2.Genus.rel.df$Order <- as.character(ps2.Genus.rel.df$Order)
ps2.Genus.rel.df$Family <- as.character(ps2.Genus.rel.df$Family)
ps2.Genus.rel.df$Genus <- as.character(ps2.Genus.rel.df$Genus)
#NA entries in Genera were renamed to "highest annotated level"_X
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Phylum <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Phylum),]$Kingdom, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Class <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Class),]$Phylum, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Order <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Order),]$Class, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Family <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Family),]$Order, "_X")
ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Genus <- paste0(ps2.Genus.rel.df[is.na(ps2.Genus.rel.df$Genus),]$Family, "_X")
# #group df by taxa and calculate mean rel. abundance
ps2.Genus.rel.df[, avg := mean(Abundance, na.rm = TRUE), by = "Genus"]
#Those not meeting these criteria were renamed as Lesser abundant genera
ps2.Genus.rel.df[(avg < 0.01 ), Genus := "Lesser Abundant Genera (n=385)"]
#Creating df with summarized lesser abundant taxa abundance
ps2.Genus.rel.df <- ps2.Genus.rel.df[, sum(Abundance), by = list(Genus, Sample, Feed, Collection, Breed, Fraction, Group.Fc, Group.Fb, Group.Fp, Group.Cb, Group.Cp, Group.Bp, Group.Fcb, Group.Fcp, Group.Cbp, Group.all, BarPlot, avg)]
colnames(ps2.Genus.rel.df)[19] <- "Abundance"
# get color codes, store in object and remember to not run that code again or different colors everytime
colcodes.genus <- distinctColorPalette(length(unique(ps2.Genus.rel.df$Genus))+9)
#plotting stacked bar chart
Genus.bar <- ggplot(data=ps2.Genus.rel.df, aes(x=BarPlot, y=Abundance*100, fill=Genus)) + geom_bar(position="stack", stat="identity", color = "black", size=0.05, width = 1)  + xlab("Samples") + ylab("Relative abundance of Genera")   + scale_fill_manual(values = colcodes.genus) + labs(fill = "Genus") + theme(axis.text.x = element_text(angle = 90)) + coord_cartesian(expand = FALSE, ylim = c(0, 100) ) + guides(fill=guide_legend(ncol=1)) + scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 5))) + facet_grid(Fraction ~ Feed , scales="free") + theme(axis.text.x = element_text(angle = 90)) + geom_vline(xintercept=c(4.5,8.5,12.5,16.5), color= "red")
Genus.bar

##arranging both plots together
bar.plot <- ggarrange(Phylum.bar, Genus.bar, ncol = 1, nrow = 2, align = "hv")
ggsave(filename = "bar-plot-l.pdf", plot = bar.plot, device = "pdf",dpi = 300, width = width, height = height, units = "in")


####COMPARING phyla and genera between six groups for significance
#PHYLUM
Phylum.rel.df <- data.table(psmelt(ps2.Phylum.rel))
#GENUS
#melt existing genus level agglomerated phyloseq object to dataframe and convert factors to character
Genus.rel.df <- data.table(psmelt(ps2.Genus.rel))
Genus.rel.df$Kingdom <- as.character(Genus.rel.df$Kingdom)
Genus.rel.df$Phylum <- as.character(Genus.rel.df$Phylum)
Genus.rel.df$Class <- as.character(Genus.rel.df$Class)
Genus.rel.df$Order <- as.character(Genus.rel.df$Order)
Genus.rel.df$Family <- as.character(Genus.rel.df$Family)
Genus.rel.df$Genus <- as.character(Genus.rel.df$Genus)
#NA entries in Genera were renamed to "highest annotated level"_X
Genus.rel.df[is.na(Genus.rel.df$Phylum),]$Phylum <- paste0(Genus.rel.df[is.na(Genus.rel.df$Phylum),]$Kingdom, "_X")
Genus.rel.df[is.na(Genus.rel.df$Class),]$Class <- paste0(Genus.rel.df[is.na(Genus.rel.df$Class),]$Phylum, "_X")
Genus.rel.df[is.na(Genus.rel.df$Order),]$Order <- paste0(Genus.rel.df[is.na(Genus.rel.df$Order),]$Class, "_X")
Genus.rel.df[is.na(Genus.rel.df$Family),]$Family <- paste0(Genus.rel.df[is.na(Genus.rel.df$Family),]$Order, "_X")
Genus.rel.df[is.na(Genus.rel.df$Genus),]$Genus <- paste0(Genus.rel.df[is.na(Genus.rel.df$Genus),]$Family, "_X")

##
###Comparison between fration
Fraction.Genus.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Genus.rel.df, group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
Fraction.Phylum.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Phylum.rel.df, group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
##Feed-wise
# Fraction.Genus.group.Bajra.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Genus.rel.df[Feed == "Bajra"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Fraction.Phylum.group.Bajra.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Phylum.rel.df[Feed == "Bajra"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Fraction.Genus.group.Jowar.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Genus.rel.df[Feed == "Jowar"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Fraction.Phylum.group.Jowar.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Phylum.rel.df[Feed == "Jowar"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Fraction.Genus.group.Makai.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Genus.rel.df[Feed == "Makai"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Fraction.Phylum.group.Makai.means1 <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = Phylum.rel.df[Feed == "Makai"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
###


###Comparison among feed
Feed.Genus.group.liquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Fraction == "Liquid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
Feed.Phylum.group.liquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Fraction == "Liquid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
Feed.Genus.group.solid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Fraction == "Solid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
Feed.Phylum.group.solid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Fraction == "Solid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Genus.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df, group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df, group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
##Collection-wise
# Feed.Genus.group.C1.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Collection == "Collection-1"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Genus.group.C2.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Collection == "Collection-2"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Genus.group.C3.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Collection == "Collection-3"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Genus.group.C4.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Collection == "Collection-4"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Genus.group.C5.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Genus.rel.df[Collection == "Collection-5"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.C1.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Collection == "Collection-1"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.C2.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Collection == "Collection-2"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.C3.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Collection == "Collection-3"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.C4.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Collection == "Collection-4"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Feed.Phylum.group.C5.means1 <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = Phylum.rel.df[Collection == "Collection-5"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
###


###Comparison among collections
Collection.Genus.group.liquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Fraction == "Liquid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
Collection.Phylum.group.liquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Fraction == "Liquid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
Collection.Genus.group.solid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Fraction == "Solid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
Collection.Phylum.group.solid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Fraction == "Solid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df, group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df, group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
##Group-wise
Collection.Genus.group.liquid.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Fraction == "Liquid"], group.by = "Genus", method = "wilcox.test", p.adjust.method = "BH"))
Collection.Phylum.group.liquid.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Fraction == "Liquid"], group.by = "OTU", method = "wilcox.test", p.adjust.method = "BH"))
Collection.Genus.group.solid.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Fraction == "Solid"], group.by = "Genus", method = "wilcox.test", p.adjust.method = "BH"))
Collection.Phylum.group.solid.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Fraction == "Solid"], group.by = "OTU", method = "wilcox.test", p.adjust.method = "BH"))
Collection.Genus.group.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df, group.by = "Genus", method = "wilcox.test", p.adjust.method = "BH"))
Collection.Phylum.group.means2 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df, group.by = "OTU", method = "wilcox.test", p.adjust.method = "BH"))
##Feed-wise
# Collection.Genus.group.Bajraliquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "BajraLiquid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.Jowarliquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "JowarLiquid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.Makailiquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "MakaiLiquid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Bajraliquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "BajraLiquid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Jowarliquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "JowarLiquid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Makailiquid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "MakaiLiquid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.Bajrasolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "BajraSolid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.Jowarsolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "JowarSolid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Genus.group.Makaisolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Genus.rel.df[Group.Fp == "MakaiSolid"], group.by = "Genus", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Bajrasolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "BajraSolid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Jowarsolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "JowarSolid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
# Collection.Phylum.group.Makaisolid.means1 <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = Phylum.rel.df[Group.Fp == "MakaiSolid"], group.by = "OTU", method = "kruskal.test", p.adjust.method = "BH"))
###

#####export the result and summarized in Excel
lst <- list("Collection.Genus.group.Bajraliquid.means1" = Collection.Genus.group.Bajraliquid.means1, "Collection.Genus.group.Bajrasolid.means1" = Collection.Genus.group.Bajrasolid.means1, "Collection.Genus.group.Jowarliquid.means1" = Collection.Genus.group.Jowarliquid.means1, "Collection.Genus.group.Jowarsolid.means1" = Collection.Genus.group.Jowarsolid.means1, "Collection.Genus.group.Makailiquid.means1" = Collection.Genus.group.Makailiquid.means1, "Collection.Genus.group.Makaisolid.means1" = Collection.Genus.group.Makaisolid.means1, "Collection.Genus.group.liquid.means1" = Collection.Genus.group.liquid.means1, "Collection.Genus.group.liquid.means2" = Collection.Genus.group.liquid.means2, "Collection.Genus.group.means1" = Collection.Genus.group.means1, "Collection.Genus.group.means2" = Collection.Genus.group.means2, "Collection.Genus.group.solid.means1" = Collection.Genus.group.solid.means1, "Collection.Genus.group.solid.means2" = Collection.Genus.group.solid.means2, "Collection.Phylum.group.Bajraliquid.means1" = Collection.Phylum.group.Bajraliquid.means1, "Collection.Phylum.group.Bajrasolid.means1" = Collection.Phylum.group.Bajrasolid.means1, "Collection.Phylum.group.Jowarliquid.means1" = Collection.Phylum.group.Jowarliquid.means1, "Collection.Phylum.group.Jowarsolid.means1" = Collection.Phylum.group.Jowarsolid.means1, "Collection.Phylum.group.Makailiquid.means1" = Collection.Phylum.group.Makailiquid.means1, "Collection.Phylum.group.Makaisolid.means1" = Collection.Phylum.group.Makaisolid.means1, "Collection.Phylum.group.liquid.means1" = Collection.Phylum.group.liquid.means1, "Collection.Phylum.group.liquid.means2" = Collection.Phylum.group.liquid.means2, "Collection.Phylum.group.means1" = Collection.Phylum.group.means1, "Collection.Phylum.group.means2" = Collection.Phylum.group.means2, "Collection.Phylum.group.solid.means1" = Collection.Phylum.group.solid.means1, "Collection.Phylum.group.solid.means2" = Collection.Phylum.group.solid.means2, "Feed.Genus.group.C1.means1" = Feed.Genus.group.C1.means1, "Feed.Genus.group.C2.means1" = Feed.Genus.group.C2.means1, "Feed.Genus.group.C3.means1" = Feed.Genus.group.C3.means1, "Feed.Genus.group.C4.means1" = Feed.Genus.group.C4.means1, "Feed.Genus.group.C5.means1" = Feed.Genus.group.C5.means1, "Feed.Genus.group.liquid.means1" = Feed.Genus.group.liquid.means1, "Feed.Genus.group.means1" = Feed.Genus.group.means1, "Feed.Genus.group.solid.means1" = Feed.Genus.group.solid.means1, "Feed.Phylum.group.C1.means1" = Feed.Phylum.group.C1.means1, "Feed.Phylum.group.C2.means1" = Feed.Phylum.group.C2.means1, "Feed.Phylum.group.C3.means1" = Feed.Phylum.group.C3.means1, "Feed.Phylum.group.C4.means1" = Feed.Phylum.group.C4.means1, "Feed.Phylum.group.C5.means1" = Feed.Phylum.group.C5.means1, "Feed.Phylum.group.liquid.means1" = Feed.Phylum.group.liquid.means1, "Feed.Phylum.group.means1" = Feed.Phylum.group.means1, "Feed.Phylum.group.solid.means1" = Feed.Phylum.group.solid.means1, "Fraction.Genus.group.Bajra.means1" = Fraction.Genus.group.Bajra.means1, "Fraction.Genus.group.Jowar.means1" = Fraction.Genus.group.Jowar.means1, "Fraction.Genus.group.Makai.means1" = Fraction.Genus.group.Makai.means1, "Fraction.Genus.group.means1" = Fraction.Genus.group.means1, "Fraction.Phylum.group.Bajra.means1" = Fraction.Phylum.group.Bajra.means1, "Fraction.Phylum.group.Jowar.means1" = Fraction.Phylum.group.Jowar.means1, "Fraction.Phylum.group.Makai.means1" = Fraction.Phylum.group.Makai.means1, "Fraction.Phylum.group.means1" = Fraction.Phylum.group.means1)
for (ii in names(lst)){
  write.table(lst[[ii]], paste(ii, ".txt", sep=""), col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}

save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/phyloseq.RData")

#####CORE microbiome from Collection-5####
# Identifying core taxa by changing detection level and prevalence
JL <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "JowarLiquid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)
JS <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "JowarSolid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)
BL <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "BajraLiquid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)
ML <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "MakaiLiquid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)
BS <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "BajraSolid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)
MS <- spread(data = Genus.rel.df[Collection == "Collection-5" & Group.Fp == "MakaiSolid",c("Sample", "Abundance", "Genus")], key = Sample, value = Abundance)

#Extracting genera with prevalence > 0.5 among groups and minimum abundance 0.005 (0.5%)
BL <- as.data.frame(BL)
rownames(BL) <- BL$Genus
BL <- as.matrix(BL[,-1])
BL.taxa <- names(which(prevalence(BL, 0.005, include.lowest=TRUE) > 0.5))
BS <- as.data.frame(BS)
rownames(BS) <- BS$Genus
BS <- as.matrix(BS[,-1])
BS.taxa <- names(which(prevalence(BS, 0.005, include.lowest=TRUE) > 0.5))
JL <- as.data.frame(JL)
rownames(JL) <- JL$Genus
JL <- as.matrix(JL[,-1])
JL.taxa <- names(which(prevalence(JL, 0.005, include.lowest=TRUE) > 0.5))
JS <- as.data.frame(JS)
rownames(JS) <- JS$Genus
JS <- as.matrix(JS[,-1])
JS.taxa <- names(which(prevalence(JS, 0.005, include.lowest=TRUE) > 0.5))
ML <- as.data.frame(ML)
rownames(ML) <- ML$Genus
ML <- as.matrix(ML[,-1])
ML.taxa <- names(which(prevalence(ML, 0.005, include.lowest=TRUE) > 0.5))
MS <- as.data.frame(MS)
rownames(MS) <- MS$Genus
MS <- as.matrix(MS[,-1])
MS.taxa <- names(which(prevalence(MS, 0.005, include.lowest=TRUE) > 0.5))

mt.core = as.data.frame(cbind(c("BL", "JL", "ML", "BS", "JS", "MS"),c("Liquid", "Liquid", "Liquid", "Solid", "Solid", "Solid")))
colnames(mt.core) = c("sets", "fraction")
core.list = list(BL = BL.taxa, JL = JL.taxa, ML = ML.taxa, BS = BS.taxa, JS = JS.taxa, MS = MS.taxa)
upset(fromList(core.list), sets = c("BL", "BS", "JL", "JS", "ML", "MS"), keep.order = TRUE,  order.by = "degree", nsets = 6, set.metadata = list(data = mt.core, plots = list(list(type = "matrix_rows", column = "fraction", colors = c(Liquid = "green", Solid = "navy"), alpha = 0.2))), query.legend = "top", queries = list(list(query = intersects, params = list("BL", "JL", "ML", "BS", "JS", "MS"), color = "orange", active = T, query.name = "Present in all samples"), list(query = intersects, params = list("BL", "JL", "ML"),  color = "red", active = T, query.name = "Liquid samples specific"), list(query = intersects, params = list("BS", "JS", "MS"), active = T, color = "blue", query.name = "Solid samples specific")))

#To check Intersection level wise names of genera
venny <- Venn(core.list)
venny@IntersectionSets


save.image("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/phyloseq.RData")


#####COG####
cog.abundance.raw <- read.delim("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/cog-abundance.txt")
cog.class.abundance.raw <- read.delim("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/cog-class-abundance.txt")

#modifying abundance file for further process
cog.class.abundance <- gather(data = cog.class.abundance.raw, key = "Sample", value = "Abundance", 2:49)
cog.class.abundance <- merge(x = cog.class.abundance, metadata, all.x = TRUE, by.x = "Sample", by.y = "Names")

#Calculating significantly different COG classes by Kruskal-Wallis test
cog.class.feed.means <- as.data.frame(compare_means(formula = Abundance ~ Feed, data = cog.class.abundance, group.by = "Class", method = "kruskal.test", p.adjust.method = "BH"))
cog.class.fraction.means <- as.data.frame(compare_means(formula = Abundance ~ Fraction, data = cog.class.abundance, group.by = "Class", method = "kruskal.test", p.adjust.method = "BH"))
cog.class.collection.means <- as.data.frame(compare_means(formula = Abundance ~ Collection, data = cog.class.abundance, group.by = "Class", method = "kruskal.test", p.adjust.method = "BH"))
cog.class.breed.means <- as.data.frame(compare_means(formula = Abundance ~ Breed, data = cog.class.abundance, group.by = "Class", method = "kruskal.test", p.adjust.method = "BH"))
#####export the result and summarized in Excel
lst.cog <- list("cog.class.feed.liquid.means" = cog.class.feed.liquid.means, "cog.class.feed.solid.means" = cog.class.feed.solid.means, "cog.class.feed.means" = cog.class.feed.means, "cog.class.fraction.means" = cog.class.fraction.means, "cog.class.collection.means" = cog.class.collection.means, "cog.class.breed.means" = cog.class.breed.means)
for (ii in names(lst.cog)){
  write.table(lst.cog[[ii]], paste(ii, ".txt", sep=""), col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}
###COG CLASS heatmap
cog.class.hm.df = as.data.frame(sapply(cog.class.abundance.raw[,-1], as.numeric), row.names = as.character(cog.class.abundance.raw$Class))
cog.class.hm.df[cog.class.hm.df == 0] <- NA
#metadata file for subset of Collection 4 and 5 is prepared, considering the order of samples from class abundance table
cog.metadata <- data.frame("Names" = colnames(cog.class.hm.df))
cog.metadata <- merge(x = cog.metadata, y = metadata[,1:5], by = "Names", all.x = TRUE)
# Create the heatmap annotation
col = list(Feed = c("Bajra" = "green", "Jowar" = "gray", "Makai" = "darkred"), Fraction = c("Liquid" = "yellow", "Solid" = "orange"), Collection = c("Collection-4" = "lightblue", "Collection-5" = "purple"), Breed = c("Bikaneri" = "deeppink", "Kachchhi" = "orchid"))
ha <- HeatmapAnnotation(Fraction = cog.metadata$Fraction, Collection = cog.metadata$Collection, Feed = cog.metadata$Feed, Breed = cog.metadata$Breed, col = col)
#create heatmap and save it
pdf("COG-class-heatmap.pdf", width = width *1.2, height = height*1.2)
Heatmap(as.matrix(log10(cog.class.hm.df)), name = "log10(Abundance)", cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, rect_gp = gpar(col = "white", lwd = 1), col = rev(viridis(50)), na_col = "white", column_title = "Samples", column_title_side = "bottom", row_title = "COG Class", row_names_side = "left", bottom_annotation = ha, heatmap_legend_param = list(at = c(-1, 0, 1, 2, 3, 4, 5, 6), legend_height = unit(6, "cm"), title_position = "leftcenter-rot", grid_width = unit(1, "cm")) )  
dev.off()

### COG NMDS plot....
COG.dist.bray = vegan::vegdist(x = as.matrix(t(cog.abundance.raw[,-1])), method = "bray")
permanova.bray.COG <- vegan::adonis(COG.dist.bray ~ Fraction + Collection + Feed + Breed, data = cog.metadata)
cat("PERMANOVA test among COGs on different parameters\n", file = "cog-permanova.txt")
capture.output(permanova.bray.COG, file = "cog-permanova.txt", append = TRUE)
COG.bray.ordination = vegan::metaMDS(comm = COG.dist.bray, distance = "bray")
COG.bray.ordination.df = as.data.frame(vegan::scores(COG.bray.ordination, display = "sites"))
COG.bray.ordination.df = merge(as.data.frame(COG.bray.ordination.df), as.data.frame(metadata), by='row.names', all.x = TRUE)
COG.NMDS.plot <- ggplot(COG.bray.ordination.df, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(NMDS1, NMDS2, colour = Collection, shape = Feed:Fraction))  + theme(aspect.ratio=1) + labs(color = "Collection", shape = "Feed-Fraction") + scale_shape_manual(values=c(0,5,2,6,3,4))  + ggtitle("NMDS plot using Bray-Curtis distance for COG ids") + theme(plot.title = element_text(size = 10)) + scale_color_manual(values = c("Purple", "darkorange1")) + stat_ellipse(geom = "polygon", aes(fill = Fraction, x = NMDS1, y = NMDS2), inherit.aes = FALSE, alpha=0.1, show.legend = FALSE)
ggsave(filename = "COG-NMDS-plot.pdf", plot = COG.NMDS.plot, device = "pdf", width = width, height = height, units = "in")  


#####CAZY ####
#import table of all GH families for plotting heatmap and checking group-wise significant differences
cazy.table.GH <- read.delim("/media/ankit/camel-rumen/BMJ-single-replicate-renamed/cazy-table-GH.txt", stringsAsFactors = FALSE)
cazy.table.hm <- cazy.table.GH[,-1]
rownames(cazy.table.hm) <- cazy.table.GH[,1]
cazy.table.hm[cazy.table.hm == 0] <- NA
cazy.table.hm <- cazy.table.hm[, order(names(cazy.table.hm))]
pdf("Cazy-GH-family-heatmap.pdf", width = width *2, height = height*2)
Heatmap(as.matrix(log10(cazy.table.hm)), name = "log10(Abundance in TPM)", cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, rect_gp = gpar(col = "white", lwd = 1), col = rev(viridis(50)), na_col = "white", column_title = "Samples", column_title_side = "bottom", row_title = "CAZy family", row_names_side = "left", bottom_annotation = ha, row_names_gp = gpar(fontsize = 7), heatmap_legend_param = list(at = c(-1, 0, 1, 2, 3, 4, 5, 6), legend_height = unit(6, "cm"), title_position = "leftcenter-rot", grid_width = unit(1, "cm")) )  
dev.off()
#checking significance
cazy.GH.abundance <- gather(data = cazy.table.GH, key = "Sample", value = "Abundance", 2:49)
cazy.GH.abundance <- merge(x = cazy.GH.abundance, metadata, all.x = TRUE, by.x = "Sample", by.y = "Names")
colnames(cazy.GH.abundance)[2] <- "CazyFamily"
cazy.GH.fraction <- compare_means(formula = Abundance ~ Fraction, data = cazy.GH.abundance, method = "kruskal.test", group.by = "CazyFamily", p.adjust.method = "BH")
cazy.GH.collection <- compare_means(formula = Abundance ~ Collection, data = cazy.GH.abundance, method = "kruskal.test", group.by = "CazyFamily", p.adjust.method = "BH")
cazy.GH.feed <- compare_means(formula = Abundance ~ Feed, data = cazy.GH.abundance, method = "kruskal.test", group.by = "CazyFamily", p.adjust.method = "BH")
#plotting significant results
GH.fraction <- ggplot(data = cazy.GH.abundance %>% filter(CazyFamily %in% as.matrix(na.omit(cazy.GH.fraction[cazy.GH.fraction$p.adj < 0.05, "CazyFamily"]))), mapping = aes(x = CazyFamily, y = log10(Abundance), fill = Fraction)) + scale_fill_manual(values = c("#E47346", "#D49DC7")) + stat_summary(fun.data = mean_sdl, geom = "bar", width = 0.5, position = position_dodge(width = 0.95, preserve = "single")) + coord_flip(expand = c(0,0)) + ylab(label = "log10(Abundance in TPM)") + ylim(0, 3.5) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(plot.margin = unit(c(0,-500,0,-500), "in"))
GH.feed <- ggplot(data = cazy.GH.abundance %>% filter(CazyFamily %in% as.matrix(na.omit(cazy.GH.feed[cazy.GH.feed$p.adj < 0.05, "CazyFamily"]))), mapping = aes(x = CazyFamily, y = log10(Abundance), fill = Feed)) + scale_fill_manual(values = c("#D2A76D", "forestgreen", "#85ABD9")) + stat_summary(fun.data = mean_sdl, geom = "bar", width = 0.5, position = position_dodge(width = 0.95, preserve = "single")) + coord_flip(expand = c(0,0))  + ylab(label = "log10(Abundance in TPM)") + ylim(0, 3.5) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank()) + theme(plot.margin = unit(c(0,-500,0,-500), "in"))
GH.collection <- ggplot(data = cazy.GH.abundance %>% filter(CazyFamily %in% as.matrix(na.omit(cazy.GH.collection[cazy.GH.collection$p.adj < 0.05, "CazyFamily"]))), mapping = aes(x = CazyFamily, y = log10(Abundance), fill = Collection)) + scale_fill_manual(values = c("#75E45E", "#8171DF")) + stat_summary(fun.data = mean_sdl, geom = "bar", width = 0.5, position = position_dodge(width = 0.95, preserve = "single")) + coord_flip(expand = c(0,0)) + ylab(label = "log10(Mean of abundance in TPM)") + theme(axis.title.y = element_blank()) + ylim(0, 3.5)
GH.barplot <- ggarrange(GH.fraction, GH.feed, GH.collection, ncol = 1, nrow = 3, align = "hv", heights = c(15,3,0.75))

ggsave(filename = "GH.barplot.pdf", plot = GH.barplot, device = "pdf", width = width *1.5, height = height*1.5, units = "in")

