#########################################################################
# R script for dada2 v1.4
# ACC deaminase primers - paired end sequencing
# Daniel Manter
# August 23, 2017
#########################################################################


################### Check and load for required packages ###############
.cran_packages <-  c("ggplot2", "gridExtra", "reshape2", "scales", "ggthemes", "devtools", 
                     "seqinr", "Biostrings", "ape", "phangorn")
new.packages <- .cran_packages[!(.cran_packages%in% installed.packages()[,'Package'])]
if (length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org', dependencies=T)

.bioc_packages <- c("phyloseq", "genomes", "msa")
new.packages <- .bioc_packages[!(.bioc_packages %in% installed.packages()[,'Package'])]
if (length(new.packages)) BiocManager::install(new.packages)


.remotes_packages <- c("metagMisc")
new.packages <- .remotes_packages[!(.remotes_packages %in% installed.packages()[,'Package'])]
if (length(new.packages)) remotes::install_github("vmikk/metagMisc")

sapply(c(.cran_packages, .bioc_packages, .remotes_packages), require, character.only=TRUE)
#########################################################################


########################### Check for upgrades ##########################
#biocLite("BiocUpgrade")
#biocLite()
#packageVersion("dada2") # must be version 1.4
#########################################################################


# NOTE: Can skip to line 246 for analysis - previous lines are for processing raw data 


# This section is for the dada2 pipeline and doesn't need to be run again
#########################################################################
if (FALSE) {
  ############################# Get data ##################################
  set.seed(100)
  pathF <- file.path(wdir, "fastq/Forward")
  pathR <- file.path(wdir, "fastq/Reverse")
  cutpathF <- file.path(pathF, "cut")
  cutpathR <- file.path(pathR, "cut")
  filtpathF <- file.path(pathF, "filtered")
  filtpathR <- file.path(pathR, "filtered")
  fastqFs <- sort(list.files(pathF, pattern="fastq"))
  fastqRs <- sort(list.files(pathR, pattern="fastq"))
  
  if (length(fastqFs) != length(fastqRs)) {
    stop("Forward and reverse files do not match.")
  } 
  
  #########################################################################
  
  
  #################### Remove Primers w/ cutadapt ##########################
  # cutadapt must be installed and in your path
  
  for(i in 1:length(fastqFs)) {
    adapter = "-g ^TSTABGCSAARCGBGAVGACTGC"
    infile = file.path(pathF, fastqFs[i])
    outfile = file.path(cutpathF, fastqFs[i])
    cmd <- paste('cutadapt', adapter, '-o', outfile, infile, sep=' ')
    system(cmd)
  }
  
  for(i in 1:length(fastqRs)) {
    adapter = "-g ^GTBACVGMGCASACSACGATRTAG"
    infile = file.path(pathR, fastqRs[i])
    outfile = file.path(cutpathR, fastqRs[i])
    cmd <- paste('cutadapt', adapter, '-o', outfile, infile, sep=' ')
    system(cmd)
  }
  #########################################################################
  
  
  
  ######################### Trim and filter ###############################
  ii <- sample(length(fastqFs), 6)
  for(i in ii) { print(plotQualityProfile(file.path(cutpathF, fastqFs[i])) + ggtitle("Fwd")) }
  
  ii <- sample(length(fastqRs), 6)
  for(i in ii) { print(plotQualityProfile(file.path(cutpathR, fastqRs[i])) + ggtitle("Rev")) }
  
  
  # Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
  filterAndTrim(fwd=file.path(cutpathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                rev=file.path(cutpathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                truncQ=2, maxEE=5, maxN=0, rm.phix=TRUE,
                compress=TRUE, verbose=TRUE, multithread=TRUE)
  #########################################################################
  
  
  ##################### Infer sequence variants ###########################
  filtFs <- list.files(filtpathF, pattern="fastq", full.names=TRUE)
  filtRs <- list.files(filtpathR, pattern="fastq", full.names=TRUE)
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) 
  sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) 
  if(!identical(sample.names, sample.namesR)) {
    stop("Forward and reverse files do not match.")
  }
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  # Learn error rates
  set.seed(100)
  errF <- learnErrors(filtFs, nreads=1e6, multithread=TRUE)
  errR <- learnErrors(filtRs, nreads=1e6, multithread=TRUE)
  
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
  }
  rm(derepF); rm(derepR)
  
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  ##########################################################################
  
  
  ################## remove chimeras & assign taxonomy #####################
  # Remove chimeras
  seqtab <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE)
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  
  # Get taxonomy from MEGAN and Blast
  tax <- read.csv('Megan.Taxonomy.txt', sep="\t")
  row.names(tax) <- tax$ID
  tax <- tax[,2:8]
  tax[tax == ""] <- "unclassified"
  
  for (i in 1:nrow(tax)){
    if (tax[i,1] == "unclassified"){
      tax[i, 1:7] <- "unclassified"
    } else if (tax[i,2] == "unclassified"){
      kingdom <- paste(tax[i,1], "unclassified", sep = " ")
      tax[i, 2:7] <- kingdom
    } else if (tax[i,3] == "unclassified"){
      phylum <- paste(tax[i,2], "unclassified", sep = " ")
      tax[i, 3:7] <- phylum
    } else if (tax[i,4] == "unclassified"){
      class <- paste(tax[i,3], "unclassified", sep = " ")
      tax[i, 4:7] <- class
    } else if (tax[i,5] == "unclassified"){
      order <- paste(tax[i,4], "unclassified", sep = " ")
      tax[i, 5:7] <- order
    } else if (tax[i,6] == "unclassified"){
      family <- paste(tax[i,5], "unclassified", sep = " ")
      tax[i, 6:7] <- family
    } else if (tax[i,7] == "unclassified"){
      genus <- paste(tax[i,6], "unclassified", sep = " ")
      tax[i, 7] <- genus
    }
  }
  
  OTU <- otu_table(seqtab, taxa_are_rows=FALSE)
  colnames(OTU) <- row.names(tax)
  
  # Write to disk
  saveRDS(seqtab, "seqtab_final.rds")
  saveRDS(tax, "tax_final.rds")
  
  # Representative Sequences
  numOtus <- ncol(seqtab)
  uniquesToFasta(seqtab, 
                 ids=paste0("OTU", seq(numOtus)),
                 fout='dada.fasta'  # do not change name
  )
  ##########################################################################
  
  ############################ phyloseq ####################################
  # Make a data.frame holding the sample data
  samples.out <- rownames(seqtab)
  sampdf <- colsplit(string=samples.out, pattern="\\.", names=c("Year", "Site", "Slope", "Strip", "Depth"))
  sampdf$Site[sampdf$Site==1] <- "Sterling"
  sampdf$Site[sampdf$Site==2] <- "Stratton"
  sampdf$Site[sampdf$Site==3] <- "Walsh"
  
  sampdf$Slope[sampdf$Slope==1] <- "Summit"
  sampdf$Slope[sampdf$Slope==2] <- "Side"
  sampdf$Slope[sampdf$Slope==3] <- "Toe"
  
  rownames(sampdf) <- samples.out
  
  # create reference tree
  ref <- readDNAStringSet("dada.fasta")
  ref.aln <- msaMuscle(ref)
  ref.aln <- msaConvert(ref.aln, type="seqinr::alignment")
  D <- dist.alignment(ref.aln, matrix='similarity')
  tre <- nj(D)
  tre <- ladderize(tre)
  plot(tre)
  
  # Construct phyloseq object (straightforward from dada2 outputs)
  ps <- phyloseq(otu_table(OTU, errorIfNULL=FALSE), 
                 sample_data(sampdf), 
                 tax_table(as.matrix(tax)),
                 tre)
  
  saveRDS(ps, 'acdS.ps.rds')
  ##########################################################################
  
  ############################# save files #################################
  dat <- psmelt(ps)
  write.csv(dat, file='otus-with-sample-data.csv')
  
  # create mothur shared file
  numOtus <- nrow(tax_table(ps))
  dat <- otu_table(ps)
  names <- c()
  for (i in seq(1,numOtus)) {
    names[i] <- paste('OTU', i, sep='')
  }
  colnames(dat) <- names
  df <- cbind(label='dada', group=row.names(dat), numOtus=numOtus, dat)
  write.csv(df, file='otu-table.shared', row.names=FALSE)
  
  dat <- tax_table(ps)
  write.csv(dat, file='tax-table.csv')
  
  dat <- colnames(otu_table(ps))
  write.csv(dat, file='rep_seqs.fasta')
  
  save.image("acdS.dada2.rdata")
  ##########################################################################
  
}
#########################################################################


# richness estimates - acdS
##########################################################################
source("myFunctions.R")

ps <- readRDS('acdS.ps.rds')

rich <- estimate_richness(ps, split=TRUE)
write.csv(rich, file='ACC.alpha.diversity.csv', row.names=TRUE)

# Alpha Diversity - graph
p <- plot_richness(ps, x="Site", measures=c("Observed", "Shannon", "InvSimpson"))
gDF <- p$data
gDF$Slope <- factor(gDF$Slope, levels=c('Summit', 'Side', 'Toe'))

colors <- RColorBrewer::brewer.pal(3, 'Set1')
g <- ggplot(gDF, aes(x=Slope, y=value, fill=Slope)) +
  geom_boxplot(alpha=1) + 
  facet_grid(variable~Site, scales='free') + 
  labs(y="Diversity Index") + 
  theme(strip.text.x=element_text(size=10, colour='blue', angle=0)) +
  theme(strip.text.y=element_text(size=10, colour='blue', angle=0)) +
  scale_fill_manual(values=colors)
g

g <- set_panel_size(p, width=unit(3,"in"), height=unit(3,"in"))
ggsave(filename='ACC.Alpha.Diversity.png', device='png', plot=g, units='in', height=8, width=11)
##########################################################################

# Abundance bar plot - classifed data
#########################################################################
P.rel <- transform_sample_counts(ps, function(x) x / sum(x))
rank <- 'Phylum'
P.glom <- tax_glom(P.rel, taxrank=rank, NArm=FALSE)
### write to csv and manually combined with qPCR (not done in R)

gDF <- read.csv('Phylum Total Abundances.csv')
gDF$Slope <- as.factor(gDF$Slope)
gDF$Site <- as.factor(gDF$Site)
gDF$Slope <- factor(gDF$Slope, levels=c('Summit', 'Side', 'Toe'))

colors <- rep(gdocs_pal()(10), times=4)

p <- ggplot(gDF, aes(x=Slope, y=Abundance, fill=Phylum) )
p <- p + geom_bar(position='stack', stat='summary', fun='mean')
p <- p + facet_wrap(~Site)
p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))
p <- p + scale_fill_manual(values=colors)
p <- p + labs(x='', y='Abundance', fill='')
p <- p + guides(fill=guide_legend(ncol=1), color=FALSE)
p <- p + theme(
  legend.text=element_text(size=14),
  legend.position = "right",
  legend.key.size = unit(0.3, "cm"),
  axis.line = element_line(size=0.5, colour = "black"),
  panel.grid.major = element_line(colour = "#d3d3d3"),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(), panel.background = element_blank(),
  axis.title=element_text(size=16),
  axis.text.x=element_text(colour="black", size = 14),
  axis.text.y=element_text(colour="black", size = 14))
p
p <- set_panel_size(p, width=unit(3, 'in'), height=unit(3, 'in'))
ggsave(filename='ACC.Phylum.Total_Abundance.png', w=12, h=5, device='png', plot=p)
#########################################################################


# acdS tree
#########################################################################
library("ggtree")
# read output from DADA2 pipeline
ps <- readRDS('acdS.ps.rds')

P.rel <- transform_sample_counts(ps, function(x) x / sum(x))
rank <- 'Genus'
P.glom <- tax_glom(P.rel, taxrank=rank, NArm=FALSE)

pdf('ACC.tree.pdf', height=6, width=6)
p <- ggtree(phy_tree(P.glom), layout='rectangular', ladderize=F)

dd <- psmelt(P.glom)
data <- merge(p$data, dd, by.x="label", by.y="OTU", rm.na=T)

data$Treatment <- data$Site
data.mean <- aggregate(Abundance ~ Treatment+x+y+node+isTip+angle+Kingdom+Phylum+Class+Order+Family+Genus, data, mean)
data.mean <- as.data.frame(data.mean %>% 
                             arrange(y) %>%
                             group_by(y) %>%
                             mutate(rank=row_number()
                             ))

cols = rep(gdocs_pal()(10), times=5)

spacing <- 0.05
data.mean$xdodge <- data.mean$x + (data.mean$rank * spacing)
data.max <- aggregate(xdodge ~ x+y+node+isTip+angle+Kingdom+Phylum+Class+Order+Family+Genus, data.mean, max)
data.max$Class <- paste('', data.max$Class, sep='')
data.max$Order <- paste('', data.max$Order, sep='')
data.max$Family <- paste('', data.max$Family, sep='')
data.max$Genus <- paste('', data.max$Genus, sep='')
data.max$myLabel <- data.max$Genus

p <- p + geom_tiplab(data=data.max, aes(x=xdodge+0.05, y=y, label=myLabel, color=Phylum), size=3)
p <- p + geom_point(data=data.mean, aes(x=xdodge, y=y, fill=Treatment, size=Abundance), shape=21, na.rm=F) 
p <- p + theme(legend.position="right")
p <- p + scale_color_manual(values=cols)
#p <- p + scale_fill_gdocs()
#p <- p + theme(panel.background=element_rect(fill='gray90', color='gray90'))
#p <- p + theme(legend.background=element_rect(fill='gray90', color='gray90'))
#p <- p + theme(legend.key=element_rect(fill='gray90', color='gray90'))
p <- p + ggplot2::xlim(0, 1.8)
p
dev.off()    
















# Rank-Abundance curve
N <- 30
barplot(sort(taxa_sums(ps), TRUE)[1:N]/nsamples(ps), las=2)

# Network
set.seed(100)
ig <- make_network(ps, dist.fun='bray', max.dist=0.6)
plot_network(ig, ps, color="Site", shape="Slope")

# perMANOVA
metadata <- as(sample_data(ps), 'data.frame')
dist <- phyloseq::distance(ps, "bray")
res <- adonis(dist ~ Site * Slope, data=metadata)
res

# PCoA
ord <- ordinate(ps, method="CAP", dist="bray", ~ Site * Slope)
p <- plot_ordination(ps, ord, color="Site", shape="Slope", title="[Parial] Constrained Analysis of Principal Coordinates")
p <- p + geom_point(size=4)
p <- p + facet_wrap(~ factor(Slope, levels=c("Summit", "Side", "Toe")))
p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))
p <- p + geom_vline(aes(xintercept=0), linetype='dashed')
p <- p + geom_hline(aes(yintercept=0), linetype='dashed')
p <- set_panel_size(p, width=unit(3,"in"), height=unit(3,"in"))
ggsave(filename='PCoA.png', device='png', plot=p, units='in', height=8, width=11)

# PCoA - species scores
df <- plot_ordination(ps, ord, type="species", justDF=TRUE)
mdf <- melt(df[, c("CAP1", "CAP2", "Class", "Order", "Family", "Genus")], 
            id=c("Class", "Order", "Family", "Genus") )
p1 <- ggplot(mdf, aes(Genus, value, color=Order)) + 
  geom_boxplot() + 
  facet_wrap(~variable, 2) + 
  theme( axis.text.x = element_text(angle = -90, vjust = 0.5) )
p1 <- set_panel_size(p1, width=unit(5,"in"), height=unit(3,"in"))
ggsave(filename='PCoA_species.scores.png', device='png', plot=p1, units='in', height=8, width=8)

# Abundance graphs
rank <- 'Phylum'
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
ps.glom <- tax_glom(ps.rel, taxrank=rank, NArm=FALSE)
setwd('C:\\Users\\Daniel.Manter\\Documents\\MyData\\AgroEcosystem\\ACC')
write.csv(file='P.acds.rel.phylum.csv', psmelt(ps.glom))


gDF <- psmelt(ps.glom)
names(gDF)[names(gDF)==rank] <- 'rank'
gDF$rank <- as.character(gDF$rank)
gDF$rank[is.na(gDF$rank)] <- "Unknown"
gDF$rank <- as.factor(gDF$rank)

col_names = c('Site', 'Slope', 'rank')
gDF[col_names] <- lapply(gDF[col_names], factor)
gDF$Slope <- factor(gDF$Slope, levels=c('Summit', 'Side', 'Toe'))

cols <- rep(gdocs_pal()(20), times=3)

p <- ggplot(gDF, aes(x=Site, y=Abundance, fill=rank, color=Slope, shape=OTU) )
p <- p + geom_bar(position='fill', stat='summary', fun.y='mean', color='grey')
p <- p + facet_wrap(~Slope)
p <- p + theme(strip.text.x=element_text(size=10, colour='blue', angle=0))
p <- p + scale_fill_manual(values=cols)
p <- p + labs(x='', y='Abundance', fill='')
p <- p + guides(fill=guide_legend(ncol=3), color=FALSE)
p <- p + theme(
  legend.position = "bottom",
  legend.key.size = unit(0.3, "cm"),
  axis.line = element_line(size=0.5, colour = "black"),
  panel.grid.major = element_line(colour = "#d3d3d3"),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(), panel.background = element_blank(),
  axis.text.x=element_text(colour="black", size = 9),
  axis.text.y=element_text(colour="black", size = 9))
p <- set_panel_size(p, width=unit(3,"in"), height=unit(3,"in"))
ggsave(filename=paste(rank, '_Abundance.png', sep=''), device='png', plot=p, units='in', height=8, width=11)




#######################################################################
P.rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
P.glom <- tax_glom(P.rel, taxrank='Genus', NArm=FALSE)


# plot tree
pdf('ACC.tree.circular.pdf', height=11, width=8)
p <- ggtree(phy_tree(P.glom), layout='circular', ladderize=F)

dd <- psmelt(P.glom)
data <- merge(p$data, dd, by.x="label", by.y="OTU", rm.na=T)

data$Treatment <- data$Site
data.mean <- aggregate(Abundance ~ Treatment+x+y+node+isTip+angle+Kingdom+Phylum+Class+Order+Family+Genus, data, mean)
data.mean <- as.data.frame(data.mean %>% 
                             arrange(y) %>%
                             group_by(y) %>%
                             mutate(rank=row_number()
                             ))

cols = rep(gdocs_pal()(10), times=5)

spacing <- 0.02
data.mean$xdodge <- data.mean$x + (data.mean$rank * spacing)
data.max <- aggregate(xdodge ~ x+y+node+isTip+angle+Kingdom+Phylum+Class+Order+Family+Genus, data.mean, max)
data.max$Class <- paste('', data.max$Class, sep='')
data.max$Order <- paste('', data.max$Order, sep='')
data.max$Family <- paste('', data.max$Family, sep='')
data.max$Genus <- paste('', data.max$Genus, sep='')
data.max$myLabel <- data.max$Genus

p <- p + geom_tiplab(data=data.max, aes(x=xdodge+0.02, y=y, label=myLabel, color=Phylum), size=3)
p <- p + geom_point(data=data.mean, aes(x=xdodge, y=y, fill=Treatment, size=Abundance), shape=21, na.rm=F) 
p <- p + theme(legend.position="right")
p <- p + scale_color_manual(values=cols)
p <- p + scale_fill_gdocs()
p <- p + theme(panel.background=element_rect(fill='gray90', color='gray90'))
p <- p + theme(legend.background=element_rect(fill='gray90', color='gray90'))
p <- p + theme(legend.key=element_rect(fill='gray90', color='gray90'))
p <- p + ggplot2::xlim(0, 1.8)
p
dev.off()    


download.file('ftp:')
