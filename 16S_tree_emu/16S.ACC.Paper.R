### install (if necessary) and load required CRAN or BioConductor packages
##################################################################################
.cran_packages <- c("ape", "ggplot2", "ggthemes", "dplyr", "tidyr", "remotes", "iNEXT",
                    "seqinr", "Biostrings", "ape", "phangorn")
new.packages <- .cran_packages[!(.cran_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

.bioc_packages <- c("phyloseq", "microbiome", "ggtree", "msa")
new.packages <- .bioc_packages[!(.bioc_packages %in% installed.packages()[,'Package'])]
if (length(new.packages)) source('http://bioconductor.org/biocLite.R')
if (length(new.packages)) BiocManager::install(new.packages)

.remotes_packages <- c("metagMisc")
new.packages <- .remotes_packages[!(.remotes_packages %in% installed.packages()[,'Package'])]
if (length(new.packages)) remotes::install_github("vmikk/metagMisc")

sapply(c(.cran_packages, .bioc_packages, .remotes_packages), require, character.only=TRUE)
##################################################################################

# functions
#########################################################################
set_panel_size <- function(p=NULL, g=ggplotGrob(p), file=NULL,
                           margin = unit(1,"mm"),
                           width=unit(4, "cm"),
                           height=unit(4, "cm")){
  
  panels <- grep("panel", g$layout$name)
  panel_index_w<- unique(g$layout$l[panels])
  panel_index_h<- unique(g$layout$t[panels])
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  
  if(getRversion() < "3.3.0"){
    
    # the following conversion is necessary
    # because there is no `[<-`.unit method
    # so promoting to unit.list allows standard list indexing
    g$widths <- grid:::unit.list(g$widths)
    g$heights <- grid:::unit.list(g$heights)
    
    g$widths[panel_index_w] <-  rep(list(width),  nw)
    g$heights[panel_index_h] <- rep(list(height), nh)
    
  } else {
    
    g$widths[panel_index_w] <-  rep(width,  nw)
    g$heights[panel_index_h] <- rep(height, nh)
    
  }
  
  if(!is.null(file))
    ggsave(file, g,
           width = convertWidth(sum(g$widths) + margin,
                                unitTo = "in", valueOnly = TRUE),
           height = convertHeight(sum(g$heights) + margin,
                                  unitTo = "in", valueOnly = TRUE))
  
  invisible(g)
}

print.fixed <- function(x) grid.draw(x)
#########################################################################


# Wrangle data into phyloseq
##################################################################################
# create a phyloseq object with all 16S data
otu <- data.frame(import_biom('AgroEco.tx.dada.biom'), check.names=F)
otu <- otu[,c(1:36)]

taxa <- read.csv('AgroEco.taxonomy.csv')
taxa <- separate(taxa, taxonomy, into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep = ";")
row.names(taxa) <- row.names(otu)
taxa <- taxa[,c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')]

meta <- read.csv('meta_data.csv')
row.names(meta) <- meta$sample_name

ps <- phyloseq(sample_data(meta),
               otu_table(otu, taxa_are_rows=T), 
               tax_table(as.matrix(taxa)))

ps_glom <- tax_glom(ps, taxrank='species')

# get acdS taxa
gibbs <- read.csv('gibbs_subset.tsv', sep='\t')
gibbs.sub <- gibbs[gibbs$sample %in% names(otu),]
acc_taxa <- unique(gibbs.sub[gibbs.sub$function. == "K01505", 'taxon'])
acc_species <- taxa[rownames(taxa) %in% acc_taxa, 'species']

# now need to subset the original phyloseq object with only acdS taxa
ps_sub <- subset_taxa(ps_glom, species %in% acc_species)

# create reference tree
ref <- readDNAStringSet("AgroEco.dada.rep_seqs.fasta")
ref.sub <- ref[acc_taxa]
ref.aln <- msaMuscle(ref.sub)
ref.aln <- msaConvert(ref.aln, type="seqinr::alignment")
D <- dist.alignment(ref.aln, matrix='similarity')
tre <- nj(D)
tre <- ladderize(tre)
plot(tre)

# create final phyloseq object
ps_final <- merge_phyloseq(ps_sub, tre)
ps_rel <- transform_sample_counts(ps_final, function(x) x / sum(x))
#########################################################################


########################################################################
P.glom <- tax_glom(ps_rel, taxrank='genus', NArm=FALSE)

# plot tree
pdf('16S.ACC.circular.tree.pdf', height=10, width=8)
p <- ggtree(phy_tree(P.glom), layout='circular', ladderize=T)

dd <- psmelt(P.glom)
data <- merge(p$data, dd, by.x="label", by.y="OTU", rm.na=T)

data$Treatment <- data$geo_loc_city
data.mean <- aggregate(Abundance ~ Treatment+x+y+node+isTip+angle+kingdom+phylum+class+order+family+genus, data, mean)
data.mean <- as.data.frame(data.mean %>% 
                             arrange(y) %>%
                             group_by(y) %>%
                             mutate(rank=row_number()
                             ))

cols = rep(gdocs_pal()(10), times=5)

spacing <- 0.25
data.mean$xdodge <- data.mean$x + (data.mean$rank * spacing)
data.max <- aggregate(xdodge ~ x+y+node+isTip+angle+kingdom+phylum+class+order+family+genus, data.mean, max)
data.max$Class <- paste('', data.max$class, sep='')
data.max$Order <- paste('', data.max$order, sep='')
data.max$Family <- paste('', data.max$family, sep='')
data.max$Genus <- paste('', data.max$genus, sep='')
data.max$myLabel <- data.max$genus

p <- p + geom_tiplab(data=data.max, aes(x=xdodge+0.05, y=y, label=myLabel, color=phylum), size=3)
p <- p + geom_point(data=data.mean, aes(x=xdodge, y=y, fill=Treatment, size=Abundance), shape=21, na.rm=F) 
p <- p + theme(legend.position="right")
p <- p + scale_color_manual(values=cols)
p <- p + scale_fill_gdocs()
#p <- p + theme(panel.background=element_rect(fill='gray90', color='gray90'))
#p <- p + theme(legend.background=element_rect(fill='gray90', color='gray90'))
#p <- p + theme(legend.key=element_rect(fill='gray90', color='gray90'))
p <- p + ggplot2::xlim(0, 1.8)
p <- p + labs(color="Phylum")
p
dev.off()    
#########################################################################


