### Stavros Giannoukakos ###
### METAGENOMICS - DOWNSTREAM ANALYSIS


args <- commandArgs(TRUE)
if (length(args) == 8) {
  # Input ASV table
  feature.table <- args[1]
  # Taxonomy classification of the ASVs
  taxonomic.classification <- args[2]
  # Phylogenetic tree of sequences
  rooted.tree <- args[3]
  # Sample metadata
  meta.data <- args[4]
  # Output directory where all analysis will be stored
  main.outdir <- args[5]
  # Groups for Differential Abundance testing
  condition_groups <- unlist(strsplit(args[6], ","))
  # Accepted protocol of lib. preparation
  accepted_protocol <- args[7]
  # Number of cores to use for the analysis
  num.cores <- as.numeric(args[8])
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}


# feature.table <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/downstream_blood/table.qza"
# taxonomic.classification <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/downstream_blood/taxonomic_classification.qza"
# rooted.tree <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/downstream_blood/rooted_tree.qza"
# meta.data <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/downstream_blood/clinical_data.txt"
# main.outdir <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/downstream_blood"
# condition_groups <- c("nonCancer", "cancer")
# accepted_protocol <- "standard"
# num.cores <- 1

set.seed(54321)
options(scipen = 999)
options(mc.cores=num.cores) 
suppressPackageStartupMessages({
    library("vegan")
    library("DESeq2")
    library("labdsv")
    library("qiime2R")
    library("phyloseq")
    library("microbiome")
    library("psadd")
    library("gplots")
    library("ggplot2")
    library("cowplot")
    library("pheatmap")
    library("RColorBrewer")
    library("tidyverse")
    library("lme4")
    library("MuMIn")
    library("htmlwidgets")
    library("dplyr")
    library("tsnemicrobiota")
    library("ggpubr")
    library("gridExtra")
    library("agricolae")
    # library("DirichletMultinomial")
    library("reshape2")
    library("magrittr")
})




##### METAGENOMICS - DOWNSTREAM ANALYSIS #####
print("PERFROMING METAGENOMICS DOWNSTREAM ANALSYSIS...")

# outdir <- file.path(main.outdir, "downstream_analysis/post-QC")
# dir.create(outdir, showWarnings = FALSE)
# setwd(outdir)

### Building the Phyloseq Object 
physeq <- qza_to_phyloseq(features=feature.table,
                          tree=rooted.tree,
                          taxonomy=taxonomic.classification,
                          metadata=meta.data)
print(summarize_phyloseq(physeq))

# Remove samples with NA in Condition and exonuclease treated samples
to_remove <- row.names(sample_data(physeq)[is.na(sample_data(physeq)$Condition) | sample_data(physeq)$Protocol != accepted_protocol , ])
physeq <- prune_samples(!(sample_names(physeq) %in% to_remove), physeq)


vegan_otu <- function(physeq) {
             OTU <- otu_table(physeq)
             if (taxa_are_rows(OTU)) {
                 OTU <- t(OTU)
                 }
             return(as(OTU, "matrix"))
             }

### Basic statistics regarding the taxonomic classification
print("Proportion of assignments at each taxonomic level:")
print(apply(tax_table(physeq)[,2:7],2,function(x){1-mean(is.na(x))}))
rm(feature.table, taxonomic.classification, rooted.tree)


#Make a data frame of read depths
data_reads <- data.frame(Reads=sample_sums(physeq))
data_reads$Samples <- rownames(data_reads)  #Add on the sample ID
#Extract the Metadata from the phyloseq object 
metadata <- data.frame(sample_data(physeq))
data_reads <- left_join(data_reads, metadata, "Samples")  #Join on the Metadata
data_reads <- data_reads[order(-data_reads$Reads, data_reads$Cond1), ]  # Order by reads

# Keeping the correct order of the attributes
data_reads$Samples <- factor(data_reads$Samples, levels=data_reads$Samples)
# data_reads$Indication <- factor(data_reads$Indication, levels=unique(data_reads$Indication))




################### GENERAL STATS ###################
#Barplot of Coverage
c1 <-  ggplot(data_reads, aes(x=Samples, y=Reads, fill=Condition)) +
       geom_bar(stat="identity", width=.4, alpha=.8) +
       scale_fill_manual("", values = colorRampPalette(brewer.pal(name="Spectral", n = 11))(length(unique(data_reads$Condition)))) +
       scale_y_continuous(breaks=seq(20000,round(max(data_reads$Reads),1), 20000)) +
       theme_light() +
       theme(axis.text.x = element_text(angle = 90, hjust = 1),
             plot.title = element_text(colour = "#a6a6a4", size=13),
             panel.grid.major.x = element_blank(),
             # panel.border = element_blank(),
             panel.border =  element_rect(colour = "gray68", size=.2, fill=NA),
             axis.ticks.y = element_blank(),
             panel.grid.minor = element_blank(),
             axis.text.y = element_text(size=8.5, color="gray45"),
             panel.grid.major.y = element_line(size=.1, color="gray78")) +
       labs(title="Coverage overview", x="Samples", y="Reads") +
       if (length(unique(data_reads$Batch > 1))) {facet_grid(. ~ Batch, scales='free', space = "free", shrink=F)}
ggsave(c1, file=paste(main.outdir,"/1.coverage-barplot.perSample.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)

# Boxplot of Per group coverage 
c2 <-  ggplot(data_reads, aes(x=Cond1, y=Reads, fill=Cond1)) +
       geom_point(aes(colour = factor(Batch)), alpha = 1, position = "jitter", shape=18, size=1.5) +
       geom_boxplot(aes(fill=Cond1), alpha = 0.8, width=0.4) +
       theme_minimal() +
       theme(plot.title = element_text(colour = "gray20", size=11),
             panel.border = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major.x = element_line(size=.1, color="gray78"),
             panel.grid.major.y = element_line(size=.1, color="gray78"),
             legend.title = element_blank(),
             # axis.text.x = element_text(angle = 0, hjust = 1),
             legend.text=element_text(size=11),
             legend.position = "right") +
       labs(title = "Per group boxplot reveiling its group's coverage",
            y = "Groups",
            x = "Reads")
ggsave(c2, file=paste(main.outdir,"/2.coverage-boxplot.perGroup.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(c1, c2)


################### TAXONOMIC ASSIGNEMENT ###################
# Output KRONA visualisation of the classification for each sample
plot_krona(physeq, output=paste(main.outdir,"/3.krona-samples",sep=""), variable="Samples", trim=T)
# Output overall KRONA visualisation of the classification
plot_krona(physeq, output=paste(main.outdir,"/4.krona-cond1",sep=""), variable="Cond1", trim=T)

# source("/Users/stavris/Desktop/Projects/metagenomics_analysis/Ranalysis/plot_hierarchy.R")
# saveWidget(plot_hierarchy(physeq = physeq, number_of_taxa_to_plot = 40), file=paste(main.outdir,"/sankey.diagram.html",sep=""))
# saveWidget(plot_hierarchy(physeq = physeq, number_of_taxa_to_plot = 50, plot_type = "sunburst"), file=paste(main.outdir,"/sunburst.diagram.html",sep=""))


### Composition plots
ps_rel_abund <- phyloseq::transform_sample_counts(physeq, function(x){x/sum(x)})
phyl <- phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
        geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
        facet_wrap(~ Cond2, scales="free_x", shrink=F) + 
        theme(panel.background = element_blank(),
              # axis.text.x = element_line(size=.1, color="gray78"),
              panel.grid.major.y = element_line(size=.3, color="gray50"),
              axis.ticks.x = element_blank()) +
        guides(fill=guide_legend(ncol=1)) +
        labs(title="Microbiome relative composition in Phylumlevel",
            x = "", 
            y = "Relative Abundance\n")
ggsave(phyl, file=paste(main.outdir, "/5.taxonomy.phylum.rel-composition.png", sep=""), width = 12, height = 8, units = "in", dpi = 1200)
rm(ps_rel_abund, phyl)


#Agglomerate to phylum-level and rename
# ps_phylum <- phyloseq::tax_glom(physeq, "Phylum")
# phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
# phyloseq::psmelt(ps_phylum) %>% 
#           ggplot(data = ., aes(x = Cond2, y = Abundance)) +
#           geom_boxplot(outlier.shape  = NA) +
#           geom_jitter(aes(color = OTU), height = 0, width = .2) +
#           guides(fill=guide_legend(ncol=1)) +
#           labs(x = "", y = "Abundance\n") +
#           facet_wrap(~ OTU, scales = "free")
#           


################### PHYLOGENETICS ###################
### Rarefraction Curves
# Rarefaction curves allow us to assess whether we've likely discovered all the microbial 'species' / ASVs in a sample.
# The lower the true species diversity, and the higher the sampling depth (number of reads), the more likely we are to 
# have saturated our species discovery curves. We can plot per-sample rates of species discovery using the 'rarecurve'
# function in vegan. By default this will plot a different line for each sample, so if you have hundreds of samples,
# these graphs can look very cluttered!
png(paste(main.outdir,"/6.rarefraction_curves.png",sep=""), width = 12, height = 6, units = "in", res = 1200)
rarecurve(t(otu_table(physeq)), col="#11698e", step=50, cex=0.37,
          main="Rarefaction Curves",
          ylab="Observed No. of Species",
          xlab="Reads")
dev.off()



################## ALPHA DIVERSITY ANALYSIS ##################
ps_rare <- rarefy_even_depth(physeq, min(data_reads$Reads), rngseed = 54321)  # Rarefying Samples
ps_rare_df <- as(sample_data(ps_rare),"data.frame")
asvs_mat <- vegan_otu(ps_rare)


alpha_div <- c("Shannon","Simpson","InvSimpson","Fisher")
microbiome_richness <- estimate_richness(ps_rare, measures=alpha_div)
# Adding the metadata to the alpha diversity metrics
microbiome_richness$Samples <- rownames(microbiome_richness)
microbiome_richness <- left_join(microbiome_richness, metadata, "Samples")


i <- 0
p <- list()
q <- list()
for (measure in alpha_div) {  print(measure)
    i <- i+1
    ### Plotting in regard of the Indication
    p1 <- ggplot(microbiome_richness %>% arrange(Cond1), aes(x=Samples, y=eval(as.name(measure)))) +
           geom_point(size=3, pch=22, aes(fill=Cond1)) + 
           scale_fill_brewer(palette="Set3") +
           theme(axis.text.x = element_text(angle = 90, hjust = 1),
                 plot.title = element_text(colour = "#a6a6a4", size=13),
                 panel.border = element_blank(),
                 panel.background = element_rect(fill="gray98"),
                 panel.grid.major.x = element_line(size=.1, color="gray78"),
                 panel.grid.major.y = element_line(size=.1, color="gray78"),
                 panel.grid.minor = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.y = element_text(size=8.5, color="gray45"),
                 legend.key = element_rect(colour = "transparent", fill = "white"),
                 legend.position="none") +
           facet_grid(. ~ Cond1, scales='free', space = "free", shrink=F) +
           # facet_wrap(.~Smoking_status) +
           labs(title=paste(measure,"Diversity",sep=" "), x="Samples",y=paste("Alpha Diversity (",measure,")",sep=""))
    # ggsave(p2, file=paste(main.outdir, paste("/6.alpha-diversity.", measure, ".cond1.squareplot.png", sep=""),sep=""), width = 12, height = 6, units = "in", dpi = 1200)
    p[[i]] <- ggplotGrob(p1)
 
    p2 <- ggplot(microbiome_richness, aes(x=Cond1, y=eval(as.name(measure)))) +
        geom_point(aes(colour = factor(Batch)), alpha = 1, position = "jitter", shape=18, size=1.5) +
        geom_boxplot(aes(fill=Cond1), alpha = 0.8, width=0.4, show.legend=FALSE) +
        scale_fill_brewer(palette="Set3") +
        theme(plot.title = element_text(colour = "#a6a6a4", size=13),
              #axis.text.x = element_text(angle = 35, hjust = 1),
              panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size=.1, color="gray78"),
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_text(size=8.5, color="gray45"),
              legend.position="bottom") +
        labs(title=paste(measure,"Diversity", sep=" "), x="Samples", y=paste("Alpha Diversity (",measure,")", sep=""))
    # ggsave(p2, file=paste(main.outdir, paste("/7.alpha-diversity.", measure, ".cond1.boxplot.png", sep=""),sep=""), width = 12, height = 6, units = "in", dpi = 1200)
    q[[i]] <- ggplotGrob(p2)
}

ggarrange(plotlist=p)
ggsave(file=paste(main.outdir, "/7.alpha-diversity.cond1.squareplot.png", sep=""), width = 16, height = 10, units = "in", dpi = 1200)
ggarrange(plotlist=q)
ggsave(file=paste(main.outdir, paste("/8.alpha-diversity.cond1.boxplot.png", sep=""),sep=""), width = 16, height = 10, units = "in", dpi = 1200)


micr <- microbiome_richness; row.names(micr) <- micr$Samples
d <- micr %>%  gather(key = metric, value = value, c("Shannon","Simpson","Fisher")) %>%
     mutate(metric = factor(metric, levels = c("Shannon","Simpson","Fisher"))) %>%
     ggplot(aes(x = Cond2, y = value)) +
     geom_boxplot(outlier.color = NA) +
     geom_jitter(aes(color = Cond2), height = 0, width = .2) +
     geom_text(aes(label=Samples), size=1.5, hjust=.2, vjust=-.8) +
     labs(x = "", y = "") +
     facet_wrap(~ metric, scales = "free") +
     theme(panel.border = element_blank(),
           panel.background = element_rect(fill="gray98"),
           panel.grid.major.x = element_line(size=.1, color="gray78"),
           panel.grid.major.y = element_line(size=.1, color="gray78"),
           panel.grid.minor = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text.y = element_text(size=8.5, color="gray45"),
           legend.key = element_rect(colour = "transparent", fill = "white"),
           legend.position="none")
ggsave(d, file=paste(main.outdir, paste("/9.alpha-diversity.cond2.boxplot.png", sep=""),sep=""), width = 12, height = 8, units = "in", dpi = 1200)


j <- 0
p <- list()
### Wilcoxon Correlation testing
for (measure in alpha_div) { print(paste("Performing Wilcoxon test on the data based on ", measure, " distance", sep=""))
    j <- j+1
    p[[j]] <- wilcox.test(eval(as.name(measure)) ~ Cond2, data = microbiome_richness, exact = FALSE, conf.int = TRUE)
}
rm(alpha_div,microbiome_richness,i,p,q,p1,p2,to_remove,measure,meta.data,num.cores,micr,d, j)



################## BETA DIVERSITY ANALYSIS - Community Statistics ################## 
# Looking at metrics of community structure - differences among samples in the 
# distribution and relative abundance of our microbial taxa
### Plot HeatMap
ps_rare_top50 <- prune_taxa(names(sort(taxa_sums(ps_rare),TRUE)[1:50]), ps_rare)
# Generate Metadata for Plotting Colours
beta.data <- data.frame(Condition=sample_data(ps_rare_top50)$Cond3)
rownames(beta.data) <- rownames(sample_data(ps_rare_top50))
# Plotting the heatmap
h1 <- pheatmap(otu_table(ps_rare_top50), 
         cluster_cols = T, 
         scale="row", 
         fontsize = 8,
         border_color = NA,
         angle_col = 90,
         annotation_col = beta.data)
ggsave(h1, file=paste(main.outdir, "/10.beta-diversity.heatmap.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
rm(ps_rare_top50, beta.data, h1)



### Ordination analysis
# Bray-Curtis
ord.nmds.bray <- ordinate(ps_rare, method="NMDS", k=2, distance='bray', trymax=50)
p1 <- plot_ordination(ps_rare, ord.nmds.bray, color="Cond4", title="NMDS plot (Bray-Curtis distance)") +
                      geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                      theme_minimal()
p2 <- plot_ordination(ps_rare, ord.nmds.bray, type = "taxa", color = "Phylum") +
                      geom_point(size=1, alpha=.5) + 
                      theme_minimal()
# ggsave(p1, file=paste(main.outdir,"/9.ordination.bray-curtis.NMDS.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
pp <- ggarrange(p1, p2, ncol = 1, nrow = 2)
ggsave(pp, file=paste(main.outdir,"/11.ordination.nmds.bray-curtis.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(ord.nmds.bray, p1, p2, pp)


# Unweighted Unifrac
ord_nmds_unifrac_unweighted <- ordinate(ps_rare, method="NMDS", distance="unifrac", weighted=FALSE)
g1 <- plot_ordination(ps_rare, ord_nmds_unifrac_unweighted, color="Cond4", title="NMDS plot (Unweighted Unifrac)") + 
                      geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                      theme_minimal()
g2 <- plot_ordination(ps_rare, ord_nmds_unifrac_unweighted, type = "taxa", color = "Phylum") +
                      geom_point(size=1, alpha=.5) + 
                      theme_minimal()
# ggsave(g1, file=paste(main.outdir,"/10.ordination.unweighted-unifrac.NMDS.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
gg <- ggarrange(g1, g2, ncol = 1, nrow = 2)
ggsave(gg, file=paste(main.outdir,"/12.ordination.nmds.unweighted-unifrac.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(ord_nmds_unifrac_unweighted, g1, g2, gg)


# Weighted Unifrac
ord_nmds_unifrac_weighted <- ordinate(ps_rare, method="NMDS", distance="unifrac", weighted=TRUE)
q1 <- plot_ordination(ps_rare, ord_nmds_unifrac_weighted, color="Cond4", title="NMDS plot (Weighted Unifrac)") +
                      geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                      # stat_ellipse(geom="polygon",aes(fill=Cond2),type="norm",alpha=0.1) + 
                      theme_minimal()
q2 <- plot_ordination(ps_rare, ord_nmds_unifrac_weighted, type = "taxa", color = "Phylum") +
                      geom_point(size=1, alpha=.5) +
                      theme_minimal()
# ggsave(q1, file=paste(main.outdir,"/11.ordination.weighted-unifrac.NMDS.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
qq <- ggarrange(q1, q2, ncol = 1, nrow = 2)
ggsave(qq, file=paste(main.outdir,"/13.ordination.nmds.weighted-unifrac.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(ord_nmds_unifrac_weighted, q1, q2, qq)


# Canonical correspondence analysis
ord_cca <- ordinate(ps_rare, "CCA")
u1 <- plot_ordination(ps_rare, ord_cca, type = "samples", color = "Cond4", title="CCA plot (Bray-Curtis)") + 
    geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
    theme_minimal()
u2 <- plot_ordination(ps_rare, ord_cca, type = "taxa", color = "Phylum") + 
    geom_point(size=3, alpha=.5) + 
    theme_minimal()
uu <- ggarrange(u1, u2, ncol = 1, nrow = 2)
ggsave(uu, file=paste(main.outdir,"/14.ordination.cca.bray-curtis.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(ord_cca, u1, u2, uu)




# T-SNE in Weighted Unifrac
tsne_res_bc <- tsne_phyloseq(ps_rare, distance='bray', perplexity = 8, verbose=0, rng_seed = 3901)
# Plot the results.
t1 <- plot_tsne_phyloseq(ps_rare, tsne_res_bc, color = 'Cond4', title='t-SNE (Bray-Curtis)') +
                         geom_point(size=3) +
                         geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                         theme_minimal() 
tsne_res_jsd <- tsne_phyloseq(ps_rare, distance='jsd', perplexity = 8, verbose=0, rng_seed = 3901)
t2 <- plot_tsne_phyloseq(ps_rare, tsne_res_jsd, color = 'Cond4', title='t-SNE (Jensen–Shannon divergence)') +
                         geom_point(size=3) +
                         geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                         theme_minimal() 

tsne_res_uu <- tsne_phyloseq(ps_rare, distance='unifrac', perplexity = 8, verbose=0, rng_seed = 3901)
t3 <- plot_tsne_phyloseq(ps_rare, tsne_res_uu, color = 'Cond4', title='t-SNE (Unweighted UniFrac)') +
                         geom_point(size=3) +
                         geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                         theme_minimal() 

tsne_res_wu <- tsne_phyloseq(ps_rare, distance='wunifrac', perplexity = 8, verbose=0, rng_seed = 3901)
t4 <- plot_tsne_phyloseq(ps_rare, tsne_res_wu, color = 'Cond4', title='t-SNE (Weighted UniFrac)') +
                         geom_point(size=3) +
                         geom_text(label=sample_data(ps_rare)$Samples, size=2, hjust=.5, vjust=-1.1) +
                         theme_minimal() 
# ggsave(t1, file=paste(main.outdir,"/13.ordination.tsne.weighted-unifrac.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)

tt <- ggarrange(t1, t2, t3, t4, ncol = 2, nrow = 2, common.legend = T, legend="bottom")
ggsave(tt, file=paste(main.outdir,"/15.ordination.tsne.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(tsne_res_bc, tsne_res_jsd, tsne_res_uu, tsne_res_wu, t1, t2, t3, t4, tt)



### Within-Group Divergence
# popvals <- as.character(unique(sample_data(ps_rare)$Cond2))  #Extract Population Info  
# div.list <- rep(list(NA), length(popvals))  ## Make an Empty List to Store Divergence Values
# reference <- apply(subset_samples(ps_rare, Cond2 == 'nonCancer'), 1, median)
# 
# 
# for(k in 1:length(popvals)){
#     div.list[[k]] <- divergence(subset_samples(ps_rare, Cond2==popvals[k]), reference)
# }
# # Expand Divergence Values from List to Data Frame 
# divergence_values <- data.frame(divergence=unlist(div.list), Sample=names(unlist(div.list)))
# #Add On Divergence Values    
# data_divergence<-left_join(ps_rare_df, divergence_values, "Samples")
# 
# #Plot Divergence By Pop 
# ggplot(data_divergence, aes(x=Indication_I, y=divergence)) + 
#        geom_violin(aes(fill=Indication_I)) + 
#        theme_minimal()
# ggsave(file=paste(main.outdir,"/10.divergence.indication_I.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)
# rm(popvals, div.list, divergence_values, data_divergence)


# Statistical Testing Using PERMANOVA
# indication_III.adonis <- adonis(asvs_mat ~ factor(Indication_III)+Age, data=ps_rare_df, permutations=10000, method="bray")
# print(soil.adonis)





################## Network Analysis ################## 
ig <- make_network(ps_rare, "taxa", distance = "jaccard", max.dist = 0.95)
p1 <- plot_network(ig, ps_rare, type="taxa", point_size = 5, label=NULL, color="Class", line_alpha = 0.05)
# ggsave(p1, file=paste(main.outdir,"/14.network-analysis.taxa.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)

ig2 <- make_network(ps_rare, "samples", distance = "jaccard", max.dist = 0.95)
p2 <- plot_network(ig2, ps_rare, type="samples", point_size = 5, label=NULL, color="Cond1", line_alpha = 0.05)
# ggsave(p2, file=paste(main.outdir,"/14.network-analysis.samples.indIII.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)

p3 <- grid.arrange(p1, p2, ncol = 1, nrow = 2)
ggsave(p3, file=paste(main.outdir,"/16.network-analysis.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(ig, ig2, p1, p2, p3)


################## Indicator Analysis and Differential Abundance Testing ################## 
### Indicator Analysis
# The ideal indicator species: 
# - is exclusively faithful to a particular group (EXCLUSIVITY) 
# - occurs in all sample units within a group(FIDELITY)

# This code will convert the Indication II to numeric clusters where cancer = 1 and NO cancer = 2
cancer_numeric <- with(sample_data(ps_rare), ifelse(Cond2=="cancer", 1, 2))
# Run the indicator analysis
cancer_indic <- indval(asvs_mat, cancer_numeric)
#Extract The Full Dataframe
cancer_indic_output <- data.frame(indval=cancer_indic$indcls, pval=cancer_indic$pval, cluster=cancer_indic$maxcls)
#Add in The Taxon Labels as a column 
#R adds an annoying 'X' to the leading edge of hexadecimal IDs that start with numbers, so we need to strip those out using 'gsub'
cancer_indic_output$taxon <- gsub("^X", "", rownames(cancer_indic_output))
#Subset to Indicator Values >0.5
cancer_indic_significant <- subset(cancer_indic_output, indval>=0.5 & pval<=0.05)
print(cancer_indic_significant)
#Annotate With Taxonomy 
ps_rare_tax <- data.frame(tax_table(ps_rare))
ps_rare_tax$taxon <- rownames(ps_rare_tax)
cancer_indic_significant <- left_join(cancer_indic_significant, ps_rare_tax, "taxon")
print(cancer_indic_significant)
rm(cancer_numeric, cancer_indic, cancer_indic_output, ps_rare_tax)


### Differential Abundance Testing Using DESeq2
# Function To Handle Results Objects and Annotate with Taxonomy
taxo<-function(resultsobject, physeqobject, alpha){
               sigtab <- resultsobject[which(resultsobject$padj<alpha), ]
               sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeqobject)[rownames(sigtab), ], "matrix"))
               colnames(sigtab)[7:12] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
               return(sigtab)
}   

# Function To Calculate Geometric Means
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
} 

# Function to parse significant data from DESeq2 results
deseqplot_data <- function(sigtab){
                  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))  # Phylum order
                  x = sort(x, TRUE)
                  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
                  # Copy Across Genus Labels and Fill in Any Unassigned
                  sigtab$Genus.long<-as.character(sigtab$Genus)
                  sigtab$Genus.long[grep("unclassified",sigtab$Genus)] <-paste0("[",as.character(sigtab$Family[grep("unclassified",sigtab$Genus)]),"]")
                  # Genus order
                  x = tapply(sigtab$log2FoldChange, sigtab$Genus.long, function(x) max(x))
                  x = sort(x, TRUE)
                  sigtab$Genus.long = factor(as.character(sigtab$Genus.long), levels=names(x))
                  return(sigtab)
}

#Fit The Model We Are Interested in Using The function in phyloseq - gets the data in a format DESeq can read
condmod <- phyloseq_to_deseq2(physeq, ~ Cond2)
# Calculate Geometric mean Counts Ourselves   
cond_ge <- apply(DESeq2::counts(condmod), 1, gm_mean)
cond_size <- estimateSizeFactors(condmod, geoMeans=cond_ge)

#Fit the DESeq model   
cond.deseq = DESeq(cond_size, test="Wald", fitType="mean")

#Extract The Results and Specify the Contrasts We Want - Here 'No Vegetation' will be used as the reference level 
cond_results<-results(cond.deseq, contrast = c("Cond2", condition_groups[2], condition_groups[1]))

#Annotate The Taxa with the taxonomy   
cond_results_taxa <- taxo(cond_results, physeq, 0.05)

#### Store The Data in a Format that We Can Plot With GGplot
cond_plot_data <- deseqplot_data(cond_results_taxa)

#####Generate Plot
dep <- ggplot(cond_plot_data, aes(x=Genus.long, y=log2FoldChange)) + 
       geom_point(shape=21, size=6, aes(fill=Phylum)) + 
       theme_classic() + 
       theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +  
       scale_fill_brewer(palette="Set2") + 
       geom_hline(yintercept = 0,linetype="dashed") + labs(x="Genus") + 
       scale_shape_manual(values=c(21, 22, 23, 24, 25))
ggsave(dep, file=paste(main.outdir,"/17.differential-abundances.analysis.png",sep=""), width = 10, height = 8, units = "in", dpi = 1200)
rm(taxo, gm_mean, deseqplot_data, condmod, cond_ge, cond_size, cond.deseq, cond_results, cond_results_taxa, cond_plot_data, dep)
