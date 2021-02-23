### Stavros Giannoukakos ###
### METAGENOMICS ANALYSIS - DECONTAM IMPLEMENTATION v.0.1

args <- commandArgs(TRUE)
if (length(args) == 3) {
  # Feature table
  metatable <- args[1]
  # Metadata file
  metadata <- args[2]
  # Output directory
  outdir <- args[3]
} else {
  cat("ERROR - The number of input arguments is not correct...\nEXITING!\n")
  quit()
}

# metatable <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/check_decontam/table.qza"
# metadata <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/check_decontam/clinical_data.txt"
# outdir <- "/Users/stavris/Desktop/Projects/metagenomics_analysis/check_decontam"


options(scipen = 999)
library("decontam")
library("phyloseq")
library("ggplot2")
library("qiime2R")
library("biomformat")



# Importing the metatable
metatableR <- qza_to_phyloseq(features=metatable, metadata=metadata)
# Noting all samples that are Negative Controls
sample_data(metatableR)$is.neg <- sample_data(metatableR)$SampleType == "Negative control"

# Inspect Library Sizes
metatableframe <- as.data.frame(sample_data(metatableR)) # Data frame of the data
# Importing the library size
metatableframe$LibrarySize <- sample_sums(metatableR)
# Sorting
metatableframe <- metatableframe[order(metatableframe$LibrarySize), ]
# Obtaining the index
metatableframe$Index <- seq(nrow(metatableframe))

# Plotting the the library sizes (i.e. the number of reads) for each sample
ggplot(data=metatableframe, aes(x=Index, y=LibrarySize, color=SampleType)) + 
       geom_point() +
       theme_minimal() +
       theme(plot.title = element_text(colour = "#a6a6a4", size=11))+
       labs(x="Index", y="Library Size (number of reads)") +
ggsave(file=paste(outdir,"/library_sizes.png",sep=""), width = 10, height = 6, units = "in", dpi = 1200)

# Decontam function
decontam_dataf <- isContaminant(metatableR, # A feature table recording the observed abundances of each ASV in each sample
                  conc = sample_data(metatableR)$Concentration, # Concentration of amplified DNA in each sample prior to sequencing
                  neg = sample_data(metatableR)$is.neg, # Logical vector where TRUE if sample is a negative control, and FALSE if not (NA entries are not included in the testing)
                  method = "auto", # Frequency, prevalence or combined will be automatically selected based on whether just conc, just neg, or both were provided
                  batch = sample_data(metatableR)$Batch, #  Contaminants identification will be performed independently within each batch
                  batch.combine = "minimum", # For each ASV the probabilities calculated in each batch are combined into a single probability that is compared to 'codethreshold' to classify contaminants 
                  threshold = 0.1,#  The probability threshold 
                  normalize = TRUE, # Normalise
                  detailed = TRUE)

print(  paste("", sep="")    )
# Obtaining the identified contaminant ASVs
microbiome_asvs <- row.names(decontam_dataf[decontam_dataf$contaminant == FALSE, ])
# Obtaining the ASV table from the metatableR 
asv_table <- as.data.frame(otu_table(metatableR))
# Final metatable without the identified contaminant ASVs
final_metatable <- asv_table[row.names(asv_table) %in% microbiome_asvs, ]
# Writing the output filtered feature table in a .biom format
write_biom(make_biom(final_metatable), biom_file=file.path(outdir, "decontam_filtered_table.biom"))
# write.table(final_metatable, file=paste(outdir,"/decontam_filtered_table.tsv", sep=""), sep="\t", row.names=T, col.names=NA, quote=F)

