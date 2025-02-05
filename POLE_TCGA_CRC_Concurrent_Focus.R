# Project: POLE-mutations in colorectal cancer
# TCGA as external validation cohort, simulation of small NGS panel (Illumina Focus Panel)
# For Illumina Focus Panel: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/ampliseq-focus-panel.html
# Author: Nic G. Reitsam
# NOTE: This R script provides a **general framework** and **guidance** for how we accessed and analyzed TCGA data for this work.
# It is **not identical** to the exact code used for generating all the figures in the publication. Certain aspects, like specific file paths, data filtering, and processing, may have been adjusted for presentation or reproducibility purposes. The code herein should help others replicate the analysis or modify it for their own use cases.
# E.g. for our own cohort, we used a customized version of somaticInteractions in analogy to the function used here. As input, we did not use a MAF but a matrix.
# For the exact steps and configurations used to generate the figures in the publication, please contact the authors.


#install necessary packages TCGAbiolinks, maftools, dplyr
#For citation: Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, Sabedot T, Malta TM, Pagnotta SM, Castiglioni I, Ceccarelli M, Bontempi G, Noushmehr H (2015). “TCGAbiolinks: An R/Bioconductor package for integrative analysis of TCGA data.” Nucleic Acids Research. doi:10.1093/nar/gkv1507, http://doi.org/10.1093/nar/gkv1507.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "maftools", "dplyr"))

library(TCGAbiolinks)
library(maftools)
library(dplyr)

#Query Colonic and Rectal Adenocarcinomas
query_coad <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

query_read <- GDCquery(
  project = "TCGA-READ",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query_coad)
GDCdownload(query_read)

mut_coad <- GDCprepare(query_coad)
mut_read <- GDCprepare(query_read)

combined_mutation_CRC <- rbind(mut_coad, mut_read)

maf_data <- read.maf(maf = combined_mutation_CRC)

#only POLE-mutated CRCs in TCGA
pole_mutated <- subsetMaf(maf = maf_data, genes = "POLE")
pole_sample_ids <- unique(pole_mutated@clinical.data$Tumor_Sample_Barcode)

#For simulating a small NGS-based panel, we need the gene list (Illumina Foucs Panel)
#gene list of focus panel illumina
genes_to_check_focus <- c("ABL1", "AKT1", "AKT3", "ALK", "AR", "AXL", "BRAF", 
                    "CCND1", "CDK4", "CDK6", "CTNNB1", "DDR2", "EGFR", 
                    "ERBB2", "ERBB3", "ERBB4", "ERG", "ESR1", "ETV1", 
                    "ETV4", "ETV5", "FGFR1", "FGFR2", "FGFR3", "FGFR4", 
                    "GNA11", "GNAQ", "HRAS", "IDH1", "IDH2", "JAK1", 
                    "JAK2", "JAK3", "KIT", "KRAS", "MAP2K1", "MAP2K2", 
                    "MET", "MTOR", "MYC", "MYCN", "NRAS", "NTRK1", 
                    "NTRK2", "NTRK3", "PDGFRA", "PIK3CA", "PPARG", 
                    "RAF1", "RET", "ROS1", "SMO")

#for comparability with our own cohort, we also consider off target calls
#here are all genes that were called in our own cohort
genes_to_check2_own_cohort <- c("ERBB2", "KRAS", "PIK3CA", "AR", "IDH1", "NRAS", "APC", "BRAF", "EGFR", "IDH2", 
                     "MAP2K1", "MET", "RET", "FGFR3", "ERBB3", "KIT", "PDGFRA", "MYC", "MTOR", "CTNNB1", 
                     "GNA11", "ALK", "AKT1", "ROS1", "HRAS", "RAF1", "MED12", "DCUN1D1", "FGFR2", "BRCA1", 
                     "TP53", "ATM", "CDH1", "KDR", "PTEN", "SMAD4", "MAP2K2", "ESR1", "GNAQ", "JAK2", 
                     "DDR2", "JAK3", "FGFR4", "JAK1", "FGFR1", "NF1", "BRCA2")

#we now take the union to get all genes (all in the panel, and all that were indeed called) and end up with having 64 genes
unique_genes <- union(genes_to_check_focus, genes_to_check2_own_cohort)
num_unique <- length(unique_genes)
print(paste("Number of unique genes:", num_unique))

# Filter for mutations in POLE-mutated samples for gene list; concurrent_mutations = all, concurrent_mutations_focus = only gene list of Focus Panel
concurrent_mutations <- subsetMaf(maf = maf_data, tsb = pole_sample_ids, genes = unique_genes)
concurrent_mutations_focus <- subsetMaf(maf = maf_data, tsb = pole_sample_ids, genes = genes_to_check_focus)

oncoplot(maf = concurrent_mutations, genes = unique_genes, top = 80)
oncoplot(maf = concurrent_mutations_focus, genes = genes_to_check_focus, top = 64)

#Check how many samples in TCGA have no concurrent mutation in genes of Focus Panel (n=2)
#If we consider also the whole gene list with all genes called in our own cohort using Focus Panel, this number drops to zero (n=0)
concurrent_samples <- unique(concurrent_mutations@data$Tumor_Sample_Barcode)
concurrent_samples_focus  <- unique(concurrent_mutations_focus@data$Tumor_Sample_Barcode)
missing_samples <- setdiff(pole_sample_ids, concurrent_samples)
missing_samples_focus <- setdiff(pole_sample_ids, concurrent_samples_focus)
print(missing_samples)
print(missing_samples_focus)

#we now just want to look at certain SNVs (also for comparability)
snv_types <- c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Silent")

# Subset maf_data object based on the Variant_Classification column
snv_data<- subsetMaf(maf_data, query = "Variant_Classification %in% snv_types")
table(snv_data@data$Variant_Classification)

#only consider the pre-defined SNVs
concurrent_SNVs_all <- subsetMaf(concurrent_mutations, query = "Variant_Classification %in% snv_types")
concurrent_SNVs_focus <- subsetMaf(concurrent_mutations_focus, query = "Variant_Classification %in% snv_types")

#setting colors: green/turquoise for mutations
mutation_colors <- c(
  'Missense_Mutation' = '#0c7068',
  'Nonsense_Mutation' = '#0c7068',
  'Splice_Site' = '#0c7068',
  'Silent' = '#0c7068',
  'Frame_Shift_Del' = '#0c7068',
  'Frame_Shift_Ins' = '#0c7068'
)

#comment: after filtering for our pre-defined SNVs (unique_genes, all), only one POLE-mutated CRC in TCGA does not show a concurrent alteration

#oncoplot with modified colors and background color
oncoplot(maf = concurrent_SNVs_all, 
         genes = unique_genes,
         top = 64,  
         colors = mutation_colors,
         fontSize = 0.5,
         bgCol = "#f5f4f4")

oncoplot(maf = concurrent_SNVs_focus, 
         genes = genes_to_check_focus,
         top = 64,  
         colors = mutation_colors,
         fontSize = 0.5,
         bgCol = "#f5f4f4") 

somaticInteractionsAllSNVs <- somaticInteractions(maf = concurrent_SNVs_all, fontSize = 0.5, showSum = FALSE, top=70)
somaticInteractionFocusSNVs <- somaticInteractions(maf = concurrent_SNVs_focus, fontSize = 0.5, showSum = FALSE, top=70)

sessionInfo()
