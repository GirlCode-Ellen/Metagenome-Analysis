

### Author: Ellen Tosolini                        
### Description: Analysis of metagenomic data  of the gut microbiome. Metagenome sequenced from stool 
# samples using whole genome sequencing (shotgun metagenomics). Data analysed is downloaded from the 
# NCBI SRA, accession ID PRJNA686821, study by Wan, Y et al (2021) investigating underdevelopment of 
# gut microbiota and bacteria species as a predictor and biomarker of autism spectrum disorder in children.


## 1. Load Data and required libraries ----------------------------------------------------------------------------------------
abundance_table <- readr :: read_delim("Abun_Wan_Y_ASD.txt") |> as.data.frame()
library(readxl)
metadata <- read_excel("Metadata_Wan_Y_ASD.xlsx") |> as.data.frame()

## Load required libraries for the analysis
library(dplyr)                        # Change and pivot data frames
library(MicrobiomeStat)               # For LINDA differential abundance testing
library(TaxSEA)                       # For taxon set enrichment analysis
library(tidyverse)                    # For data tidying/transformation
library(ggplot2)                      # For making plots to visualise data/results
library(ggpubr)                       # For arranging plots into a grid


## 2. Tidy the metadata table and define the case and control groups ------------------------------------------------
## Make sample id (the "run id") the row names (name each row in the metadata by the corresponding SRR number)
rownames(metadata) <- metadata$Run
##Change name of column 'Sample Name' to more compatible name  (remove requirement for backwards apostrophes)
colnames(metadata)[33] <- "sample_name"
## Extract 't' or 'a' label from sample name column to define two groups
substr(metadata$sample_name,1,1)


## 3. Tidy the abundance table ---------------------------------------------------------------------------------
## check abundance table for duplicated clade names and remove them
apply(abundance_table[,-1],2,sum) 
abundance_table <- abundance_table[!grepl("t__", abundance_table$clade_name),] #remove the type
apply(abundance_table[,-1],2,sum) 


## Change the row and column names of the abundance table
# Make the bacteria species name our row names
species_names <- sapply(abundance_table$clade_name, function(x) strsplit(x, "\\|")[[1]][7]) # extract 7th element - species
species_names <- gsub("s__", "", species_names)  # Remove 's__'
print(species_names)


# Apply species name into our abundance table
abundance_table <- abundance_table |> select(-clade_name) 
abundance_table = as.data.frame(abundance_table)
rownames(abundance_table) <- species_names

# Remove met_output suffix from SRR numbers in column names of abundance table
colnames(abundance_table) <- gsub("_met_output", "", colnames(abundance_table))


## 4. Filter out bacteria species with low abundance --------------------------------------------------------
#Test to ensure each bacteria species is present in 0.5%+ of reads in 3+ samples (greater than 2)
table(apply(abundance_table>0.5,1,sum)>2)
abundance_table = abundance_table[apply(abundance_table>0.5,1,sum)>2,]
#Now table only contains rows with bacteria species that meet above criteria


## 5. Compare metadata and abundance table to ensure both include the same samples
## Remove samples without any phenotype label (denoted in the sample name column with t or a)
metadata <- metadata[!(metadata$sample_name==""),]
metadata$phenotype <- metadata$sample_name
metadata$phenotype <- substr(metadata$phenotype,1,1)

## Make sure both datasets have the same samples
samples2keep <- intersect(metadata$Run, colnames(abundance_table))
abundance_table <- abundance_table[,samples2keep, drop = FALSE] 
metadata <- metadata[metadata$Run %in% samples2keep,]
abundance_table <- abundance_table[,metadata$Run] 


## 6. Test for differential abundance using LINDA
## RUN LINDA on the dataset
linda_res <- linda(abundance_table, metadata, formula = '~phenotype')
linda_results <- linda_res$output$phenotypet

## Extract LINDA ranks
linda_ranks <- linda_results$log2FoldChange
names(linda_ranks) <- rownames(linda_results)

## 7. Plot differential abundance results
metadata$Blautia_massiliensis <- as.numeric(abundance_table["Blautia_massiliensis",])
metadata$Agathobaculum_butyriciproducens <- as.numeric(abundance_table["Agathobaculum_butyriciproducens",])
metadata$Dorea_longicatena <- as.numeric(abundance_table["Dorea_longicatena",])
metadata$Faecalibacterium_prausnitzii = as.numeric(abundance_table["Faecalibacterium_prausnitzii",])
metadata$Group <- metadata$phenotype
metadata$Group <- gsub("t", "Control", metadata$Group)
metadata$Group <- gsub("a", "Autistic", metadata$Group)

fig_p1 <- ggplot(metadata,aes(x=Group,y = Blautia_massiliensis, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme_classic()+
  scale_color_manual(values=c("Autistic"="firebrick4","Control"="deepskyblue4"))+
  stat_compare_means(comparisons = list(c("Autistic","Control")))+
  ylab("B. Massiliensis Rel Abun (%)")+
  guides(color="none")
fig_p1

fig_p2 <- ggplot(metadata,aes(x=Group,y = Agathobaculum_butyriciproducens, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme_classic()+
  scale_color_manual(values=c("Autistic"="firebrick4","Control"="deepskyblue4"))+
  stat_compare_means(comparisons = list(c("Autistic","Control")))+
  ylab("A. butyriciproducens Rel Abun (%)")+
  guides(color="none")
fig_p2

fig_p3 <- ggplot(metadata,aes(x=Group,y = Dorea_longicatena, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme_classic()+
  scale_color_manual(values=c("Autistic"="firebrick4","Control"="deepskyblue4"))+
  stat_compare_means(comparisons = list(c("Autistic","Control")))+
  ylab("D. longicatena Rel Abun (%)")+
  guides(color="none")
fig_p3


fig_p4 <- ggplot(metadata,aes(x=Group,y=Faecalibacterium_prausnitzii,color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("Autistic"="firebrick4","Control"="deepskyblue4"))+
  stat_compare_means(comparisons = list(c("Autistic","Control")))+
  ylab("F. prausnitzii Rel Abun (%)")+
  guides(color="none")
fig_p4

# Arrange the 4 plots into a grid
grid.arrange(fig_p1,fig_p2,fig_p3,
             fig_p4,nrow=2)



## 7. Run taxon set enrichment analysis
## Run TaxSEA for LINDA
linda_taxsea_results <- TaxSEA(taxon_ranks = linda_ranks,lookup_missing = FALSE)
linda_metabo_results <- linda_taxsea_results$Metabolite_producers
linda_disease <- linda_taxsea_results$Health_associations
linda_BugSig <- linda_taxsea_results$BugSigDB
linda_BacDive <- linda_taxsea_results$BacDive_bacterial_physiology
linda_GBM <- linda_taxsea_results$Gut_Brain_Modules_VallesColomer2019
linda_all <- linda_taxsea_results$All_databases

#Inspect results tables manually as well

## Any taxsea results with FDR less than 0.1
if (any(linda_all$FDR < 0.1)) {
  print(linda_all$taxonSetName[linda_all$FDR < 0.1])
} else {
  print("No FDR values are less than 0.1")
}
##None


