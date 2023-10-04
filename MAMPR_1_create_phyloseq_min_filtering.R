###########################################################################
############### CREATE + EDIT PHYLOSEQ OBJECT ##################
###########################################################################


##### Load packages

library(tidyverse)
library(stringr)
library(qiime2R)
library(vegan)
library(geiger)
library(picante)
library(phyloseq) # used to handle the 16S dataset
library(ape) # dependency of phyloseq
library(readr) # for reading in files
library(ggplot2) # plotting tools
library(ggpubr) # plotting tools
library(directlabels) # label plots
library(ggrepel) # label plots
library(dplyr) # manipulate data
library(reshape2) # manipulate data
library(DESeq2) # normalization
library(microbiome)
library(foreach)
library(ggsci)
library(speedyseq)
library(ggtext)


#Build phyloseq object

metadata<-read_tsv("metadata.txt") #read metadata file (tsv)

#Edit metadata file to have full titles for graphing
metadata <- metadata %>%
  mutate(
    Nice_Stage = case_when(stage_score == "veg" ~ "Vegetative",
                           stage_score == "flow_only" ~ "Flowering",
                           stage_score == "mature_siliques" ~ "Ripe Siliques",
                           stage_score == "young_siliques" ~  "Unripe Siliques",
                           stage_score == "Mock" ~  "Mock",
                           stage_score == "NC" ~  "Negative Control"),
    Nice_Tissue = case_when(Tissue == "rose" ~ "Rosette",
                            Tissue == "root" ~ "Root",
                            Tissue == "caul" ~ "Cauline Leaves",
                            Tissue == "silinew" ~  "Siliques",
                            Tissue == "silisen" ~  "Mature Siliques",
                            Tissue == "flow" ~  "Flowers",
                            Tissue == "stem" ~  "Stem",
                            Tissue == "Mock" ~  "Mock",
                            Tissue == "NC" ~  "Negative Control", 
                            Tissue == "soil" ~ "Soil")
  )


## Import taxonomy, otu table, and rooted tree objects from QIIME2 ##

#Taxonomy table, 16S
bactaxonomy <- read_qza("TAX_MAMPR_2018_16S_P1-P14_dada2_nonspike_noec.qza") #Import 16S taxonomy table
bactaxtable<-bactaxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus")) #Separate out by taxon level

#Taxonomy table, ITS
funtaxonomy <- read_qza("TAX_ITS_P1-P14_nonspike.qza") #Import ITS taxonomy table
funtaxtable<-funtaxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus"))

## OTU_table ##
#16S
bacseqvars <- read_qza("TAB_MAMPR_2018_16S_P1-P14_dada2_nonspike_noec.qza")$data #Import OTU table

#ITS
funseqvars <- read_qza("TAB_ITS_P1-P14_nonspike.qza")$data #Import OTU table


## Import tree (rooted) ##
#16S
bactree<-read_qza("SEQ_MAMPR_2018_16S_P1-P14_dada2_nonspike_noec_rooted-tree.qza") #Import tree object

#ITS
funtree<-read_qza("SEQ_ITS_P1-P14_nonspike_rooted-tree.qza") #Import OTU table


#Make the phyloseq object
#16S
STS <-phyloseq::phyloseq(otu_table(bacseqvars, taxa_are_rows = T), 
                         phy_tree(bactree$data), 
                         tax_table(as.data.frame(bactaxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
                         sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sampleid"))
)


#ITS
ITS <-phyloseq::phyloseq(otu_table(funseqvars, taxa_are_rows = T), 
                         phy_tree(funtree$data), 
                         tax_table(as.data.frame(funtaxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
                         sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sampleid"))
)

##### SELECT ITS OR 16S FOR DOWNSTREAM #######

ps_choice <- STS

#################################
##### Edit phyloseq object ######
#################################

#Factorize important metadata#

sample_data(ps_choice)$Nice_Stage <- factor(sample_data(ps_choice)$Nice_Stage, levels = c("Vegetative","Flowering", "Unripe Siliques", "Ripe Siliques", "Mock", "NC"))
sample_data(ps_choice)$Nice_Tissue <- factor(sample_data(ps_choice)$Nice_Tissue, levels=c("Soil", "Root","Rosette","Stem","Cauline Leaves","Flowers","Siliques","Mature Siliques", "Mock", "Negative Control"))
sample_data(ps_choice)$genotype <- factor(sample_data(ps_choice)$genotype, levels=c("col0", "efr", "lyk4", "lore", "fls2", "Mock", "NC", "soil"))
sample_data(ps_choice)$MiSeq_Run <- factor(as.character(as.numeric(sample_data(ps_choice)$MiSeq_Run)))
sample_data(ps_choice)$STS_MiSeq_Run <- factor(as.character(as.numeric(sample_data(ps_choice)$STS_MiSeq_Run)))
sample_data(ps_choice)$ITS_MiSeqRun <- factor(as.character(as.numeric(sample_data(ps_choice)$ITS_MiSeqRun)))
sample_data(ps_choice)$Indiv_ID <- factor(as.character(sample_data(ps_choice)$Indiv_ID))
sample_data(ps_choice)$weeks_old <- factor(as.character(as.numeric(sample_data(ps_choice)$weeks_old)))
sample_data(ps_choice)$plate_id <- factor(sample_data(ps_choice)$plate_id, levels=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14"))
sample_data(ps_choice)$Tissue <- factor(sample_data(ps_choice)$Tissue, levels = c("root", "rose", "stem", "caul", "flow", "silinew", "silisen", "NC", "Mock", "soil"))
sample_data(ps_choice)$stage_score <- factor(sample_data(ps_choice)$stage_score, levels = c("veg", "flow_only", "young_siliques", "mature_siliques"))


# Add combintation of factors categories 
sample_data(ps_choice)$Stage_Tissue <- paste0(as(sample_data(ps_choice)$stage_score, "character"), "_", as(sample_data(ps_choice)$Tissue, "character"))
sample_data(ps_choice)$Stage_Tissue_Genotype <- paste0(as(sample_data(ps_choice)$stage_score, "character"), "_", as(sample_data(ps_choice)$Tissue, "character"), "_", as(sample_data(ps_choice)$genotype, "character"))
sample_data(ps_choice)$Tissue_Genotype <- paste0(as(sample_data(ps_choice)$Tissue, "character"), "_", as(sample_data(ps_choice)$genotype, "character"))


#Simplify ASV names
taxa_names(ps_choice) <- paste0("ASV", seq(ntaxa(ps_choice)))


##Remove samples of the incorrect genotype##

#Conservative - remove all samples with heterozygous contamination and/or without confirmation of homo mut
wp_conservative <- c("wild_type", "hetero_mut","hetero_contam", "no_band")
ps_con <- subset_samples(ps_choice, !mutant_confirmation %in% wp_conservative)
ps_con <-prune_taxa(taxa_sums(ps_con)>0, ps_con)


#Merge samples spread across multiple tubes
derep = merge_samples(ps_con, "Derep_ID") #Abundance values are summed


#Retain dataframe through merge
#OTU table summed in the merge step above, but this alters values in metadata file. 
#IMPORTANT: This step retains data that is the same for samples spread across multiple tubes (stage, harvest date, tissue type etc) but will NOT correctly retain data that is different for tubes, i.e. sample weight

#Create a dereplicated template to retain sample variables
sample_data(ps_con)$Future_rownames <- rownames(sample_data(ps_con)) #Save original rownames column, not used in future dataframe
df <- data.frame(sample_data(ps_con)) #Turn sample data into a dataframe object

#Remove all duplicated rows, but lose numbers
rmreps <- subset(df, !duplicated(Derep_ID)) # Remove rows with duplicate values in Unique_indiv (it will retain 1 of each Unique_indiv values)
sorted <- rmreps[order(rmreps$Derep_ID),] #Order dataframe by Derep_ID
rownames(sorted) <- sorted$Derep_ID 

sample_data(derep) <- sorted #Add retained metadata file to phyloseq object

#Reorient OTU table which is toggled during merge, so that taxa_are_rows = TRUE
otu_table(derep) <- t(otu_table(derep))




#### Amplicon specific filtering  ####

#16S
ps_kingdom <- subset_taxa(derep, Kingdom == "d__Bacteria") #Remove non-bacteria
ps_kingdom  <- subset_taxa(ps_kingdom, Order != " o__Chloroplast") #Remove chloroplasts

#OR

#ITS
ps_kingdom <- subset_taxa(derep, Kingdom == "k__Fungi") #Remove non-fungi

### Minimum filters for dataset ####
nolow <- prune_samples(sample_sums(ps_kingdom)>500, ps_kingdom) #Remove samples with 500 or less reads
norare <- prune_taxa(taxa_sums(nolow)>9, nolow) #Remove ASVs with less than 10 reads overall

exp_only <- subset_samples(norare, Tissue != "Mock" & Tissue != "NC" & Tissue != "silisen") #Retain only experimental samples (all well-sampled plant tissues and soil)
exp_only <- prune_taxa(taxa_sums(exp_only)>0, exp_only) #Remove any ASVs lost by this filtering

plant_only <- subset_samples(norare, Tissue != "Mock" & Tissue != "NC" & Tissue != "soil" & Tissue != "silisen") #Retain only well-sampled plant tissues
plant_only <- prune_taxa(taxa_sums(plant_only)>0, plant_only)  #Remove any ASVs lost by this filtering






