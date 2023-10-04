##################################################
############### BETA DIVERSITY ##################
##################################################

#Assumes loaded packages and phyloseq objects from MAMPR_1 script
#phyloseq objects:
ps_con
exp_only
plant_only



####################################################
#### Create absolute abundance phyloseq object #####
####################################################

#The merge step to generate the phyloseq objects exp_only and plant_only messes up the quantitation for the spike read ratios (stored in metadata)
#Use unmerged ps_con

#Merge samples for absolute abundance
aa <- ps_con

#Add percentage spikes present in each sample
sample_data(aa)$ITS_percent_spikes <- sample_data(aa)$ITS_Spike/(sample_data(aa)$ITS_Spike + sample_data(aa)$ITS_nonspike)
sample_data(aa)$STS_percent_spikes <- sample_data(aa)$Spike_Reads_16S/(sample_data(aa)$Spike_Reads_16S + sample_data(aa)$Nonspike_Reads_16S)

#Marker specific spike filtering, threshold based on Tkacz 2018
aa_ITS <- subset_samples(aa, ITS_percent_spikes >= 0.2 & ITS_percent_spikes <= 0.8) #20% threshold 
aa_STS <- subset_samples(aa, STS_percent_spikes >= 0.2 & STS_percent_spikes <= 0.8) #20% threshold 

#dereplicate ITS and 16S phyloseq objects
derep_aa_ITS = merge_samples(aa_ITS, "Derep_ID")
derep_aa_STS = merge_samples(aa_STS, "Derep_ID")

#Create a dereplicated template to retain sample variables
sample_data(aa_STS)$Future_rownames <- rownames(sample_data(aa_STS)) #can change dataframe for ITS samples
df <- data.frame(sample_data(aa_STS)) #Turn sample data into a dataframe object

#Remove all duplicated rows, but lose numbers
rmreps <- subset(df, !duplicated(Derep_ID)) # Remove rows with duplicate values in Unique_indiv (it will retain 1 of each Unique_indiv values)
sorted <- rmreps[order(rmreps$Derep_ID),] #Order dataframe by Derep_ID
rownames(sorted) <- sorted$Derep_ID 

#Fix spike data numbers
fixed_spikes <- df %>%
  group_by(Derep_ID) %>%
  summarize(ITS_summed_spikes = sum(ITS_Spike), ITS_summed_nonspikes = sum(ITS_nonspike), 
            STS_summed_spikes = sum(Spike_Reads_16S), STS_summed_nonspikes = sum(Nonspike_Reads_16S)) #Add total spike and nonspike reads within individual grouped sample

spike_adj_df <- full_join(sorted, fixed_spikes, by = "Derep_ID") #Add adjusted spike info to merged dataset
spike_adj_df$STS_factor_adj <- median(spike_adj_df$STS_summed_spikes)/spike_adj_df$STS_summed_spikes  #create factor to adjust OTU table to absolute abundance

rownames(spike_adj_df) <- spike_adj_df$Derep_ID #correct rownames
head(spike_adj_df) #sanity check, look good?
sample_data(derep_aa_STS) <- spike_adj_df #Replace dataframe of derep with adjusted dataframe

#Reorient OTU table which is toggled during merge, so that taxa_are_rows = TRUE
otu_table(derep_aa_STS) <- t(otu_table(derep_aa_STS))

#clean up to only plant reads
ps_kingdom_L <- subset_taxa(derep_aa_STS, Kingdom == "d__Bacteria") #Remove non-bacteria
ps_kingdom_L  <- subset_taxa(ps_kingdom_L, Order != " o__Chloroplast") #Remove chloroplasts

#Minimal filtering
aa_nolow <- prune_samples(sample_sums(ps_kingdom_L)>500, ps_kingdom_L) #Remove samples with 500 or less reads
aa_norare <- prune_taxa(taxa_sums(aa_nolow)>9, aa_nolow) #Remove ASVs with less than 10 reads overall

aa_exp_only <- subset_samples(aa_norare, Tissue != "Mock" & Tissue != "NC" & Tissue != "silisen")
aa_exp_only <-prune_taxa(taxa_sums(aa_exp_only) >0, aa_exp_only)
aa_plant_only <- subset_samples(aa_norare, Tissue != "Mock" & Tissue != "NC" & Tissue != "soil" & Tissue != "silisen")



##########################
#### Transformations #####
##########################

#
##
### Repeat rarefy: needed for both raw data or rarefied absolute abundance adjusted

#Selected depths from MAMPR_2_alpha_diversity: 16S = 1380, ITS = 751
#Can use raw data set (exp_only) or absolute abundance adjusted dataset (aa_exp_only)

#Generate rarefied table
exp_only <- prune_taxa(taxa_sums(exp_only)>0, exp_only) #Remove ASVs not present in dataset
#Get a rarefied table for future manipulation
rare_example <- otu_table(rarefy_even_depth(exp_only, sample.size=1380, trimOTUs =FALSE, #can set sampling depth here
                                            replace = FALSE, verbose = TRUE, rngseed = 8)) 

#Set up blank matrix for repeat rarefying
dim(rare_example) #Determine dimensions of rarefied object
total_matrix <- matrix(0, nrow = 19573, ncol = 960) #Create blank OTU table with above dimensions

#Repeat rarefy
seed_set <- sample.int(1e6, 100) # generate set of seeds for sampling
for (y in seed_set){
  rare_obj <- otu_table(rarefy_even_depth(exp_only, sample.size=1380, trimOTUs =FALSE,
                                          replace = FALSE, verbose = TRUE, rngseed = y)) #extract rarefied otu table for each seed
  total_matrix <- total_matrix + rare_obj #Iteratively add matrices to get single total score
}

rare_sum_matrix <- total_matrix #Copy rarefied matrix
mean_ot <- rare_sum_matrix/100 #Average sums to get representative OTU table

#Save rarefied matrix since it takes a real long time to run
saveRDS(mean_ot, "Rep_Rare_100x_Mean_OTU_table_exponly_1380depth.rds")


#If pre-run OTU table, load repeat rarefied table into phyloseq object for downstream
mean_ot <- readRDS("Rep_Rare_100x_Mean_OTU_table_exponly_1380depth.rds")


#Replace OTU table in phyloseq object with new rarefied matrix
rare_phylo <- rarefy_even_depth(exp_only, sample.size=1380, trimOTUs =FALSE,
                                replace = FALSE, verbose = TRUE, rngseed = 8) #Create phyloseq object template

otu_table(rare_phylo) <- otu_table(mean_ot, taxa_are_rows = T) #Replace OTU table with new mean rarefied OTU table

# rare_phylo is complete for raw adjusted
# to actually adjust rarefied data set to absolute counts, continue to next step

#
##
### ABSOLUTE ABUNDANCE ADJUSTMENTS

# Adjust separately for ALR and rarefied phyloseqs


#RAREFIED ABSOLUTE ADJUSTMENT
aa_choice <- aa_rare_phylo #rarefied dataset
aa_ps <- filter_taxa(aa_choice, function(x) sum(x >= 1) > 0, TRUE) #One or more read in at least 1 sample after rarefying

#Extract 16S adjustment factor vector, use sample_data(aa_ps)$ITS_factor_adj for ITS
sample_vec <- sample_data(aa_ps)$STS_factor_adj 

#Transform with OTU table absolute abundance factor
abs_mat <- sweep(as.matrix(otu_table(aa_ps)), 1 + taxa_are_rows(aa_ps), 
                 sample_vec, FUN = '*') 

aa_final_rare <- aa_ps #copy original phyloseq object
otu_table(aa_final_rare) <- otu_table(abs_mat, taxa_are_rows = TRUE) #Replace OTU table with absolute adjusted matrix





#ALR ADJUSTMENT
aa_ps  <- aa_plant_only #NOT rarefied dataset

#Extract transformation vector
sample_vec <- sample_data(aa_ps)$STS_factor_adj #16S adjustment

#Transform OTU table with absolute abundance factor
abs_mat <- sweep(as.matrix(otu_table(aa_ps)), 1 + taxa_are_rows(aa_ps), 
                 sample_vec, FUN = '*')

# Log transform for "robust" ALR - already adjusted by constant in above step
alr_mat <- log(abs_mat) #Transforms zeros to infinite
alr_mat[is.infinite(alr_mat)] <- 0 # Convert zeros back to zeros

aa_alr <-  aa_ps #copy phyloseq object
otu_table(aa_alr) <- otu_table(alr_mat, taxa_are_rows = TRUE) #Add abundance adjusted, log-transformed, zero-controlled ALR

aa_alr #ALR phyloseq object ready for downstream




#CLR
plant_only #take untransformed plant_only through core filtering





#Amalgamate taxonomy at a higher level

# Select phyloseq object
glom_choice <- rare_phylo

#Combine by branche length
glom_trial <-  tree_glom(glom_choice, 0.05) #Change number to reflect taxonomic collapse number

plot_tree(glom_trial, "treeonly", nodeplotblank, label.tips="taxa_names") #View tree
sort(table(as.data.frame(tax_table(glom_trial))$Genus)) #View collapse on taxonomic level. A lot of collapse within genera level with 0.05, a few genera also collapse (maybe NAs with IDed?)
ntaxa(glom_trial) #Number of taxa left after glom
length(unique(as.data.frame(tax_table(glom_trial))$Genus)) #Number of unique genera left after glom

#chosen glom (set numbers to chosen depth)
glommed_ps <- tree_glom(glom_choice, 0.05)


#glom by taxa
#First, set all unclassified and uncultured to NA in taxonomy dataframe
tt <- as.data.frame(tax_table(plant_only))
ts <- data.frame(lapply(tt, function(x) {
  gsub(".*uncultured*.", NA, x)
}))
rownames(ts) <- rownames(tt)

glom_ps <- plant_only #copy phyloseq object
tax_table(glom_ps) <- as.matrix(ts) #Add dataset with NA designation for unclassified ranks to phyloseq object

genus_glom <- tax_glom(glom_ps, taxrank = "Genus") #Glom by taxonomic rank


#ASV losses for classification at different taxa levels
dim(subset(ts, is.na(ts$Genus) == TRUE)) 
5009/20181 
#Genus = 24.8% loss

dim(subset(ts, is.na(ts$Family) == TRUE)) 
1273/20181
#Family = 6.31 % loss

dim(subset(ts, is.na(ts$Order) == TRUE)) 
78/20181
#Order = 0.387 % loss

dim(subset(ts, is.na(ts$Class) == TRUE)) 
29/20181 
#Class = 0.144 % loss

dim(subset(ts, is.na(ts$Phylum) == TRUE)) 
#Phylum = complete 






##########################
##### Core selection #####
##########################

#options to run through cores:
rare_phylo #Rarefied dataset




ps_x <- subset_samples(rare_phylo, Tissue != "soil") #Select phyloseq and remove soil
ps_x <- subset_samples(ps_x,  Future_rownames %in% normal_list) 


#Core option 1, CORE A: prevalence + abundance filters apply to entire plant data set
core_a <- core(ps_x , detection = 1/200, prevalence = (4/nsamples(ps_x))) #At least 0.005% RA in 4 samples across total dataset
core_a <- prune_samples(sample_sums(core_a)>0, core_a) #Remove lost ASVs


##Core option 2, CORE B: prevalence + abundance filters apply to tissue x stage subsets

stage_tissue_list <- unique(sample_data(ps_x)$Stage_Tissue) #extract vector of all tissue x stage subsets

all_core <- c() #empty vector to accumulate saved ASVs

#find species that are at at least 1% RA in 20% of each tissue x stage subset and append to list
for (st in stage_tissue_list){
  b <- subset_samples(ps_x, Stage_Tissue == st) #subset phyloseq object by tissue x stage subset
  core_b <- core_members(b, detection = 1/100, prevalence = 2/10) #filter ASVs within subset to reach 1% RA in 20% of tissues and create a vector of retained ASVs
  all_core <- c(all_core, core_b) #append vector to dataset vector
}

core_b <- prune_taxa(taxa_names(ps_x) %in% all_core, ps_x) #retain only ASVs that are meet thresholds

### Core option 3, Plant Indicators; ASVs associated with plants/plant tissues as opposed to surrounding soil

library(indicspecies) #Load package

#Select phyloseq object. Must have soil retained.
mp <- filter_taxa(rare_phylo, function(x) sum(x >= 1) > 4, TRUE) #One or more read in at least 5 samples after rarefying

# Convert ASV counts to presence/absence (P/A) by replacing all values > 0 with 1
OTUCounts <- as.data.frame(otu_table(mp))
OTUCounts[OTUCounts>0] <- 1

# Replace OTU table with P/A table
otu_table(mp) <- otu_table(as.matrix(OTUCounts),taxa_are_rows=TRUE)

# Make categories for soil vs. rhizosphere vs. phyllosphere.
sample_type <- data.frame(sample_data(mp))
sample_type$Sphere <- ifelse(sample_type$Tissue=="root","Rhizo", 
                             ifelse(sample_type$Tissue == "soil", "Soil", "Phylo"))
sample_type$IsPlant <- ifelse(sample_type$Tissue == "soil", "No", "Yes")

category_list <- c("Tissue","Sphere","stage_score") #,"IsPlant" #Create list of categories of indicator species

#Find plant indicators
output_file_path <- "STS_Indicator_Species" #create file to save indicator species

for(group in category_list){ # For each grouping variable
  print(group)
  group_list <- as.character(sample_type[order(sample_type[,group]),][,group]) # Order samples by group and get the names of the groups to which the ordered samples belong
  group_order <- as.character(rownames(sample_type[order(sample_type[,group]),])) # Get the order of samples arranged by group
  species_matrix <- t(otu_table(mp)) # Transform the ASV matrix so that ASVs are columns and samples are rows
  sample_type_ordered <- species_matrix[match(group_order,rownames(species_matrix)),] # Order the rows (samples) by group
  
  # Permutation test to find species significantly associated with sample group:
  prefsign = signassoc(sample_type_ordered, cluster=group_list, mode=0, alternative = "two.sided",control = how(nperm=9999)) 
  df <- as.data.frame(prefsign) # Put permutation test results into dataframe
  indicator_ASVs <- df[df$psidak < 0.01,] # Select results with desired p-value threshold
  assign(paste(group,"Indicators",sep="_"),indicator_ASVs) # Store results
  write.table(indicator_ASVs, file=paste(output_file_path,paste(group,"Indicators",sep="_"),sep=""), sep="\t", row.names=TRUE,quote=FALSE) # Save results
}

# Make dataframe for storing results
VennFrame <- NULL
VennFrame <- data.frame(matrix(ncol = length(category_list), nrow = length(colnames(species_matrix))),stringsAsFactors = FALSE)
VennFrame  <- setNames(data.frame(VennFrame), category_list)
rownames(VennFrame) <- colnames(species_matrix)

for(group in category_list){
  VennFrame[,group] <- ifelse(rownames(VennFrame) %in% rownames(get(paste(group,"Indicators",sep="_"))),1,0)
}

# ASVs significantly associated with tissue type: 460
colSums2(as.matrix(VennFrame[,"Tissue"]))
colSums2(as.matrix(VennFrame[,"Sphere"]))

# The "best" column indicates the tissue with the strongest association by the number of the column. 
# Add the name of the tissue for which each ASV is an indicator in the "TopTissue" column.
Tissue_Indicators$TopTissue <- NA
Tissue_Indicators[Tissue_Indicators$best==1,]$TopTissue <- colnames(Tissue_Indicators)[1]
Tissue_Indicators[Tissue_Indicators$best==2,]$TopTissue <- colnames(Tissue_Indicators)[2]
Tissue_Indicators[Tissue_Indicators$best==3,]$TopTissue <- colnames(Tissue_Indicators)[3]
Tissue_Indicators[Tissue_Indicators$best==4,]$TopTissue <- colnames(Tissue_Indicators)[4]
Tissue_Indicators[Tissue_Indicators$best==5,]$TopTissue <- colnames(Tissue_Indicators)[5]
Tissue_Indicators[Tissue_Indicators$best==6,]$TopTissue <- colnames(Tissue_Indicators)[6]
Tissue_Indicators[Tissue_Indicators$best==7,]$TopTissue <- colnames(Tissue_Indicators)[7]

#Sphere Indicators
Sphere_Indicators$TopSphere <- NA
Sphere_Indicators[Sphere_Indicators$best==1,]$TopSphere <- colnames(Sphere_Indicators)[1]
Sphere_Indicators[Sphere_Indicators$best==2,]$TopSphere <- colnames(Sphere_Indicators)[2]
Sphere_Indicators[Sphere_Indicators$best==3,]$TopSphere <- colnames(Sphere_Indicators)[3]


# Make dataframe for storing results
VennFrame2 <- NULL
VennFrame2 <- data.frame(matrix(ncol = length(unique(Tissue_Indicators$TopTissue)), nrow = length(colnames(species_matrix))),stringsAsFactors = FALSE)
VennFrame2  <- setNames(data.frame(VennFrame2), unique(Tissue_Indicators$TopTissue))
rownames(VennFrame2) <- colnames(species_matrix)

for(tissue in unique(Tissue_Indicators$TopTissue)){
  # print(tissue)
  VennFrame2[,tissue] <- ifelse(rownames(VennFrame2) %in% rownames(Tissue_Indicators[Tissue_Indicators$TopTissue==tissue,]),1,0)
  # table(VennFrame[,group])
}

#Plant ASV list (either associated with a specific plant tissue or plant sphere)
plant_associated_ASVs <- unique(c(subset(Sphere_Indicators, TopSphere != "Soil")$Taxa_name, subset(Tissue_Indicators, TopTissue != "soil")$Taxa_name))

#Filter phyloseq object for indicator speices with plant associated ASVs
mp_po <- subset_samples(mp, Tissue != "soil")
mp_po <- prune_taxa(taxa_names(mp_po) %in% (plant_associated_ASVs), mp_po)




#####################################
##### Beta diversity statistics #####
#####################################


##### Bray-Curtis and Jaccard #######


#Rarefied, Absolute abundance rarefied

beta_choice <- core_b
beta_df <- data.frame(sample_data(beta_choice))

#Distances on rarefied dataset

##Phyloseq 
bray_dist <- phyloseq::distance(beta_choice, method="bray")
jaccard_dist <- phyloseq::distance(beta_choice, method="jaccard")

bray_perm <- adonis2(bray_dist ~ STS_MiSeq_Run/plate_id + Tissue*stage_score*genotype,
        beta_df,
        permutations=999)
bray_perm


jaccard_perm <- adonis2(jaccard_dist ~ STS_MiSeq_Run/plate_id + Tissue*stage_score*genotype,
        beta_df,
        permutations=999)
jaccard_perm


###### Weighted UniFrac: 16S only #######

#phyloseq's weighted unifrac is flawed (https://github.com/joey711/phyloseq/issues/956). Use QIIME2's Unifrac
#Convert all phyloseq objects to qiime2 objects


#define function to take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
# function modified from cedwardson4, phyloseq-extras https://github.com/cedwardson4/phyloseq-extras/blob/master/phyloseq2QIIME2.R
phyloseq2qiime2<-function(physeq){
  library(phyloseq)
  library(biomformat)
  library(ape)
  library(Biostrings)
  library(dada2)
  if(packageVersion("biomformat") < "1.7") {
    stop("This will only work with biomformat version > 1.7")
  }
  ps_name <-deparse(substitute(physeq))
  taxa_are_rows_logical<-taxa_are_rows(physeq)
  #write OTU table to biom file
  if(is.null(access(physeq,"otu_table"))==FALSE){
    if(taxa_are_rows_logical==TRUE) {
      otu<-as(otu_table(physeq),"matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    } else if (taxa_are_rows_logical==FALSE) {
      otu<-t(as(otu_table(physeq),"matrix"))
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    }
  }
  #write sample data (metadata) to tsv
  if(is.null(access(physeq,"sam_data"))==FALSE){
    caro_edit <- sample_data(physeq)
    caro_edit$sampleid <- rownames(caro_edit) #Preserve sample names of metadata
    caro_edit <- as_tibble(caro_edit)%>%
      relocate(sampleid) #Move sample id column to the first column of dataframe
    write.table(caro_edit,file=paste0(ps_name,"_sample-metadata.txt"), 
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
  }
  #write taxonomy table to qiime2 formatted taxonomy
  if(is.null(access(physeq,"tax_table"))==FALSE){
    tax<-as(tax_table(physeq),"matrix")
    tax_cols <- colnames(tax)
    tax<-as.data.frame(tax)
    tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
    for(co in tax_cols) tax[co]<-NULL
    write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
    print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
  }
  #write phylogenetic tree to newick formwat
  if(is.null(access(physeq,"phy_tree"))==FALSE){
    if(is.rooted(phy_tree(physeq))==TRUE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
      print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
    } else if (is.rooted(phy_tree(physeq))==FALSE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
      print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
    }      
  }
  #write representative sequences to fasta format
  if(is.null(access(physeq,"refseq"))==FALSE){
    writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))      
  } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))    
  }
}


phyloseq2qiime2(beta_choice) #convert phyloseq object into QIIME2 objects

#####   switch to qiime2 environment in Unix   #####
##### IN QIIME2 COMMAND LINE WILL NOT RUN IN R #####

#Make qiime artifact for OTU Table
qiime tools import \
--input-path beta_choice_features-table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path TAB_beta_choice.qza #Set name you want here

#Make qiime artifact for phylogenetic tree
qiime tools import \
--input-path rare_core_a_beta_choice.newick \ #Ready for QIIME2 file
--output-path ROOTED_TREE_beta_choice.qza \ #Save as QIIME2 artifact
--type 'Phylogeny[Rooted]'

#Run calculate weighted Unifrac distances
qiime diversity beta-phylogenetic \
--i-table TAB_beta_choice.qza \
--i-phylogeny ROOTED_TREE_beta_choice.qza \
--p-metric weighted_unifrac \
--o-distance-matrix weighted_unifrac_distance_matrix_beta_choice.qza



#####          switch to R        #####
##### IN R WILL NOT RUN IN QIIME2 #####

dist <- read_qza('distance_matrix.qza')$data #Import distance matrix from QIIME2
beta_df <- data.frame(sample_data(beta_choice)) #extract metadata from original phyloseq object

perm <- adonis2(idist ~ STS_MiSeq_Run/plate_id + Tissue*stage_score*genotype,
                 beta_df,
                 permutations=999)
perm






##Pairwise post hoc tests for Ripe Siliques
#Rarefied, Absolute abundance rarefied


ms_ph<- subset_samples(core_b, stage_score == "mature_siliques") #Only select Ripe Siliques stage
ms_ph <- subset_samples(ms_ph, genotype == "col0" | genotype == "lore") #select genotype, swap out lore for "lyk4", "efr", "fls2"

beta_df <- data.frame(sample_data(ms_ph))

#Statistics
##Phyloseq 
bray_dist <- phyloseq::distance(ms_ph, method="bray")

bray_perm <- adonis2(bray_dist ~ Tissue*genotype,
                     beta_df,
                     permutations=999)
bray_perm








###### Aitchison #######

#function to convert phyloseq object into object for vegan distances
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

#ALR
 #ALR prepped phyloseq object
beta_choice <- aa_alr_core
beta_df <- data.frame(sample_data(beta_choice)) #save metadata
vegan_beta <- psotu2veg(beta_choice) #Convert to vegan compatible object

alr_dist <- vegdist(vegan_beta, method = "euclidean") #Calculate Aitchison distance for ALR

#PERMANOVA
adonis2(alr_dist ~ STS_MiSeq_Run/plate_id + Tissue*stage_score*genotype,
        beta_df,
        permutations=999)

#CLR
beta_choice <- plant_only_core
beta_choice <- microbiome::transform(beta_choice, "compositional") #transform to relative abundance
beta_df <- data.frame(sample_data(beta_choice)) #save metadata

vegan_beta <- psotu2veg(beta_choice) #transform to vegan compatible object
clr_dist <- vegdist(vegan_beta, method = 'robust.aitchison') #Run robust CLR/transformation and distance calculation

#PERMANOVA
adonis2(clr_dist ~ STS_MiSeq_Run/plate_id + Tissue*stage_score*genotype,
        beta_df,
        permutations=999)









###################################
##### Beta diversity graphics #####
###################################

y <- core_b #select phyloseq object
ord <- phyloseq::ordinate(y, method = "PCoA", dist = bray_dist) #ordinate graph

#Create dataframe for in depth mapping

meta <-data.frame(sample_data(y)) %>%
  rownames_to_column("SampleID") #extract metadata

#create dataframe with vectors combined with metadata
dt <- data.frame(ord$vectors) %>%
  rownames_to_column("SampleID")  %>% #vector frames
  left_join(meta) %>% # join ordination vectors with sample data
  mutate(
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



#Centroids of groups (demo = Tissue x Genotype subsets)
centroid_all <- dt %>%
  dplyr::group_by(Tissue_Genotype) %>%
  summarize(across(where(is.numeric), mean)) %>%
  separate(Tissue_Genotype, c("Tissue", "Genotype")) %>%
  mutate(
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


dt$genotype <- factor(dt$genotype, levels = c("col0", "efr", "fls2", "lore", "lyk4")) #Set order of genotypes
dt$Nice_Tissue <- factor(dt$Nice_Tissue, levels = c("Root", "Rosette", "Stem", "Cauline Leaves", "Flowers", "Siliques")) #Set order of Tissues
axis_x = "Axis.1" #Choose Axis 1 to view ordination 
axis_y = "Axis.2"  #Choose Axis 2 to view ordination 


genotype_labels <- c("WT", "*efr*", "*fls2*", "*lore*", "*lyk4*") #set labels for genotype (italics)
#Generate ggplot object with sample points
sample_plot <- ggplot() +
  theme_bw() +
  xlab(paste("Axis 1: ", round(100*ord$values$Relative_eig[1]),"%")) + #Select plot axis x label
  ylab(paste("Axis 2: ", round(100*ord$values$Relative_eig[2]),"%")) + #Select plot axis y label
  
  #plot data
  geom_point(
    data = dt,
    cex = 1.5, 
    aes_string(x=axis_x, y=axis_y, color= "genotype" , shape = "Nice_Tissue" )) +
  scale_color_discrete(labels = genotype_labels,
                       type = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) + #set colors
  scale_shape_manual(values = c(15, 16, 17, 3, 8, 11)) + #set shapes
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_markdown(size = 12))+
  guides(fill = guide_legend(override.aes = list(size=22))) +
  labs(shape = "Tissue", color = "Genotype")
sample_plot


#Add group centroids to plot
s <- sample_plot + geom_point(data = centroid_all,
                              aes_string(x = axis_x, y = axis_y, color = "Genotype", shape = "Nice_Tissue"),
                              size = 4) 

s


#Add ellipses
s + stat_ellipse(data=subset(dt, Nice_Tissue == "Root" | Nice_Tissue == "Rosette"),
                 aes_string(x=axis_x, y=axis_y, color="genotype", shape = "Nice_Tissue"),
                 type = "t", level = 0.85)








##################################################
############### BETA DISPERSION ##################
##################################################

perm_dispr <- data.frame() #Create dataframe to store results

rare_phylo <- subset_samples(rare_phylo, Tissue != "soil") #Remove soil from phyloseq object
rare_abs <- subset_samples(rare_abs, Tissue != "soil") #Remove soil from phyloseq object

phylo_list <- list("Rarefied" = rare_phylo, "AbsAbunRarefied" = rare_abs) #named list of phyloseq objects
index_list <- c("bray", "jaccard") #diversity index list

#rarefied phyloseqs list
for (ps in names(phylo_list)){
  #Apply core filters
  
  min_filt <- filter_taxa(phylo_list[[ps]], function(x) sum(x >= 1) > 0, TRUE) #Filter extremely rare taxa
  
  #Core A
  core_a <- core(phylo_list[[ps]] , detection = 1/200, prevalence = (4/nsamples(phylo_list[[ps]]))) #At least 0.005% RA in 4 samples
  core_a <- prune_samples(sample_sums(core_a)>0, core_a)
  
  #Core B
  stage_tissue_list <- unique(sample_data(phylo_list[[ps]])$Stage_Tissue) #extract all tissue x stage subsets
  all_core <- c() #empty list for loop
  
  #find species that are at at least 1% RA in 20% of each tissue + time subset and append to list
  for (st in stage_tissue_list){
    b <- subset_samples(phylo_list[[ps]], Stage_Tissue == st)
    core_b_list <- core_members(b, detection = 1/100, prevalence = 2/10)
    all_core <- c(all_core, core_b_list)
  }
  
  core_b <- prune_taxa(taxa_names(phylo_list[[ps]]) %in% all_core, phylo_list[[ps]])
  
  #Indicator core
  indic_core <- subset_samples(phylo_list[[ps]], Tissue != "soil")
  indic_core <- prune_taxa(taxa_sums(indic_core) >= 1, indic_core)
  indic_core <- prune_taxa(taxa_names(indic_core) %in% (plant_associated_ASVs), indic_core)
  
  #Create list of core phyloseq objects
  core_list <- list("Minimal_filtering" = min_filt, "Core_A" = core_a,  "Core_B" = core_b, "Indicator_Core" = indic_core) 
  
  for(core in names(core_list)){
    b_df <- as.data.frame(as.matrix(sample_data(core_list[[core]])))
    
    #calculate beta dispersion for overall dataset
    for(index in index_list){
      bdist <- phyloseq::distance(core_list[[core]], method=index)
      beta_dis <- betadisper(d = bdist, group = b_df$genotype, type = "median")
      beta_perm <- permutest(beta_dis, perumtations = 999)
      dispr_df <- data.frame(beta_perm$tab)
      dispr_df$test <- "Overall Genotype"
      dispr_df$core <- core
      dispr_df$index <- index
      dispr_df$phyloseq <- ps
      perm_dispr <- rbind(perm_dispr, dispr_df) #Add to data frame
      
      
      # Do dispersion on subsets of data
      #Within tissue
      for(x in tissue_list){
        ss <- subset_samples(core_list[[core]], Tissue == x)
        ss <- prune_taxa(taxa_sums(ss) > 0, ss)
        ss_df <- as.data.frame(as.matrix(sample_data(ss)))
        idist <- phyloseq::distance(ss, method=index)
        beta_ss <- betadisper(d = idist, group = ss_df$genotype, type = "median")
        perm <- permutest(beta_ss, perumtations = 999)
        perm_df <- data.frame(perm$tab)
        perm_df$test <- paste0("Within_",x)
        perm_df$core <- core
        perm_df$index <- index
        perm_df$phyloseq <- ps
        perm_dispr <- rbind(perm_dispr, perm_df)
      }
      #Within Stage
      for(x in stage_list){
        ss <- subset_samples(core_list[[core]], stage_score == x)
        ss <- prune_taxa(taxa_sums(ss) > 0, ss)
        ss_df <- as.data.frame(as.matrix(sample_data(ss)))
        idist <- phyloseq::distance(ss, method=index)
        beta_ss <- betadisper(d = idist, group = ss_df$genotype, type = "median")
        perm <- permutest(beta_ss, perumtations = 999)
        perm_df <- data.frame(perm$tab)
        perm_df$test <- paste0("Within_",x)
        perm_df$core <- core
        perm_df$index <- index
        perm_df$phyloseq <- ps
        perm_dispr <- rbind(perm_dispr, perm_df)
      }
      #within Stage x Tissue subsets
      for(x in stage_tissue_list){
        ss <- subset_samples(core_list[[core]], Stage_Tissue == x)
        ss <- prune_taxa(taxa_sums(ss) > 0, ss)
        ss_df <- as.data.frame(as.matrix(sample_data(ss)))
        idist <- phyloseq::distance(ss, method=index)
        beta_ss <- betadisper(d = idist, group = ss_df$genotype, type = "median")
        perm <- permutest(beta_dis, perumtations = 999)
        perm_df <- data.frame(perm$tab)
        perm_df$test <- paste0("Within_",x)
        perm_df$core <- core
        perm_df$index <- index
        perm_df$phyloseq <- ps
        perm_dispr <- rbind(perm_dispr, perm_df)
      }
    }
  }
}

write.table(perm_dispr, file = "BETADISPER_all_varieties_16S.csv", sep = ",")

