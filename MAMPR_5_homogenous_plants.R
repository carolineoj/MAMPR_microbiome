####################################################
############### HOMOGENOUS PLANTS ##################
####################################################

#phyloseq object preprocessing
rare_po <- subset_samples(rare_phylo, Tissue != "soil") #Remove soil
rare_abs <- subset_samples(rare_abs, Tissue != "soil") 

ps <- filter_taxa(rare_po, function(x) sum(x >= 1) > 0, TRUE) #Filter extremely rare taxa

#
##
###Statisics

#define list of comparision subsets
list_of_tis_lists <- list(c("root", "rose", "stem", "caul", "flow", "silinew"), 
                          c("rose", "stem", "caul", "flow", "silinew"),
                          c("rose", "stem", "caul", "flow"),
                          c("stem", "caul", "flow"),
                          c("root", "rose"),
                          c("rose", "stem"),
                          c("caul", "stem"),
                          c("flow", "stem"),
                          c("silinew", "stem"))

#define beta diversity metrics to test
dist_list <- c("bray", "jaccard")

#Define stages to test
stage_list <- c("veg", "flow_only", "young_siliques", "mature_siliques")


#Create dataframe for output
homogeneous_df <- data.frame()
phylo_list <- list("Rarefied" = rare_po, "AbsAbunRarefied" = rare_abs) #named list of phyloseq objects
index_list <- c("bray", "jaccard") #diversity index list

#rarefied phyloseqs list
for (ps in names(phylo_list)){
  
  #Define filtering/core level
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
  
  
  for(core in names(core_list)){ #for each filtering level
    for(tis_list in list_of_tis_lists){ #for set of tissues we could examine
      rots <- subset_samples(core_list[[core]], Tissue %in% tis_list) #E.g. select only root and rosette samples in one iteration
      ID_list <- sample_data(rots)$Indiv_ID #list of all individual plant IDs
      n_occur <- data.frame(table(ID_list)) #Table of frequencies of all IDs, i.e. number of tissues each individiual has in the dataset
      complete <- n_occur[n_occur$Freq == length(tis_list),] #List of all plants that have the complete set of tissues represented in the dataset
      whole_list <- complete$ID_list #Create list
      
      sub_ps <- subset_samples(rots, Indiv_ID %in% whole_list) #Subset phyloseq to have complete samples
      sub_ps <- prune_taxa(taxa_sums(sub_ps)>0, sub_ps) #Remove lost taxa
      
      for(d in dist_list){
        sub_dist <- as.matrix(phyloseq::distance(sub_ps, method = d)) #Create distance object
        wsd <- sub_dist[str_detect(row.names(sub_dist), paste(whole_list, collapse = '_|')), 
                        str_detect(row.names(sub_dist), paste(whole_list, collapse = '_|'))] #subset matrix by the complete sample list
        ad <- as.dist(wsd) #Convert matrix into dist object for betadisper
        data <- data.frame(sample_data(subset_samples(ps, Indiv_ID %in% whole_list & Tissue %in% tis_list))) # #subset sample data frame
        bd <- betadisper(ad, group = data$Indiv_ID) #Make betadisper object
        
        #overall distances (not separated by stage)
        dist_to_med <- bd$group.distances #Distances to Individual medians
        overall_df <- data.frame(dist_to_med) #Dataframe of all individual distances (not subsetted by timepoint)
        overall_df$genotype <- unlist(lapply(row.names(overall_df), function (x) strsplit(as.character(x), "-")[[1]][1])) 
        names(overall_df)[names(overall_df) == "dist_to_med"] <- "Dist_to_median" #alter dataframe to target genotype
        ko <- kruskal.test(Dist_to_median ~ genotype, overall_df) #Effect of genotype
        
        #Effect of stage
        overall_df$Indiv_ID <- row.names(overall_df)
        ifelse("rose" %in% tis_list == TRUE, rose_df <- subset(data, Tissue == "rose"), rose_df <- subset(data, Tissue == "stem"))
        ifelse("rose" %in% tis_list == TRUE, comp_tis <- "rose", comp_tis <- "stem")
        
        #Select single row per sample with rosette weight
        df_with_dist <- full_join(overall_df, rose_df, by = "Indiv_ID") #Add metadata
        df_with_dist$Dist_to_median <- as.numeric(df_with_dist$Dist_to_median) #Convert to numeric
        df_with_dist$Final_weight_mg <- as.numeric(df_with_dist$Final_weight_mg) #Convert to numeric
        k_stage <- kruskal.test(Dist_to_median ~ stage_score, df_with_dist)
        
        #Rosette size vs distance
        ct_overall <- cor.test(df_with_dist$Final_weight_mg, df_with_dist$Dist_to_median, method=c("spearman"))
        
        
        
        #Rosette distance within stage
        for(stage in stage_list){
          df_stage <- subset(df_with_dist, stage_score == stage)
          if (nrow(df_stage) == 0) next
          ct_stage <- cor.test(df_stage$Final_weight_mg, df_stage$Dist_to_median, method=c("spearman"))
          
          
          #within stage differences for genotype
          stage_dist <- dist_to_med[names(dist_to_med) %in% stage] #subset dataframe by stage
          df <- data.frame(stage_dist) #make stage specific dataframe
          df$genotype <- unlist(lapply(row.names(df), function (x) strsplit(as.character(x), "-")[[1]][1])) 
          names(df)[names(df) == "stage_dist"] <- "Dist_to_median" #alter dataframe to target genotype
          ks <- kruskal.test(Dist_to_median ~ genotype.x, df_stage)
          row_vec <- cbind(ps, core, length(tis_list), paste0(tis_list, collapse = ""), k_stage$p.value, 
                           comp_tis, ct_overall$p.value, ct_stage$p.value, 
                           ko$p.value, ks$p.value, d, stage)
          homogeneous_df <- rbind(homogeneous_df, row_vec) #Add data to dataframe
        }
      }
    }
  }
}

#add column names to distance data frame
hom_df_colnames <- c("Phyloseq_Object", "Core_used", "Number_of_Tissues", "Involved_Tissues", "Stage_homogenity_pval", "Weighed_Tissue", 
                     "All_weight_homogenity_correlation_pval", "Stage_specific_weight_homogenity_cor_pval", 
                     "Genotype_homo_overall_pval", "Genotype_homo_by_stage_pval","Distance_Metric", "Stage")

colnames(homogeneous_df) <- hom_df_colnames

write.csv(homogeneous_df, 
          "Homogeneous_plants_16S_rarefied_multimetric.csv")



tis_list <- c("rose", "stem", "caul", "flow")
d = "bray"




#
##
### Graphics
aerial_nr <- c("rose", "stem", "flow", "caul") #select subset to graph -includes rosettes, stems, flowers and cauline leaves
tis_list <- aerial_nr
target_ps <- core_b

rots <- subset_samples(target_ps, Tissue %in% tis_list) #Select only root and rosette samples
ID_list <- sample_data(rots)$Indiv_ID #list of all individual IDs
n_occur <- data.frame(table(ID_list)) #Table of frequencies of all IDs
complete <- n_occur[n_occur$Freq == length(tis_list),] #List of all samples that have both root and rosettes for individual IDs
whole_list <- complete$ID_list #List of whole plants

sub_ps <- subset_samples(rots, Indiv_ID %in% whole_list) #Subset to only relevant samples
sub_ps <- prune_taxa(taxa_sums(sub_ps)>0, sub_ps)
sub_dist <- as.matrix(phyloseq::distance(sub_ps, method = "bray")) #Create distance object
wsd <- sub_dist[str_detect(row.names(sub_dist), paste(whole_list, collapse = '_|')), 
                str_detect(row.names(sub_dist), paste(whole_list, collapse = '_|'))] #subset matrix by the complete sample list
ad <- as.dist(wsd) #Convert matrix into dist object for betadisper
#data <- data.frame(sample_data(subset_samples(sub_ps, Indiv_ID %in% whole_list & Tissue %in% tis_list))) # #subset sample data frame

data <- data.frame(sample_data(sub_ps))
bd <- betadisper(ad, group = data$Indiv_ID) #Make betadisper object

#overall distances (not separated by stage)
dist_to_med <- bd$group.distances #Distances to Individual group medians
overall_df <- data.frame(dist_to_med) #Dataframe of all individual distances (not subsetted by timepoint)
overall_df$genotype <- unlist(lapply(row.names(overall_df), function (x) strsplit(as.character(x), "-")[[1]][1])) 
names(overall_df)[names(overall_df) == "dist_to_med"] <- "Dist_to_median" #alter dataframe to target genotype

#Effect of stage
overall_df$Indiv_ID <- row.names(overall_df)
ifelse("rose" %in% tis_list == TRUE, rose_df <- subset(data, Tissue == "rose"), rose_df <- subset(data, Tissue == "stem"))
ifelse("rose" %in% tis_list == TRUE, comp_tis <- "rose", comp_tis <- "stem")

#Select single row per sample with rosette weight
df_with_dist <- full_join(overall_df, rose_df, by = "Indiv_ID") #Add metadata
df_with_dist$Dist_to_median <- as.numeric(df_with_dist$Dist_to_median) #Convert to numeric
df_with_dist$Final_weight_mg <- as.numeric(df_with_dist$Final_weight_mg) #Convert to numeric
df_with_dist <- df_with_dist %>%
  mutate(
    Nice_Stage = case_when(stage_score == "veg" ~ "Vegetative",
                           stage_score == "flow_only" ~ "Flowering",
                           stage_score == "mature_siliques" ~ "Ripe Siliques",
                           stage_score == "young_siliques" ~  "Unripe Siliques"))
df_with_dist$Nice_Stage <- factor(df_with_dist$Nice_Stage, levels = c("Vegetative","Flowering", "Unripe Siliques", "Ripe Siliques"))

m <-  ggplot(data = df_with_dist, aes(x = Nice_Stage, y = Dist_to_median, fill = Nice_Stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(cex = 1)+
  scale_fill_viridis_d(option = "C") +
  
  #significance bars for stage
  geom_signif(y_position = c(0.6, 0.65, 0.58), xmin = c(1, 1, 2), xmax = c(2, 3, 3),
              annotation = c("*", "NS", "NS"), tip_length = 0.05) + #Change annotation for ITS data
  
  #formatting for the graph
  theme_bw() +
  xlab("Stage") +
  ylab("Mean distance to median \n (within individual)") +
  theme(legend.position="none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_markdown()) +
  guides(fill = guide_legend(override.aes = list(size=22))) +
  ylim(0.3,0.675)

m #view plot

df_with_dist$genotype.x <- factor(df_with_dist$genotype.x, levels = c("col0", "efr", "fls2", "lore", "lyk4"))
m_geno <-  ggplot(data = df_with_dist, aes(x = genotype.x, y = Dist_to_median, fill = genotype.x))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(cex = 1)+
  scale_fill_viridis_d()+
  #sig bars for genotype 
  ggsignif::geom_signif(y_position = 0.58, xmin = 0.7, xmax = 5.3,
                        annotation = "NS", tip_length = 0) +
  theme_bw() +
  xlab("Genotype") +
  ylab("Mean distance to median \n (within individual)") +
  scale_x_discrete(
    labels = genotype_labels)+
  theme(legend.position="none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_markdown()) +
  guides(fill = guide_legend(override.aes = list(size=22))) +
  ylim(0.3,0.675)


m_geno


fig <- ggarrange(m_geno, m_ITS_geno,
                 m, m_ITS,
                 labels = c("A", "B", "C", "D"))

annotate_figure(fig, left = text_grob("Mean distance to individual median", rot = 90, vjust = 1, size = 14))


