##################################################
########### DIFFERENTIAL ABUNDANCE ###############
##################################################

#Assumes packages and objects from previous scripts MAMPR loaded


#ANCOM-BC2
library(ANCOMBC)

#Select phyloseq object
#Need to use raw counts for ANCOM-BC2 for bias estimation
da_ps <- plant_only 

#Apply core filters (pick 1/3)
#Core A
da_core_a <- core(da_ps, detection = 1/200, prevalence = (4/nsamples(da_ps))) #At least 0.005% RA in 4 samples across total dataset
da_core_a <- prune_samples(sample_sums(core_a)>0, core_a) #Remove lost ASVs

#Core B
stage_tissue_list <- unique(sample_data(da_ps)$Stage_Tissue) #extract vector of all tissue x stage subsets
all_core <- c() #empty vector to accumulate saved ASVs
#find species that are at at least 1% RA in 20% of each tissue x stage subset and append to list
for (st in stage_tissue_list){
  b <- subset_samples(da_ps, Stage_Tissue == st) #subset phyloseq object by tissue x stage subset
  core_b <- core_members(b, detection = 1/100, prevalence = 2/10) #filter ASVs within subset to reach 1% RA in 20% of tissues and create a vector of retained ASVs
  all_core <- c(all_core, core_b) #append vector to dataset vector
}
da_core_b <- prune_taxa(taxa_names(da_ps) %in% all_core, da_ps) #retain only ASVs that are meet thresholds


#Indicator species core
da_indic_core <-  prune_taxa(taxa_names(da_ps) %in% plant_associated_ASVs, da_ps)


#Relevant covariables: Categorical - Stage, Tissue, Genotype
# Group variable of interest = genotype
# Reference group = col0

#Make sure reference group is coded correctly
sample_data(da_ps)$genotype <- factor(sample(da_ps)$genotype, 
                                      levels = c("col0", "efr", "lyk4", "lore", "fls2")) #First group is reference group


#Also run with struc_zero = TRUE, then manually examine if any trends look real or are an artifact of very low representation overall
output = ancombc2(data = da_ps, assay_name = NULL,
                  tax_level = NULL, fix_formula = "Tissue + stage_score + genotype",
                  p_adj_method = "holm", prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                  group = "genotype", struc_zero = FALSE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))



# Examine list of Structural zeros
tab_zero = output$zero_ind
dim(tab_zero)
head(tab_zero)

#Loops for subsets
tissue_list <- unique(sample_data(da_ps)$Tissue) #Just within Tissue
tissue_stage_list <- unique(sample_data(da_ps)$Stage_Tissue) #within stage and tissue
stage_list <- unique(sample_data(da_ps)$stage_score) #Just within Tissue

dif_abun_by_tis <- data.frame()

#run ANCOM function
for(x in tissue_stage_list){
  ips <- subset_samples(da_ps, Stage_Tissue == x) #subset to only consider a particular Stage x Tissue subset
  ips <- core(ips, detection = 1/100, prevalence = 2/10) #core microbiome for subset
  dif = ancombc2(data = ips, assay_name = NULL, #ANCOM-BC2 function
                 tax_level = NULL, fix_formula = "stage_score + genotype",
                 p_adj_method = "holm", prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                 group = "genotype", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                             nrow = 2, 
                                                             byrow = TRUE),
                                                      matrix(c(-1, 0, 1, -1),
                                                             nrow = 2, 
                                                             byrow = TRUE)),
                                      node = list(2, 2),
                                      solver = "ECOS",
                                      B = 100))
  isdif <- dif$res %>%
    dplyr::select(taxon, contains("genotype")) %>%
    filter(diff_genotypeefr == 1 | diff_genotypelyk4 == 1 | diff_genotypelore == 1 | diff_genotypefls2 == 1) #Create dataframe of positive hits
  ifelse(nrow(isdif) == 0, print(paste0(x,": Nothing to see here!")), print(paste0(x,": Investigate further!")))
  ifelse(nrow(isdif) == 0, next, isdif$Tissue_tested <- x)
  dif_abun_by_tis <- rbind(dif_abun_by_tis, isdif) #Append to dataframe of positive hits
}


for(x in stage_list){
  ips <- subset_samples(da_ps, stage_score == x)
  ips <- prune_taxa(taxa_names(ips) %in% all_core, ips)
  dif = ancombc2(data = ips, assay_name = NULL,
                 tax_level = NULL, fix_formula = "Tissue + genotype",
                 p_adj_method = "holm", prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                 group = "genotype", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                             nrow = 2, 
                                                             byrow = TRUE),
                                                      matrix(c(-1, 0, 1, -1),
                                                             nrow = 2, 
                                                             byrow = TRUE)),
                                      node = list(2, 2),
                                      solver = "ECOS",
                                      B = 100))
  isdif <- dif$res %>%
    dplyr::select(taxon, contains("genotype")) %>%
    filter(diff_genotypeefr == 1 | diff_genotypelyk4 == 1 | diff_genotypelore == 1 | diff_genotypefls2 == 1)
  ifelse(nrow(isdif) == 0, print(paste0(x,": Nothing to see here!")), print(paste0(x,": Investigate further!")))
  ifelse(nrow(isdif) == 0, next, isdif$Tissue_tested <- x)
  dif_abun_by_tis <- rbind(dif_abun_by_tis, isdif)
}



#### Targeted Differential Abundance ######
#Target fungal (ITS) communities in Ripe Siliques for lore and col0 plants

#subset to specific timepoint and genotypes
targeted_dif_abun <- subset_samples(core_b_ITS, stage_score == "mature_siliques") #subset to only Ripe Siliques developmental stage
targeted_dif_abun <- subset_samples(targeted_dif_abun, genotype == "col0" | genotype == "lore") #Only consider lore and col0 genotypes
targeted_dif_abun <- prune_taxa(taxa_sums(targeted_dif_abun) > 0, targeted_dif_abun) #Remove lost taxa


#dif abun for tissue and tp specific core
stage_tissue_list <- unique(sample_data(targeted_dif_abun)$Stage_Tissue) #extract all tissue x stage subsets
all_core <- c() #empty list for loop

#find species that are at at least 1% RA in 20% of each tissue + time subset and append to list
for (st in stage_tissue_list){
  b <- subset_samples(targeted_dif_abun, Stage_Tissue == st)
  core_b_ITS <- core_members(b, detection = 1/100, prevalence = 2/10)
  all_core <- c(all_core, core_b_ITS)
}

core_b_ITS <- prune_taxa(taxa_names(targeted_dif_abun) %in% all_core, targeted_dif_abun)

#Try different thresholds of abundance/prevalence
strict_core <- core(targeted_dif_abun, detection = 1/100, prevalence = 5/10)

#Run ANCOM BC function (try core_b_ITS, strict_core, and also try struc_zero= TRUE)
output = ancombc2(data = strict_core, assay_name = NULL,
                  tax_level = NULL, fix_formula = "Tissue + genotype",
                  p_adj_method = "holm", prv_cut = 0.1, lib_cut = 0, s0_perc = 0.05,
                  group = "genotype", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))


#Structural zeros
tab_zero = output$zero_ind
write.table(tab_zero, "Structural_Zeros_fungal_lore_col0_sig_tp.csv", sep = ",")





#DESeq2
library(DESeq2)

core_b_ITS
diagdds <- phyloseq_to_deseq2(core_b_ITS, ~genotype) #convert to deseq2 object 

# calculate geometric means prior to estimate size factors, robust to zeros
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

#Estimate size factors
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

#Examine results
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rootvegprev)[rownames(sigtab), ], "matrix"))
head(sigtab)

loreres = results(diagdds, contrast = c("genotype", "lore", "col0"))
loreres = loreres[order(loreres$padj, na.last=NA), ]
alpha = 0.01
loresigtab = loreres[(loreres$padj < alpha), ]

sigtab
