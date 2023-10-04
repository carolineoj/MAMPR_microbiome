##################################################
############### ALPHA DIVERSITY ##################
##################################################

#Assumes loaded packages and phyloseq objects from MAMPR_1 script
#phyloseq objects:
exp_only
plant_only



###########################
## Choose sampling depth ##
###########################


#Create and view arefaction curves to select subsampling depth
#Pick phyloseq object
depth_ps <- exp_only

otu_tab <- t(abundances(depth_ps)) #Extract OTU table
sample_data(depth_ps)$Sample <- rownames(sample_data(depth_ps)) #Create a column with rownames in data frame
mdata <- as_tibble(sample_data(depth_ps)) #sample data into a tibble

#create rarefaction curve with vegan
rarecurve_data <- vegan::rarecurve(otu_tab, 
                                   step = 50, label = FALSE,
                                   cex = 0.6,
                                   xlim= c(0,2000),
                                   ylim= c(0,500))

#Add metadata to rarefaction curve data
q = map_dfr(rarecurve_data, bind_rows) %>% bind_cols(Sample = rownames(otu_tab),.) %>% 
  pivot_longer(-Sample) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name,"N", ""))) %>%
  select(-name) 
q$info <- q$Sample
q <- full_join(q,
               mdata %>% dplyr::select(Sample, stage_score, Tissue, genotype),
               by = "Sample") #full data set with rarefaction and relevant metadata

#subset dataset to only view a subset of rarefaction data
r <- subset(q, Tissue != "soil")

#create visualiztion (ggplot object) of rarefaction plot
rare_plot <- r %>%
  ggplot(aes(x=n_seqs, y=value, group = Sample, color = Tissue)) + 
  xlab("Number of reads") +
  ylab("Number of ASVs") +
  geom_line() +
  xlim(0,2000) +
  ylim(0,500)  

#graphics settings for nice plot formats
graphic_presets <- function(myplot){
  myplot +
    scale_color_frontiers() +
    scale_fill_frontiers() +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))}

#View rarefation curve 
graphic_presets(rare_plot)


# Manually examine samples that will be discarded
sample_data(depth_ps)$Sample_Depth_alpha <- sample_sums(depth_ps) #Create column with sampling depth of each sample
dim(sample_data(depth_ps)) #Total samples present in dataset (passed minimal filtering in MAMPR_1) ITS: 920 samples, 16S: 975 smaples

lowcov <- subset_samples(depth_ps, Sample_Depth_alpha < 1500) # Phyloseq object of samples lost at 1500 cut off
vlowcov <- subset_samples(depth_ps, Sample_Depth_alpha < 750) # Phyloseq object of samples lost at 750 cut off

#View which combinations of Tissue x Genotype are lost with each cut off
table(sample_data(lowcov)$Tissue_Genotype)
lowcov_rel <- as.data.frame(as.matrix(sample_data(lowcov))) %>%
  select(Tissue, genotype, stage_score, Sample_Depth_alpha)
lowcov_rel[order(lowcov_rel$Sample_Depth_alpha),]


#ITS depths: vlow: 751, low: 1010
#16S depths: 1380 to not lose all of the lyk4 flowers,1004 to improve losses of fls2 silinew 

######################################################
#######  Selected depths: 16S = 1380, ITS = 751 ###### 
######################################################





###############################################
## Repeat rarefy + calculate Alpha diversity ##
###############################################



#Repeat rarefy for OTU table + Alpha Div metrics (Shannon, Pielou, Faith's)

#Repeat rarefy for alpha diversity metrics
# 1) Rarefy table
#2) Calculate Alpha diversity metrics 
#3) Combine alpha diversity metrics into table and take the mean
#Set for 16S

# Import weighted.faith function by Nate Swenson lefse_0.5 https://github.com/NGSwenson/lefse_0.5/blob/master/R/weighted.faith.R
source("weighted.faith_function.R") 

df_alpha <- data.frame() #Create empty dataframe for alpha diversity measures
seed_set <- sample.int(1e6, 100) # generate set of seeds for sampling

#Calculate alpha diversity metrics over 100 iterations
for (y in seed_set){
  #Rarefy phyloseq table
  x <- rarefy_even_depth(exp_only, sample.size=1380, trimOTUs =FALSE, #Set sample.size to rarefy depth of choice
                         replace = FALSE, verbose = TRUE, rngseed = y) 
  #extract Shannon and Pielou's evenness
  alpha <- microbiome::alpha(x, index = c("diversity_shannon", "evenness_pielou"))
  shannon <- alpha$diversity_shannon #Shannon diversity
  pielou <-alpha$evenness_pielou # Pielou evenness
  md <- sample_data(x)
  
  #Faith's phylogenetic distance in picante 
  #NOTE: LONG STEP! DO NOT RUN FOR ITS
  phy <- phy_tree(x)
  comm <- as.data.frame(t(otu_table(x)))
  faith <- pd(comm, phy)
  raw_faith <- faith$PD
  Alpha_adj_faith = raw_faith/num_species
  num_species <- faith$SR
  #Create alpha diversity dataframe with the iterations alpha diversity measures
  x_alpha <- cbind(md, Shannon = shannon, Pielou = pielou, Faith = raw_faith, 
                   Alpha_adj_faith = Alpha_adj_faith, Number_of_species = num_species)
  #Add each iteration of alpha diversity to a dataframe. 
  df_alpha <- rbind(df_alpha, x_alpha)
}

#Save alpha diveristy dataframe
write.table(df_alpha, file = "Alpha_Diversity_repeat_rarefy_16S.csv", sep = ",")


################
## Statistics ##
################


#Prep data frame for ANOVA
df_alpha <- read.table("Alpha_Diversity_repeat_rarefy_16S.csv",
                       sep = ",", header = TRUE, fill = TRUE)

df_alpha <- subset(df_alpha, Tissue != "soil") #Remove soil from analysis

#Create dataframe that calculates the mean and standard deviation of alpha diversity for each sample
STS_reprare_means <- df_alpha %>%
  group_by(Derep_ID) %>%
  summarize(Nice_Stage = Nice_Stage, Nice_Tissue = Nice_Tissue, genotype = genotype, Stage_Tissue = Stage_Tissue, #retain important metadata
            Shannon_mean = mean(Shannon), Shannon_sd = sd(Shannon), #Shannon diversity mean and standard deviation
            Pielou_mean = mean(Pielou), Pielou_sd = sd(Pielou), #Pielou's evenness mean and standard deviation
            Faith_mean = mean(Faith), Faith_sd = sd(Faith), #Faith's evenness mean and standard deviation
            Alpha_adj_faith_mean = mean(Alpha_adj_faith), Alpha_adj_faith_sd = sd(Alpha_adj_faith),
            Richness_mean = mean(Number_of_species)
  )

STS_reprare_uni <- distinct(STS_reprare_means) #Remove duplicate rows

STS_rru_po <- subset(STS_reprare_uni, Nice_Tissue != "Soil") #Make sure that soil is removed


#Baseline ANOVAs 
#For explanatory factor, can substitute Pielou_mean (Pielou's evenness) with Shannon_mean (Shannon Diversity) and/or Alpha_adj_faith_mean (richness adjusted Faith's Phylogenetic diversity)
my_anova <- aov(Pielou_mean ~ Nice_Tissue * Nice_Stage * genotype, data = STS_rru_po) #Run ANOVA with stage, tissue and genotype as interacting factors
summary(my_anova)
my_anova <- aov(Pielou_mean ~ Stage_Tissue * genotype, data = STS_rru_po) #Run ANOVA with each Stage and Tissue combinations interacting with genotype as a factor (no missing groups)
summary(my_anova)

#relevant subsets, to make sure no missing groups in ANOVA
po_rr <- subset(STS_rru_po, Nice_Tissue == "Root" | Nice_Tissue == "Rosette") #only roots and rosettes across all time
po_r <-  subset(STS_rru_po, Nice_Tissue != "Siliques" & Nice_Stage != "Vegetative") #three stages, no siliques (Vegetative only has roots and rosettes, Flowering doesn't have siliques)
po_reprod <- subset(STS_rru_po, Nice_Stage == "Unripe Siliques" | Nice_Stage == "Ripe Siliques") #only final two stages, all tissue types present

#Base ANOVA for subsets 
#subsitute diversity metric (Faith_mean) as above, also change subset in data section to reflect evaluated subset.
my_anova <- aov(Faith_mean ~ Nice_Tissue * Nice_Stage * genotype, data = po_rr) 
summary(my_anova)

#ANOVA permutations for p values as per Manly (2007), based on scripts by David Howell (https://www.uvm.edu/~statdhtx/StatPages/Permutation%20Anova/PermTestsAnova.html)

alpha_frame <- STS_rru_po #define dataframe to run anova on 

#create specific dataframe with relevant data for permuations. for the div_col value, substitute "alpha_frame$Pielou_mean" with "alpha_frame$Shannon_mean" and/or "alpha_frame$Faith_mean"
af <- data.frame(div_col <-  alpha_frame$Pielou_mean,
                 stages <-  alpha_frame$Nice_Stage,
                 genos <-  alpha_frame$genotype,
                 tissues <- alpha_frame$Nice_Tissue)

colnames(af) <- c("div_col", "stages", "genos", "tissues") #Name columns in anova dataframe
af <- subset(af, tissues != "soil") #double check that soil is removed


# Create standard Anova on these data
#Define each factor and interaction F values
mod1 <- lm(div_col ~tissues* stages * genos, data = af) #run anova
ANOVA <- summary(aov(mod1)) #store results
print(ANOVA) 
#Saving F of real data value for future use, selects in order of apperance in the ANOVA table
Ftissue <-  ANOVA[[1]]$"F value"[1] #selects F value for first factor on ANOVA table (tissue alone)
Fstage <-  ANOVA[[1]]$"F value"[2] #selects F value for 2nd factor on ANOVA table (stage alone)
Fgeno <-  ANOVA[[1]]$"F value"[3]
Fintts <- ANOVA[[1]]$"F value"[4]
Fintgt <-  ANOVA[[1]]$"F value"[5]
Fintgs <-  ANOVA[[1]]$"F value"[6]
Finttsg <- ANOVA[[1]]$"F value"[7]

# Now start resampling
nreps <- 999 #define number of resampling iterations
#Set up space to store F values for each factor/interaction as calculated .
FT <- numeric(nreps)   
FS <- numeric(nreps)  
FG <- numeric(nreps)
FTS <- numeric(nreps)
FGT <- numeric(nreps)
FGS <- numeric(nreps)
FTSG <- numeric(nreps)

# The first F of 999
FT[1] <- Ftissue          
FS[1] <- Fstage
FG[1] <- Fgeno
FTS[1] <- Fintts
FGT[1] <- Fintgt
FGS[1] <- Fintgs
FTSG[1] <- Finttsg

#FREE PERMUTATIONS. Permuting loop.
for (i in 2:nreps) {
  
  new_div <- sample(af$div_col, length(af$div_col)) #Shuffle diversity scores from dataset randomly
  af$new_div <- new_div #add this shuffled score into the dataframe
  
  mod2 <- lm(new_div ~ tissues* stages * genos, data = af) #run ANOVA with the shuffled diversity metrics
  b <- summary(aov(mod2)) #Store ANOVA values
  #Store F values from shuffled ANOVA
  FT[i] <- b[[1]]$"F value"[1]
  FS[i] <- b[[1]]$"F value"[2]
  FG[i] <- b[[1]]$"F value"[3]
  FTS[i] <-  b[[1]]$"F value"[4]
  FGT[i] <- b[[1]]$"F value"[5]
  FGS[i]<- b[[1]]$"F value"[6]
  FTSG[i]<- b[[1]]$"F value"[7]
}

#Calculate the probability that the resampled/perumted F values exceed the actual F value from the dataset
probT <- length(FT[FT >= Ftissue+ .Machine$double.eps ^0.5])/nreps #.Machine$double.eps ^0.5 adds a small pseudocount
probS <- length(FS[FS >= Fstage + .Machine$double.eps ^0.5])/nreps       
probG  <-  length(FG[FG >= Fgeno + .Machine$double.eps ^0.5])/nreps
probTS <- length(FTS[FTS >= Fintts + .Machine$double.eps ^0.5])/nreps
probGT <- length(FGT[FGT >= Fintgt + .Machine$double.eps ^0.5])/nreps
probGS <- length(FGS[FGS >= Fintgs + .Machine$double.eps ^0.5])/nreps
probTSG <- length(FTSG[FTSG >= Finttsg + .Machine$double.eps ^0.5])/nreps

#View perumuted p values and manually add to ANOVA table. 
probT
probS
probG
probTS
probGT
probGS
probTSG





################
## Graphics ##
################

alpha_po <- STS_reprare_means

alpha_po$genotype <- as.character(alpha_po$genotype) #Change genotype class

#Alpha div graphics: overall

alpha_po <- alpha_po %>% 
  mutate(Nice_Genotype = case_when(genotype == "col0" ~ "WT",
                                   TRUE ~ genotype)) #Edit col0 coding to show WT

alpha_po$Nice_Genotype <- factor(alpha_po$Nice_Genotype, levels = c("WT", "efr", "fls2", "lore", "lyk4")) #Set order for genotype
alpha_po$Nice_Tissue <- factor(alpha_po$Nice_Tissue, 
                               levels = c("Root", "Rosette", "Stem", "Cauline Leaves", "Flowers", "Siliques")) #Set graphing order for Tissues

alpha_po$Nice_Stage <- factor(alpha_po$Nice_Stage, levels = c("Vegetative", "Flowering", "Unripe Siliques", "Ripe Siliques")) #Set graphing order for Stages


genotype_labels <- c("WT", "*efr*", "*fls2*", "*lore*", "*lyk4*") #set labels for genotype (italics)



div_plot <- ggplot(alpha_po, aes(x = Nice_Genotype, 
                                  y = Alpha_adj_faith_mean, #Change to reflect diversity metric of choice
                                  fill = Nice_Genotype
)) + 
  geom_boxplot() + #create box plot
  geom_jitter(cex = 0.2) + # add data points
  #facet_wrap(~Nice_Stage) + #to view genotype by stage
  #facet_grid(rows = vars(Nice_Stage), cols = vars(Nice_Tissue), switch = "y") + #to view genotype by stage by tissue
  
  #formatting
  theme_bw() + 
  xlab("Genotype") + 
  ylab("Faith's PD \n (scaled by richness)") +
  scale_colour_viridis_d()+
  scale_fill_viridis_d() + 
  scale_x_discrete(labels = genotype_labels)+
  theme(legend.position="none",
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        axis.text.x = element_markdown()) +
  
  #y axis limits
  #ylim(0,6.25)+ #Shannon
  #ylim(0,1.1) + #Pielou
  ylim(0,0.5) + #Alpha adj PD
  
  #manually add significance annotations
  geom_signif(y_position = c(0.48), xmin = c(0.7), xmax = c(5.3), annotation = c("NS"), tip_length = 0)

div_plot #View plot


#combined 16S and ITS plot
alpha_plots <- ggarrange(plotlist = list(div_plot, div_plot_ITS), labels = c("A", "B"))

alpha_plots
