#############################################################
############### MICROBIAL LOAD AND FITNESS ##################
#############################################################



#########################
######### LOAD ##########
#########################


#Load analysis: use ps_kingdom_L phyloseq object generated in MAMPR_3 script
ps_kingdom_L

#Adjust total load THIS DEPENDS ON IF ps_kingdom_L IS AN ITS OR 16S BASED PHYLOSEQ OBJECT
sample_data(ps_kingdom_L)$ITS_sample_sums_adj <- sample_data(ps_kingdom_L)$ITS_factor_adj * sample_sums(ps_kingdom_L)
#OR
sample_data(ps_kingdom_L)$STS_sample_sums_adj <- sample_data(ps_kingdom_L)$STS_factor_adj * sample_sums(ps_kingdom_L)

#Get dataframe with raw microbial load data
load_df <- sample_data(ps_kingdom_L)

# Plant only
po <- subset(load_df, Tissue != "NC" & Tissue != "Mock" & Tissue != "silisen"  & Tissue != "soil")
df <- as.data.frame(as.matrix((po)))

#Make sure load data is numeric
df$ITS_sample_sums_adj <- as.numeric(df$ITS_sample_sums_adj) 
df$STS_sample_sums_adj <- as.numeric(df$STS_sample_sums_adj) 




#
##
### Statistics

multivariate_load <- aov(log10(STS_sample_sums_adj) ~ Tissue*stage_score*genotype, 
                         data=df)
summary(multivariate_load)




#
##
### Graphics
load_choice <- df

#Graphs: Overall + subsets
load_choice <- load_choice %>% 
  mutate(Nice_Genotype = case_when(genotype == "col0" ~ "WT",
                                   TRUE ~ genotype)) #change col0 label to wild type label

#Set order of factors for graphing
load_choice$Nice_Genotype <- factor(load_choice$Nice_Genotype, levels = c("WT", "efr", "fls2", "lore", "lyk4"))
load_choice$Nice_Tissue <- factor(load_choice$Nice_Tissue, 
                                  levels = c("Root", "Rosette", "Stem", "Cauline Leaves", "Flowers", "Siliques"))
load_choice$Nice_Stage <- factor(load_choice$Nice_Stage, levels = c("Vegetative", "Flowering", "Unripe Siliques", "Ripe Siliques"))


genotype_labels <- c("WT", "*efr*", "*fls2*", "*lore*", "*lyk4*") #set labels for genotype (italics)

div_plot <- ggplot(load_choice, aes(x = genotype, 
                                    y = log10(STS_sample_sums_adj), #change for ITS dataset
                                    fill = genotype)) + 
  geom_boxplot() +
  theme_bw() + 
  geom_jitter(cex = 0.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  #facet_wrap(~Nice_Tissue)+ #option to subset by Tissue
  #facet_grid(vars(Nice_Stage), vars(Nice_Tissue)) + #option to subset by tissue x stage
  theme(legend.position="right") + 
  xlab("Genotype") + 
  ylab("log10 Spike \n Normalized Reads") +
  scale_x_discrete(
    labels = genotype_labels)+
  theme(legend.position="none",
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size=14),
        axis.text.x = element_markdown()) +
  ylim(3,5) +
  geom_signif(
    y_position = c(4.9), xmin = c(0.7), xmax = c(5.3),
    annotation = c("NS"), tip_length = 0)
div_plot

#Combine plots from 16S and ITS
alpha_plots <- ggarrange(plotlist = list(div_plot, div_plot_ITS),
                         labels = c("A", "B"))

alpha_plots







#######################
####### FITNESS #######
#######################


#Fitness information is contained in the metadata file

fitness_df <- subset(metadata, Tissue == "rose") #Only select one tissue per individual plant to avoid falsely inflating sampling

#Set orders, numeric, and factorize factors
fitness_df$Nice_Stage <- factor(fitness_df$Nice_Stage, levels = c("Vegetative", "Flowering", "Unripe Siliques", "Ripe Siliques"))
fitness_df$genotype <- factor(fitness_df$genotype, levels = c("col0", "efr", "fls2", "lore", "lyk4"))
fitness_df$Total_siliques <- as.numeric(fitness_df$Total_siliques)
fitness_df$Final_weight_mg <- as.numeric(fitness_df$Final_weight_mg)
fitness_df$weeks_old_f <- factor(fitness_df$weeks_old)
fitness_df$weeks_old <- as.numeric(fitness_df$weeks_old)
fitness_df$Day_old <- as.numeric(fitness_df$Day_old)
fitness_df <- subset(fitness_df, Final_weight_mg > 0)

fitness_df <- fitness_df %>%
  mutate(wt_mutant = case_when(genotype == "col0" ~ "WT",
                               TRUE ~ "mutant"))


#
##
### Statistics 
#Get sample size


#All samples
#look for genotype effect within stage

#Fitness by rosette dry weight (all stages: test "veg", "flow_only", "young_siliques" and "mature_siliques")
#alternatively could test by weeks old value
kruskal.test(Final_weight_mg ~ genotype, subset(fitness_df, stage_score == "veg" )) 

#Fitness by total silique count (two stages: test "young_siliques" and "mature_siliques")
#alternatively could test by weeks old value
kruskal.test(log(Total_siliques) ~ genotype, subset(fitness_df, stage_score == "mature_siliques"))


#difference between stages?

#rosette dry weight change by stage 
#alternatively could test by weeks old value
kruskal.test(Final_weight_mg ~ stage_score, fitness_df) 

#silique count by stage (two stages: test "young_siliques" and "mature_siliques")
#alternatively could test by weeks old value
kruskal.test(log(Total_siliques) ~ stage_score, subset(fitness_df, stage_score == "young_siliques" | stage_score == "mature_siliques"))



#
##
### Graphics

#create list for significance comparisions 
comparison_list <- list(c("Unripe Siliques", "Ripe Siliques"))


#Silique count box plots
silique_boxplots <- ggplot(subset(fitness_df, Nice_Stage == "Unripe Siliques" | Nice_Stage == "Ripe Siliques"),
                           aes(x=Nice_Stage, y= log10(Total_siliques), fill = genotype)) + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), aes(x=Nice_Stage, y= log10(Total_siliques), fill = genotype), size = 1) +
  xlab("Developmental stage") +
  ylab("Silique count (log10)") +
  theme_classic() +
  #For viridis color pallette
  scale_fill_discrete(labels = genotype_labels,
                      type = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
    theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_markdown(size = 12))+
  labs(fill = "Genotype") +
  geom_signif(comparisons = comparison_list,
              map_signif_level = TRUE) +
  #log transformed signficance levels
  geom_signif(
    y_position = c(3.5, 3.5), xmin = c(0.7, 1.7), xmax = c(1.3, 2.3),
    annotation = c("NS", "NS"), tip_length = 0)

silique_boxplots


#rosette weight box plots
dryweight_boxplots <- ggplot(fitness_df, aes(x=Nice_Stage, y= Final_weight_mg, fill = genotype)) + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), aes(x=Nice_Stage, y= Final_weight_mg, fill = genotype), size = 1) +
  xlab("Developmental stage") +
  ylab("Rosette dry weight (mg)") +
  theme_classic() +
  #For viridis color pallette
  scale_fill_discrete(labels = genotype_labels,
                      type = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_markdown(size = 12))+
  labs(fill = "Genotype") +
  geom_signif(
    y_position = c(32, 50, 50, 50), xmin = c(0.7, 1.7, 2.7, 3.7), xmax = c(1.3, 2.3, 3.3, 4.3),
    annotation = c("NS", "NS", "NS", "NS"), tip_length = 0) 

dryweight_boxplots

#combine plots for fitness
fitness_plots <- ggarrange(plotlist = list(silique_boxplots, dryweight_boxplots),
                           labels = c("A", "B"),
                           common.legend = TRUE, legend = "right", widths = c(0.6, 1))

fitness_plots

a_fitness <- annotate_figure(fitness_plots, top = text_grob("Genotype effect on fitness", size = 18))
a_fitness








#Correlation between rosette weight and silique count
#Graphics and statistics
rose_sili_cor <- ggplot(data = subset(fitness_df, stage_score == "young_siliques" | stage_score == "mature_siliques"), 
                        aes(x=Final_weight_mg, y=log10(Total_siliques), color = weeks_old_f, group = weeks_old_f)) + 
  geom_point() +
  theme_classic() +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor(method = "pearson") + #pearson correlation statistics 
  scale_shape_manual(values=c(3, 5, 7)) + 
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  labs(color = "Age (weeks)") +
  xlab("Rosette dry weight") +
  ylab("Silique count (log10)")
rose_sili_cor

#DATA CORRELATIONS FOR LOG10 SIIQUES AND ROSETTE DRY WEIGHT
# WEEK 31: R = 0.55, p = 0.0017
# WEEK 32: R = 0.81, p < 0.001
# WEEK 33: R = 0.67, p = 0.07




