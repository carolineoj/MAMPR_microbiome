# MAMPR_microbiome
Scripts used in MAMP-detecting PRRs have little effect on endophytic microbiome assembly in A. thaliana

Various scripts used to analyze the endophytic microbiomes of A. thaliana PRR mutants in the field. The order of the scripts is linear; later scripts build on earlier scripts.

Order of scripts:

Preprocessing
1) MAMPR_cutadapt.py #uses cutadapt to remove primers. In python.
2) MAMPR_dada2.R #Run DADA2 to denoise and infer amplicon sequence variants (ASVs). In R.
3) MAMPR_QIIME2_spike_processing_taxonomy.txt #Processes sequences by identifying and removing spikes and classifying sequences. In QIIME2 command line environment

Analysis (Almost all in R, some partially in QIIME2 command line environment)
1) MAMPR_1_create_phyloseq_min_filtering.R #Filter and process data set
2) MAMPR_2_alpha_diversity.R #Rarefaction and alpha diversity analysis
3) MAMPR_3_beta_diversity_and_ dispersion.R #Core filtering, beta diversity and dispersion
4) MAMPR_4_differential_abundance.R #Differential abundance
5) MAMPR_5_homogenous_plants.R #Degree of within-plant similarity of different tissues.
6) MAMPR_6_load_and_fitness.R #Microbial load (estimated by spike ratio) and early fitness analysis




