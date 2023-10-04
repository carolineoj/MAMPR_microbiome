### Cutadapt ###

import os
import sys
import re
import glob


#Edit command_part1 based on 16S or ITS

#ITS
#command_part1 =  '~/.local/bin/cutadapt -e 0.15 -a "CTTGGTCATTTAGAGGAAGTAA;required...GCATCGATGAAGAACGCAGC" -A "GCTGCGTTCTTCATCGATGC;required...TTACTTCCTCTAAATGACCAAG" -m 100:100'

#16S
command_part1 = '~/.local/bin/cutadapt -e 0.2 -g AACMGGATTAGATACCCKG -G ACGTCATCCCCACCTTCC -m 100:100' 

#Set input and output paths

folderpath = 'ITS_FASTQs/MAMPRfield_ITS_FASTQs_P11-P14' #input folder

outputfolder = 'ITS_final_analysis/Cutadapt_trimmed' #trimmed output folder

untrimmedfolder = 'ITS_final_analysis/Cutadapt_untrimmed' #untrimmed output folder


os.chdir(folderpath)           #Move into folder
R1_files = glob.iglob('*_R1_001.fastq.gz') #Select all read 1 files
for R1_file in R1_files:
    R2_file = R1_file.partition("_R")[0]+'_R2_001.fastq.gz' #select paired read 2 file
    print('Running the following command:\n '+command_part1+' -o '+outputfolder+'/trimmed_'+R1_file+' -p '+outputfolder+'/trimmed_'+R2_file+' --untrimmed-output '+untrimmedfolder+'/untrimmed_'+R1_file+' --untrimmed-paired-output '+untrimmedfolder+'/untrimmed_'+R2_file+' '+R1_file+' '+R2_file+'')
    os.system(command_part1+' -o '+outputfolder+'/trimmed_'+R1_file+' -p '+outputfolder+'/trimmed_'+R2_file+' --untrimmed-output '+untrimmedfolder+'/untrimmed_'+R1_file+' --untrimmed-paired-output '+untrimmedfolder+'/untrimmed_'+R2_file+' '+R1_file+' '+R2_file+'\n')
# exit()