library(dada2)
library(purrr)
library(stringr)

path <- "Cutadapt_trimmed/P11-P14"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

plotQualityProfile(fnFs[1:2])

#Filter and trim
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), #truncLen=c(216,149), #use truncLen for 16S based on FIGARO recomendations (different for each run). No truncating for ITS
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

#Learn error rates

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)



#merge paired end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)



#Save as an RDS object
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "ITS_P11-P14_seqtab.rds") # CHANGE ME to where you want sequence table saved

### Repeat above procces for each run


#load and merge sequence tables
st1 <- readRDS("ITS_final_analysis/ITS_P1-P5_seqtab.rds")
st2 <- readRDS("ITS_final_analysis/ITS_P6-P10_seqtab.rds")
st3 <- readRDS("ITS_final_analysis/ITS_P11-P14_seqtab.rds")
st.all <- mergeSequenceTables(st1, st2, st3)

#Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
st_copy <- seqtab

cnm <- collapseNoMismatch(seqtab, minOverlap = 50)
cnm_copy <- cnm

table(nchar(getSequences(seqtab)))

#Edit sample names

target_names <- str_split(row.names(cnm), "_") %>% map_chr(`[`, 2)
row.names(cnm) <- target_names

#Write tables
#Frequency ASV table
write.table(t(cnm), "ITS_final_analysis/ITS_P1-P14_seqtab_cnm_name_correct_dada2.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#Fastas
uniquesToFasta(cnm, fout='ITS_final_analysis/ITS_P1-P14_seqtab_cnm_name_correct_dada2.fna', ids=colnames(cnm))



