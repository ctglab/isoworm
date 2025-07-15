### Create and load vector of packages
setwd("/home/runner/work/isoworm/isoworm/workflow/R")
source("config_all.R")
lapply(my_packages, require, character.only = TRUE) 
##### open all the bed files ##################################################
samples_dir <- list.dirs(polyA_bam_dir, recursive = FALSE)
dataframes_peak <- list()
for (i in samples_dir){
  setwd(i)
  in.list <- list.files(path= i, recursive = T, full.names = T)
  j <- grep("_Aligned.sortedByCoord.out.bam$",in.list)
  name <- paste0("quant_",gsub(".*/(SRR[0-9]+).*", "\\1",in.list[j]))
  gal <- readGAlignments(in.list[j])
  coverage <- coverage(gal)
  chr7_peaks <- slice(coverage$chr7, lower = 1)
  cutoff <- chr7_peaks[(start(chr7_peaks) > 140719336 - 50) & (end(chr7_peaks) < 140719437 + 50)]
  if (length(cutoff) == 1) { 
    cutoff <- median(cutoff) * 0.05
    chr7_peaks <- slice(coverage$chr7, lower = cutoff)
    filter_80 <- subset(chr7_peaks, width(chr7_peaks) >= 80)
    dataframes_peak[[name]] <- assign(name, filter_80)
  }
}

final_list_transcript_1 <-list()
final_list_transcript_2  <- list()
final_list_transcript  <- list()
for (i in 1:length(dataframes_peak)){
  name <- paste0("peaks_",names(dataframes_peak)[[i]])
  final_list_transcript_1[[name]]  <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_unique_region_transcript_1) & (end(dataframes_peak[[i]]) <= start_unique_region_transcript_1)))
  final_list_transcript_2[[name]]   <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_unique_region_transcript_2) & (end(dataframes_peak[[i]]) <= start_unique_region_transcript_2)))
  final_list_transcript[[name]] <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_transcript) & (end(dataframes_peak[[i]]) <= start_transcript)))
}

# Crea un dataframe con le informazioni estratte
bed_list_transcript_1 <-list()
for (i in 1:length(final_list_transcript_1)){
  if (length(final_list_transcript_1[[i]]) == 0) {
    next  # Skip this iteration
  }
  name <- paste0("bed_", names(final_list_transcript_1)[[i]])
  bed_list_transcript_1[[name]] <- assign(name, data.frame(seqnames = "chr7", starts = start(final_list_transcript_1[[i]]), ends = end(final_list_transcript_1[[i]])))
}

bed_list_transcript_2 <-list()
for (i in 1:length(final_list_transcript_2)){
  if (length(final_list_transcript_2[[i]]) == 0) {
    next  # Skip this iteration
  }
  name <- paste0("bed_", names(final_list_transcript_2)[[i]])
  bed_list_transcript_2[[name]] <- assign(name, data.frame(seqnames = "chr7", starts = start(final_list_transcript_2[[i]]), ends = end(final_list_transcript_2[[i]])))
}

final_polyA_transcript_1 <- bedtoolsr::bt.multiinter(bed_list_transcript_1)
print(final_polyA_transcript_1)
final_polyA_transcript_1
final_polyA_filtered_transcript_1 <- subset(final_polyA_transcript_1, (V3 - V2) >= 30 & V4 >= 50)
final_polyA_filtered_transcript_1 <- subset(final_polyA_filtered_transcript_1, select = 1:5)
colnames(final_polyA_filtered_transcript_1) <- c("chr7","start","end","total of samples","samples names")

final_polyA_transcript_2  <- bedtoolsr::bt.multiinter(bed_list_transcript_2)
final_polyA_filtered_transcript_2 <- subset(final_polyA_transcript_2, (V3 - V2) >= 30 & V4 >= 50)
final_polyA_filtered_transcript_2 <- subset(final_polyA_filtered_transcript_2, select = 1:5)
colnames(final_polyA_filtered_transcript_2) <- c("chr7","start","end","total of samples","samples names")

setwd(final_output)
write.csv(x=final_polyA_filtered_transcript_1, file="polyA_filtered_unique_regiontranscript_1.csv", row.names = FALSE)
write.csv(x=final_polyA_filtered_transcript_2,  file="polyA_filtered_unique_regiontranscript_2.csv", row.names = FALSE)
