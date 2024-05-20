### Create and load vector of packages
setwd("/Volumes/HD2/home/maurizio/R/finale")
source("config_all.R")
lapply(my_packages, require, character.only = TRUE) 
##### open all the bed files ##################################################
samples_dir <- list.dirs(polyA_bam_dir, recursive = FALSE)
dataframes_peak <- list()
for (i in samples_dir){
  setwd(i)
  in.list <- list.files(path= i, recursive = T, full.names = T)
  j <- grep("_small_Aligned.sortedByCoord.out.bam$",in.list)
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

final_list_ref <-list()
final_list_x1  <- list()
final_list_BRAF  <- list()
for (i in 1:length(dataframes_peak)){
  name <- paste0("peaks_",names(dataframes_peak)[[i]])
  final_list_ref[[name]]  <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_3utr_ref) & (end(dataframes_peak[[i]]) <= start_3utr_ref)))
  final_list_x1[[name]]   <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_3utr_X1) & (end(dataframes_peak[[i]]) <= start_3utr_X1)))
  final_list_BRAF[[name]] <- assign(name, subset(dataframes_peak[[i]], (start(dataframes_peak[[i]]) >= end_BRAF) & (end(dataframes_peak[[i]]) <= start_BRAF)))
}

# Crea un dataframe con le informazioni estratte
bed_list_ref <-list()
for (i in 1:length(final_list_ref)){
  if (length(final_list_ref[[i]]) == 0) {
    next  # Skip this iteration
  }
  name <- paste0("bed_", names(final_list_ref)[[i]])
  bed_list_ref[[name]] <- assign(name, data.frame(seqnames = "chr7", starts = start(final_list_ref[[i]]), ends = end(final_list_ref[[i]])))
}

bed_list_x1 <-list()
for (i in 1:length(final_list_x1)){
  if (length(final_list_x1[[i]]) == 0) {
    next  # Skip this iteration
  }
  name <- paste0("bed_", names(final_list_x1)[[i]])
  bed_list_x1[[name]] <- assign(name, data.frame(seqnames = "chr7", starts = start(final_list_x1[[i]]), ends = end(final_list_x1[[i]])))
}

final_polyA_ref <- bedtoolsr::bt.multiinter(bed_list_ref)
final_polyA_filtered_ref <- subset(final_polyA_ref, (V3 - V2) >= 30 & V4 >= 50)
final_polyA_filtered_ref <- subset(final_polyA_filtered_ref, select = 1:5)
colnames(final_polyA_filtered_ref) <- c("chr7","start","end","total of samples","samples names")

final_polyA_x1  <- bedtoolsr::bt.multiinter(bed_list_x1)
final_polyA_filtered_x1 <- subset(final_polyA_x1, (V3 - V2) >= 30 & V4 >= 50)
final_polyA_filtered_x1 <- subset(final_polyA_filtered_x1, select = 1:5)
colnames(final_polyA_filtered_x1) <- c("chr7","start","end","total of samples","samples names")

setwd(final_output)
write.csv(x=final_polyA_filtered_ref, file="polyA_filtered_3UTR220.csv", row.names = FALSE)
write.csv(x=final_polyA_filtered_x1, file="polyA_filtered_3UTR204.csv", row.names = FALSE)
