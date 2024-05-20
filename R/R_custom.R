### Create and load vector of packages
setwd(final_output)
source("config_all.R")
lapply(my_packages, require, character.only = TRUE) 

### Create a list with all txt files containing the samples for each tissues
### Extract all samples from the txt list files for each tissues
setwd(txt_files_dir)
txt_files      <- list()
samples_files  <- list()
count = 1
for(j in samples_typologies){
  txt_files     <- c(txt_files, list(assign(paste0(j, "_files"), list.files(full.names = F , recursive =F, pattern= j))))
  samples_files <- c(samples_files,list(assign(paste0(j, "_samples"),lapply(txt_files[[count]],readLines))))
  count = count + 1 
}

files_results_dir <- list.dirs(results_dir, recursive = FALSE)
all_results_dir   <- list.dirs(files_results_dir, recursive = FALSE)
bam_samples       <- all_results_dir[grepl("/bam",all_results_dir, ignore.case = TRUE)]
ballgown_samples       <- all_results_dir[grepl("/ballgown",all_results_dir, ignore.case = TRUE)]
bam_directories <- bam_samples[sapply(bam_samples, function(x) any(grepl(paste(samples_typologies, collapse = "|"), x)))]
ballgown_directories <- ballgown_samples[sapply(ballgown_samples, function(x) any(grepl(paste(samples_typologies, collapse = "|"), x)))]

## Create a list with all the total reads for each samples inside the results directory
library_reads_number <- list()
for (y in 1:length(samples_files)){
  library_reads_number <- list()
  for (i in 1:length(samples_files[[y]])){
    for(j in samples_files[[y]][[i]]){
      setwd(bam_directories[y])
      current_sample <- j
      in.list <- list.files(path= j, recursive = T, full.names = T)
      library_reads  <- read.csv(in.list[1],sep="\t")[c(4),]
      library_reads_number <-  na.omit(append(library_reads_number,as.numeric(library_reads[2][,1])))
      assign(paste0(samples_typologies[y],"_reads_counts"), library_reads_number)
    }
  }
} 

### Extract reads count from stringtie final file for REF and X1 (fragments) for each cell line typology
## Create a list with all the count reads for each samples (X1 and ref) inside the results directory  (for each group)
read_count_x1 <- list()
read_count_ref <- list()
for (y in 1:length(samples_files)){
  read_count_x1  <- list()
  read_count_ref <- list()
  for (i in 1:length(samples_files[[y]])){
    for(j in samples_files[[y]][[i]]){
      setwd(ballgown_directories[y])
      current_sample <- j
      in.list <- list.files(path= j, recursive = T, full.names = T)
      exons_data <- fread(in.list[1])
      read_count_x1  <- na.omit(append(read_count_x1, as.numeric(exons_data$rcount[1])))
      read_count_ref <- na.omit(append(read_count_ref,as.numeric(exons_data$rcount[3])))
      assign(paste0(samples_typologies[y],"_X1_counts"),  read_count_x1)
      assign(paste0(samples_typologies[y],"_ref_counts"), read_count_ref)
    }
  }
} 

### Calculate FPKM for X1 and Ref for each tissue typology
FPKM_Ref <-list()
FPKM_X1  <-list()
count = 0
for (y in 1:length(samples_files)){
  FPKM_Ref <-list()
  FPKM_X1  <-list()
  count = 0
  for (i in 1:length(samples_files[[y]])){
    for(j in 1:length(samples_files[[y]][[i]])){
      count = count + 1
      FPKM_X1  <- append(FPKM_X1, 10^9* get(paste0(samples_typologies[y],"_X1_counts"))[[count]]/(length_x1*get(paste0(samples_typologies[y],"_reads_counts"))[[count]]))
      FPKM_Ref <- append(FPKM_Ref,10^9* get(paste0(samples_typologies[y],"_ref_counts"))[[count]]/(length_ref*get(paste0(samples_typologies[y],"_reads_counts"))[[count]]))
      assign(paste0(samples_typologies[y],"_X1_FPKM"),  FPKM_X1)
      assign(paste0(samples_typologies[y],"_ref_FPKM"), FPKM_Ref)  
    }
  }
}

### Create the BRAF expression list of dataframe tissues
BRAF_results <- list()
count = 1
for(i in samples_typologies){
  assign(paste0(i,"_BRAF_results"), do.call(rbind, Map(data.frame, BRAF_ref=get(paste0(i,"_ref_FPKM")), BRAF_X1=get(paste0(i,"_X1_FPKM")))))
  BRAF_results[count] <- list(get(paste0(i,"_BRAF_results")))
  colnames(BRAF_results[[count]]) <- c("BRAF Reference", "BRAF X1")
  count = count + 1 
}
db_ratio <- BRAF_results
for(i in 1:length(label_plots)){
  colnames(db_ratio[[i]]) <- c("BRAF Reference", "BRAF X1")
  db_ratio[[i]]$ratio <- log2(((db_ratio[[i]]$`BRAF X1`+0.01)/(db_ratio[[i]]$`BRAF Reference`+0.01)))
  db_ratio[[i]]       <- db_ratio[[i]] %>% 
    pivot_longer(
      cols = `BRAF Reference`:ratio,
      names_to = "Isoforms",
      values_to = "value"  
    )
}

db_ratio <- rbindlist(db_ratio)
db_ratio <- subset(db_ratio,Isoforms=="ratio")
db_ratio$group <- NA
db_ratio$IDS   <- NA
tissues_numbers  <- c()
all_objects <- ls()
pattern<- "_BRAF_results"
varnames <- grep(pattern, all_objects, value = TRUE)
sapply(varnames, get)

for (i in 1:length(label_plots)){
  tissues_numbers  <- c(tissues_numbers,as.integer(count(get(varnames[i]))))
}

count  = 1
count2 = tissues_numbers[1]
for (i in 1:length(label_plots)){
  db_ratio$group[count:count2] <- label_plots[i]
  db_ratio$IDS[count:count2] <- unlist(samples_files[[i]])
  count  = tissues_numbers[i] + count
  count2 = tissues_numbers[i+1] + count -1
}


# Change  automatically color by groups
# Opening the graphical device
bp_ratios <- ggplot(db_ratio, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() +
  #scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  ylim(0,10) +
  geom_jitter(position=position_jitter(0.2), size=2.0, alpha=0.9, pch=21, fill = "black") +
  theme_bw() +
  theme(text=element_text(size=28, face = "bold"),
        axis.text=element_text(size=22, face = "bold"),
        axis.title.x =  element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  #scale_fill_manual(values=group.colors_plot) +
  ylab(expression(paste(italic("BRAF-204/BRAF-220"), " (IsoWorm Custom)"))) +
  xlab("Tissue Typology") 
bp_ratios
# Closing the graphical device
setwd(final_output)
pdf("ratio_BRAF.pdf", width = 21.5,height = 11)
bp_ratios
dev.off() 