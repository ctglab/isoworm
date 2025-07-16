source("/home/runner/work/isoworm/isoworm/workflow/R/config_all.R")
lapply(my_packages, require, character.only = TRUE)
library(data.table)
### Create a list with all txt files containing the samples for each tissues line 
### Extract all samples from the txt list files for each tissues line typology  
setwd(txt_samples)
txt_files_dir <- list()
samples_files <- list()
count = 1
for(j in samples_typologies){
  txt_files_dir    <- c(txt_files_dir , list(assign(paste0(j, "_files"), list.files(full.names = F , recursive =F, pattern= j))))
  samples_files    <- c(samples_files,list(assign(paste0(j, "_samples"),lapply(txt_files_dir[[count]],readLines))))
  count = count + 1 
}

dir_salmon  <- list.dirs(results_dir, recursive = FALSE)
subdir_dir_salmon <- list.dirs(dir_salmon, recursive = FALSE)
salmon_samples <- subdir_dir_salmon[grepl("/salmon", subdir_dir_salmon, ignore.case = TRUE)]
print(salmon_samples)
### Extract tpm count for transcript_1 and transcript_2 for each tissues line typology
## Create a list with all the tpm count for each samples (transcript_1 and transcript_2) inside the results directory  (for each group)
tpm_transcript_2_list   <- list()
tpm_transcript_1_list  <- list()
tpm_transcript_list <- list()
coldata_list  <- list()
for (y in 1:length(samples_files)){
  tpm_transcript_2   <- list()
  tpm_transcript_1  <- list()
  tpm_transcript <- list()
  coldata <- data.frame(files = character(), names = character(), stringsAsFactors = FALSE)
  for (i in 1:length(samples_files[[y]])){
    for(j in samples_files[[y]][[i]]){
      setwd(salmon_samples[y])
      current_sample <- j
      print(j)
      in.list <- list.files(path= '.', recursive = T, full.names = T)
      print(in.list)
      print(paste("Attempting to read file:", in.list[11]))
      transcript_data <- fread(in.list[11])
      transcript_2_tpm  <- subset(transcript_data, Name == name_transcript_1)
      transcript_1_tpm <- subset(transcript_data, Name == name_transcript_2)
      transcript <- subset(transcript_data, Name %in% transcript_ids)
      transcript <- subset(transcript, select = c("Name", "TPM"))
      tpm_transcript <- append(tpm_transcript,transcript)
      tpm_transcript_2  <- append(tpm_transcript_2,transcript_2_tpm$TPM)
      tpm_transcript_1  <- append(tpm_transcript_1,transcript_1_tpm$TPM)
      assign(paste0(samples_typologies[y],"_tpm_transcript"), tpm_transcript)
      assign(paste0(samples_typologies[y],"_tpm_transcript_2"), tpm_transcript_2)
      assign(paste0(samples_typologies[y],"_tpm_transcript_1"), tpm_transcript_1)
      transcript_data <- in.list[11]
      files <- file.path(salmon_samples[y],transcript_data)
      coldata <- rbind(coldata, data.frame(files = files, names = current_sample, stringsAsFactors = FALSE))
    }
  }
  assign(paste0(samples_typologies[y],"_tpm_transcript"), 
  divide_in_sottoliste(get(paste0(samples_typologies[y],"_tpm_transcript"))))
  assign(paste0(samples_typologies[y],"_merged"), 
  Reduce(merge_by_target_id, get(paste0(samples_typologies[y],"_tpm_transcript"))))
  coldata_list[[paste0("coldata", y)]] <- coldata
}

#### ratiooo
transcript_results <- list()
count = 1
for(i in samples_typologies){
  assign(paste0(i,"_transcript_results"), do.call(rbind, Map(data.frame, transcript_transcript_1=get(paste0(i,"_tpm_transcript_1")), transcript_transcript_2=get(paste0(i,"_tpm_transcript_2")))))
  transcript_results[count] <- list(get(paste0(i,"_transcript_results")))
  names(transcript_results[count]) <- i
  colnames(transcript_results[[count]]) <- c("transcript 1", "transcript 2")
  count = count + 1 
}

for(i in 1:length(samples_typologies)){
  transcript_results[[i]]$ratio <- ((transcript_results[[i]]$`transcript 2`+0.01)/(transcript_results[[i]]$`transcript 1`+0.01))
  transcript_results[[i]] <- transcript_results[[i]] %>% 
    pivot_longer(
      cols = `transcript Reference`:ratio,
      names_to = "Isoforms",
      values_to = "value"  
    )
}

db_ratio <- rbindlist(transcript_results)
db_ratio <- subset(db_ratio,Isoforms=="ratio")
db_ratio$group <- NA
db_ratio$IDS   <- NA
tissues_numbers  <- c()
for (i in 1:length(samples_typologies)){
  tissues_numbers  <- c(tissues_numbers,nrow(get(paste0(samples_typologies[i],"_transcript_results"))))
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
bp_ratios <- ggplot(db_ratio, aes(x=group, y=log2(value), fill=group)) + 
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
  scale_fill_manual(values=group.colors_plot) +
  ylab(expression(paste(italic("transcript-204/transcript-220"), " (IsoWorm Salmon)"))) +
  xlab("Tissue Typology") 

# Closing the graphical device
setwd(final_output)
pdf("ratio_salmon.pdf", width = 21.5,height = 11)
bp_ratios
dev.off() 

### Salmon Pie charts
for (y in 1:length(samples_files)){
  df_name <- paste0(samples_typologies[y], "_merged")
  df <- get(df_name)
  df$means <- rowMeans(df[, -1], na.rm = TRUE)
  assign(df_name, df)
}

ens <- fread("transcripts-Summary-Homo_sapiens_Gene_Summary_ENSG00000157764.csv")          
for (y in 1:length(samples_files)){
  n_max <- ncol(get(paste0(samples_typologies[y],"_merged")))
  temp <- get(paste0(samples_typologies[y],"_merged"))
  new <- c("")
  new <- merge(temp,ens, by = "Name")
  assign(paste0(samples_typologies[y],"_plot"), ggplot(new[,c("ID","means")], aes(x="", y=means, fill=ID)) +
           geom_bar(stat="identity", width=1) +
           coord_polar("y", start=0) +
           labs(x = NULL, y = NULL, fill = NULL) +
           theme_classic() +
           geom_label_repel(data =subset(new,ID == "transcript-204" | ID== "transcript-220"),
                            aes(x="", y=means, label = ID),
                            size = 4.5, show.legend = FALSE,
                            fontface = "italic") +
           theme(axis.line = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 plot.title = element_text(hjust = 0.5),
                 legend.text = element_text(face = "italic")) +
           scale_fill_manual(values = group.colors_pie_chart) +
           ggtitle(label_plots[y])  # Use sorted label_ccle vector here
  )
}


#Next, we may draw our plot and table horizontally side-by-side to each other:
plot_final <- NULL
for (sample_name in samples_typologies) {
  plot_name <- paste0(sample_name, "_plot")
  if (exists(plot_name)) {
    if (is.null(plot_final)) {
      plot_final <- get(plot_name)
    } else {
      plot_final <- plot_final + get(plot_name)
    }
  }
}

plot_final <- plot_final + plot_layout(guides = "collect")

setwd(final_output)
pdf("pie_charts.pdf", width = 10 ,height = 10)
plot_final
dev.off() 

## total expression level salmon
## build linked Txome
makeLinkedTxome(indexDir=indexDir,source="GENCODE",
                organism="Homo sapiens",release="43",genome="GRCh38",
                fasta=fasta, gtf=gtf)

se_list <- list()
for (i in 1:length(coldata_list)){
  coldata_list[[i]] <- coldata_list[[i]][!duplicated(coldata_list[[i]]$names), ]
  se <- tximeta(coldata_list[[i]], useHub=FALSE)
  assign(paste0(samples_typologies[i], "_se"), se)
  se_list <- c(se_list,se)
}

# report quantification for genes
for (i in 1:length(se_list)){
  gse <- summarizeToGene(se_list[[i]], countsFromAbundance="lengthScaledTPM")
  gene_level_TPM <- gse@assays@data@listData[["abundance"]] %>% as.data.frame() %>% tibble::rownames_to_column(var="genes")
  gene_level_TPM <- subset(gene_level_TPM, genes == name_gene)
  gene_level_TPM <- as.data.frame(t(as.matrix(gene_level_TPM)))
  gene_level_TPM$groups <- label_plots[[i]]
  gene_level_TPM <- gene_level_TPM[-1, ]
  colnames(gene_level_TPM) <- c("transcript","groups")
  assign(paste0(samples_typologies[i], "_gene_level_TPM"), gene_level_TPM)
}

objects <- ls()
pattern <- "_gene_level_TPM"
selected_objects <- grep(pattern, objects, value = TRUE)
combined_df <- do.call(bind_rows, lapply(selected_objects, get))

plot_total <- ggplot(combined_df,aes(x=groups, y=log2(as.numeric(transcript)+0.01), fill=groups)) + 
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), size=2.0, alpha=0.9, pch=21,fill="black") +
  ylim(0,6) +
  theme_bw () +
  theme(text=element_text(size=28, face = "bold"),
        axis.text=element_text(size=22, face = "bold"),
        axis.title.x =  element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_fill_manual(values =group.colors_plot) +
  ylab("Log2(TPM + 0.01)")

# Closing the graphical device
setwd(final_output)
pdf("total_salmon.pdf", width = 23.5, height = 13)
plot_total
dev.off() 








