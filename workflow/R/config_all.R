## file:config_all.R
## values
length_x1   <- 140725145-140719327
length_ref  <- 140732564-140730665
length_BRAF <- 9807	
# 3UTR ref
end_3utr_ref   <- 140730665 - 100
start_3utr_ref <- 140734900 + 100
# 3UTR X1
end_3utr_X1   <- 140719327  - 100
start_3utr_X1 <- 140726516  + 100
## BRAF
end_BRAF   <- 140719327 - 100
start_BRAF <- 140924929 + 100
BRAF_ids <- c("ENST00000288602.11","ENST00000469930.2","ENST00000479537.6","ENST00000496384.7",
              "ENST00000497784.2","ENST00000642228.1","ENST00000642272.1","ENST00000642808.1",
              "ENST00000642875.1","ENST00000643356.1","ENST00000643790.1","ENST00000644120.1",
              "ENST00000644650.1","ENST00000644905.1","ENST00000644969.2","ENST00000645443.1",
              "ENST00000646334.1","ENST00000646427.1","ENST00000646730.1","ENST00000646891.2",
              "ENST00000647434.1")

## libraries
my_packages <- c("dplyr","ggplot2","tidyr","ggsignif","ggpubr","patchwork","stringr",
                 "ggrepel","bedtoolsr","tibble","GenomicAlignments","tximeta",
                 "tximport")

## Sample typologies
samples_typologies <-  c("lung_geo")

## Labels for the graphic 
label_plots <- c("lung_geo")

## Colors for the graphics
group.colors_plot <- c()
group.colors_pie_chart <- c( "#708090", "#E6E6E6", "#D9D9D9", "#008bbf", 
                             "#BFBFBF", "#B3B3B3", "#A6A6A6", "#999999",
                                         "#8C8C8C", "#808080", "#737373", "#666666",
                                         "#595959", "#4D4D4D", "#FFD700", "#333333",
                                         "#262626", "#1A1A1A", "#0D0D0D", "#f28e6e", 
                                         "#cfcfc4")

## directories
txt_samples   <- "/home/runner/work/isoworm/isoworm/test_data"
results_dir   <- "/home/runner/work/isoworm/isoworm/results"
polyA_bam_dir <- "/home/runner/work/isoworm/isoworm/results/polyA/bam"
final_output  <- "/home/runner/work/isoworm/isoworm/workflow/R"
## total salmon
indexDir <- file.path("/home/runner/work/isoworm/isoworm/test_data/salmon_index_v43")
fasta    <- file.path("/home/runner/work/isoworm/isoworm/test_data/chr7_transcripts.fa")
gtf      <- file.path("/home/runner/work/isoworm/isoworm/test_data/chr7.gtf")

## functions
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Funzione ricorsiva per dividere la lista in sottoliste da 2 elementi
divide_in_sottoliste <- function(lst) {
  if (length(lst) <= 1) {
    return(list(lst))  # Se rimane solo un elemento, restituisci una lista con quell'elemento
  } else {
    current_pair <- lst[1:2]  # Prendi i primi due elementi
    remaining_list <- lst[-(1:2)]  # Rimuovi i primi due elementi dalla lista
    
    if (length(remaining_list) > 0) {
      # Chiamata ricorsiva solo se la parte rimanente della lista non è vuota
      recursive_result <- divide_in_sottoliste(remaining_list)
      
      # Combina la coppia corrente con il risultato ricorsivo
      final_result <- c(list(current_pair), recursive_result)
      
      return(final_result)
    } else {
      return(list(current_pair))
    }
  }
}

## merge target
merge_by_target_id <- function(df1, df2) {
  return(merge(df1, df2, by="Name", all = TRUE))
}