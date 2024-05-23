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
my_packages <- list(
  "dplyr",
  "ggplot2",
  list(package = "tidyr", source = "cran"),
  list(package = "ggsignif", source = "cran"),
  list(package = "ggpubr", source = "cran"),
  list(package = "patchwork", source = "cran"),
  list(package = "stringr", source = "cran"),
  list(package = "ggrepel", source = "cran"),
  list(package = "GenomicAlignments", source = "bioconductor"),
  list(package = "bedtoolsr", source = "github:nanxstats/bedtoolsr"),
  list(package = "tximport", source = "bioconductor"),
  list(package = "tximeta", source = "bioconductor"),
  list(package = "tibble", source = "cran"))

my_packages_2 <- c("dplyr","ggplot2","tidyr","ggsignif","ggpubr","patchwork","stringr",
                 "ggrepel","GenomicAlignments","bedtoolsr","tximport","tximeta","tibble")

## Sample typologies
samples_typologies <-  c("Lung")

## Labels for the graphic 
label_plots <- c("Lung")

## Colors for the graphics
group.colors_plot <- c()
group.colors_pie_chart <- c()

## directories
txt_samples <- "test_data"
results_dir <- "results"
polyA_bam_dir <- ""
final_output <- "workflow/R"
## total salmon
indexDir <- file.path("results/salmon_index_v43/")
fasta <- file.path("results/salmon_index_v43/")
gtf   <- file.path("results/salmon_index_v43/")

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

install_if_not_present <- function(packages) {
  install_from_cran <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  install_from_bioconductor <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    }
  }
  
  install_from_github <- function(repo) {
    pkg <- tools::file_path_sans_ext(basename(repo))
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github(repo)
    }
  }
  
  for (pkg_info in packages) {
    if (is.character(pkg_info)) {
      install_from_cran(pkg_info)
    } else if (is.list(pkg_info)) {
      pkg <- pkg_info$package
      source <- pkg_info$source
      if (source == "bioconductor") {
        install_from_bioconductor(pkg)
      } else if (source == "github") {
        install_from_github(pkg)
      }
    }
  }
}
