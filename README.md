![Testing](https://github.com/ctglab/isoworm/workflows/Testing/badge.svg)
[![Journal](https://img.shields.io/badge/Journal-DOI-blue?)](https://doi.org/10.1002/1878-0261.70043)
[![GitHub Website](https://img.shields.io/website-up-down-green-red/http/monip.org.svg)](https://ctglab.github.io/)
[![GitHub issues](https://img.shields.io/github/issues-raw/ctglab/isoworm)](https://github.com/ctglab/isoworm/issues)
![GitHub open pull requests](https://img.shields.io/github/issues-pr-raw/ctglab/isoworm)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ctglab/isoworm) 
![GitHub last commit](https://img.shields.io/github/last-commit/ctglab/isoworm)
![GitHub contributors](https://img.shields.io/github/contributors/ctglab/isoworm)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
![GitHub forks](https://img.shields.io/github/forks/ctglab/isoworm)
![GitHub Repo stars](https://img.shields.io/github/stars/ctglab/isoworm)
![GitHub watchers](https://img.shields.io/github/watchers/ctglab/isoworm.svg)

# What is IsoWorm
IsoWorm, is a Snakemake pipeline developed to quantify isoforms expression levels in large RNA-seq datasets (paired-end short-reads). 
The pipeline consists of a series of interconnected modules that perform various stages of data analysis. 
It starts with a txt file containing SRA IDs, while the indications about the RNAseq library type, or a BAM file, and the references files (FASTA and GTF) are in the snakemake config file. 
The custom module of IsoWorm could be used to specifically analyse isoforms (in our case study, BRAF), using custom gtf files to quantify isoform-specific genomic regions. 
The quantification is made through Stringtie and all the plots are generated with R language. Conversely, the Salmon module of IsoWorm was used to quantify all the isoforms annotated in Ensembl db (our reference). 
An R script generates pie charts for isoform expression. IsoWorm implent also a module for single-end reads tp identifies polyA sites using custom R scripts, starting form Quant Seq 3' REV sequencing data.

## Cite 
> *Podda MS, Tatoni D, Mattei G, Magi A, D'Aurizio R, Poliseno L. Landscape of BRAF transcript variants in human cancer. Mol Oncol. 2025 May 25. doi: 10.1002/1878-0261.70043. Epub ahead of print. PMID: 40415485.*   
> Published in *Molecular Oncology*, 2024


## Getting Started

### Input Files and Configuration
To run this Snakemake workflow, you need a .txt file that lists your sample names. This file is used to assign the '{file}' wildcard within the workflow, which links each rule's output to the correct sample.

For example, if your .txt file contains the following sample names (column):

Sample1  
Sample2  
Sample3  

Then your FASTQ files must follow this naming convention:

Sample1_1.fastq.gz or Sample_1_R1.fastq
Sample1_2.fastq.gz or Sample_1_R2.fastq
Sample2_1.fastq.gz or Sample_2_R1.fastq
Sample2_2.fastq.gz or Sample_2_R2.fastq
Sample3_1.fastq.gz or Sample_3_R1.fastq 
Sample3_2.fastq.gz or Sample_3_R2.fastq

Change the directory where you need to download te fastq from SRA in the first line of 'download_fastq_script.sh'.
You can start directly from the FASTQ files without downloading them or cleaning them, as long as their filenames match the entries in the .txt file as in the example.

#### Configuration
Edit the config.yaml file to set the appropriate parameters for your isoform analysis.

#### Warning R plots
If you want to generate plots in R, make sure to update the appropriate settings in the R configuration file. Please note that in the R scripts used to generate the plots, the custom module is set by default to display only two transcripts in the plot, as well as in the boxplot of the Salmon module (but not in the pie chart). If you wish to plot more than two isoforms, you will also need to modify the corresponding R script. This limitation applies only to the plotting step and does not affect the quantification part of the pipeline, which is managed by Snakemake.

### Config_final.yml mandatory parameters
####  top level directories
- `workflow_type: ""  ` - options: "polyA_module", "salmon_module", "custom_module", "custom_and_salmon_modules", "singlecell_module"
- `sourcedir: ` -  your output directory
- `refdir:    ` -  your gtf fasta and all reference files directory
- `sampledir: ` -  your txt samples files directory
- `envsdir:   ` - your envs files directory
- `workflow:  ` - your workflow (.smk) files directory
- `samples: ` -  your txt file containig the sra samples or samples name here!
  
#### reference files, genome indices and data
- `stargenomedir, GRCh38.primary_assembly.genome:` - directory for STAR genome and Single cell index files
- `fasta: GRCh38.primary_assembly.genome:` -  genome fasta reference file for STAR and Single cell
- `fasta_salmon: GRCh38.primary_assembly.genome:` -  transcript fasta reference for salmon
- `gtf: GRCh38.primary_assembly.genome:` - gtf file for all transcripts
- `gtf_personal: GRCh38.primary_assembly.genome:` - gtf file customize for your transcript of interest 


### Output
#### polyA modules
- `SAindex` - star index
- `{sample_name}_SE_small_Aligned.sortedByCoord.out.bam` - sliced bam of you gene of interest (BRAF in our case study), single end
- `polyA_filtered_3UTR204.csv` - peaks for poly A in BRAF-204 UTRs
- `polyA_filtered_3UTR220.csv` - peaks for poly A in BRAF-220 UTRs
#### salmon modules
- `salmon_index` - salmon index
- `quant.sf` - all transcript quantified by salmon
- `ratio_salmon.pdf` -  box plots ratio between of two isoforms of interest
- `pie_charts.pdf` - pie charts expressions values of all our isoforms of interest
- `total_salmon.pdf` - total expression levels of our gene of interest
#### custom modules
- `SAindex` - star index
- `{file}_small_Aligned.sortedByCoord.out.bam` - sliced bam of you gene of interest (BRAF in our case study)
- `ratio_BRAF.pdf` -  box plots ratio between our two isoforms of interest

#### single cell module (BETA)
- `filtered_data.h5ad` - filtered count from a single single cell sample
- `umap_plot.png`- umap plot with differents cell clusters
- `qc_plots` - quality controll plot

### Dependencies
- [miniconda](https://conda.io/miniconda.html) - install it according to the [instructions](https://conda.io/docs/user-guide/install/index.html).
- [snakemake](https://anaconda.org/bioconda/snakemake) install using `conda`.
- The rest of the dependencies are automatically installed using the `conda` feature of `snakemake`.

### Installation

Clone the repository:

```bash
git clone https://github.com/ctglab/isoworm
```

### Usage

Edit `config_final.yml` to set the input datasets and parameters, edit `config.R` to set the input datasets and parameters for R if you want the plots and edit `download_fastq_script.sh` with the directory where you want to download your fastqs, then issue:
```bash
snakemake -s snakefile_final.smk --use-conda --rerun-incomplete --core 2 -k
```
