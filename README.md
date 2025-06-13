[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
![Testing](https://github.com/mspodda/isoworm/workflows/Testing/badge.svg)
[![GitHub issues](https://img.shields.io/github/issues-raw/mspodda/isoworm)](https://github.com/mspodda/isoworm/issues)
![GitHub open pull requests](https://img.shields.io/github/issues-pr-raw/mspodda/isoworm)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/mspodda/isoworm) 
![GitHub last commit](https://img.shields.io/github/last-commit/mspodda/isoworm)
![GitHub contributors](https://img.shields.io/github/contributors/mspodda/isoworm)
[![GitHub Website](https://img.shields.io/website-up-down-green-red/http/monip.org.svg)](https://ctglab.github.io/)
![GitHub forks](https://img.shields.io/github/forks/mspodda/isoworm)
![GitHub Repo stars](https://img.shields.io/github/stars/mspodda/isoworm)
![GitHub watchers](https://img.shields.io/github/watchers/mspodda/isoworm.svg)

# What is Isoworm
IsoWorm, is a Snakemake pipeline developed to quantify isoforms expression levels in large RNA-seq datasets (paired-end short-reads). 
The pipeline consists of a series of interconnected modules that perform various stages of data analysis. 
It starts with a txt file containing SRA IDs, while the indications about the RNAseq library type, or a BAM file, and the references files (FASTA and GTF) are in the snakemake config file. 
The custom module of IsoWorm could be used to specifically analyse isoforms (in our case study, BRAF), using custom gtf files to quantify isoform-specific genomic regions. 
The quantification is made through Stringtie and all the plots are generated with R language. Conversely, the Salmon module of IsoWorm was used to quantify all the isoforms annotated in Ensembl db (our reference). 
An R script generates pie charts for genen isoform expression. IsoWorm implent also a module for single-end reads tp identifies polyA sites using custom R scripts, starting form Quant Seq 3' REV sequencing data.
An extra module quantify gene expression on single cell data


## Getting Started

### Input

The input files and parameters are specified in `config_final.yml`, and for R plots and script in config file for R:

####  top level directories
- `workflow_type: ""  ` - options: "polyA_module", "salmon_module", "custom_module", "custom_and_salmon_modules", "singlecell_module"
- `sourcedir: ` -  your output directory
- `refdir:    ` -  your gtf fasta and all reference files directory
- `sampledir: ` -  your txt samples files directory
- `envsdir:   ` - your envs files directory
- `workflow:  ` - your workflow (.smk) files directory
- `samples: ` -  your txt file containig the sra samples here!
  
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
- `ratio_salmon.pdf` -  box plots ratio between our two isoforms of interest
- `pie_charts.pdf` - pie charts expressions values of all our isoforms of interest
- `total_salmon.pdf` - total expression levels of our gene of interest
#### custom modules
- `SAindex` - star index
- `{file}_small_Aligned.sortedByCoord.out.bam` - sliced bam of you gene of interest (BRAF in our case study)
- `ratio_BRAF.pdf` -  box plots ratio between our two isoforms of interest

### single cell module (BETA)
- `filtered_data.h5ad` - filtered count from a single single cell sample
- `umap_plot.png`- umap plot with differents cell clusters
- `qc_plots` - quality controll plot
- 
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

Edit `config.yml` to set the input datasets and parameters, edit `config.R` to set the input datasets and parameters for R and edit `script.sh` with the directory where you want to download your fastqs, then issue:
```bash
snakemake -s snakefile_final.smk --use-conda --rerun-incomplete --core 2 -k
```
