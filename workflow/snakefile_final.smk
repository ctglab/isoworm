from util.varsub import varsub

## config file
configfile: "config/config_final.yml"
varsub(config)

# Check the workflow type
workflow_type = config['workflow_type']

## Determine which genome reference you would like to use
## here we are using GRCh38.primary_assembly.genome.fa + relative gtf/gff3
## depending on the freeze variable inside the config file, the appropriate
## references and data files will be chosen from the config

## If you have modified the R script section that generates the plots, please remove the '#' characters related to the plot in the rule_all part

freeze = config['freeze']

## read list of samples, one per line
with open(config['datadirs']['samples']) as f:
    SAMPLES = f.read().splitlines()

# Add conditions to execute rules based on workflow type (If you have modified the R script section that generates the plots, please remove the '#' characters related to the plot in the rule_all part)
if workflow_type == "polyA_module":
    rule all:
        input:
            starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex",
            small_bam_SE = expand(config['datadirs']['bam'] + "/{file}/" + "{file}_SE_small_Aligned.sortedByCoord.out.bam", file = SAMPLES),
            tab_zip =  expand(config['datadirs']['bam'] + "/{file}/" + "{file}_SE_SJ.out.tab.gz", file = SAMPLES),
            txt = expand(config['datadirs']['bam'] + "/{file}/"+"{file}_SE_reads_count_transcript.txt", file = SAMPLES),
            polyA_transcript_2 = expand(config['datadirs']['r'] + "polyA_filtered_3UTR_transcript_2.csv", file = SAMPLES),
            polyA_transcript_1 = expand(config['datadirs']['r'] + "polyA_filtered_3UTR_transcript_1.csv", file = SAMPLES)
elif workflow_type == "salmon_module":
    rule all:
        input:
            salmon_index = config['refdir'] + "/salmon_index_v43/" + 'ctable.bin',
            salmon_sf =  expand(config['datadirs']['salmon'] + "/{file}/"  + 'quant.sf', file = SAMPLES),
            #boxplot_salmon = expand(config['datadirs']['r'] + "ratio_salmon.pdf", file = SAMPLES),
            #piecharts_salmon = expand(config['datadirs']['r'] + "pie_charts.pdf", file = SAMPLES),
            #total_salmon = expand(config['datadirs']['r'] + "total_salmon.pdf", file = SAMPLES)
elif workflow_type == "custom_module":
    rule all:
        input:
            starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex",
            small_bam = expand(config['datadirs']['bam'] + "/{file}/" + "{file}_small_Aligned.sortedByCoord.out.bam", file = SAMPLES),
            gtf_stringtie = expand(config['datadirs']['stringtie'] + "/{file}/" + "{file}_transcriptome.gtf", file = SAMPLES),
            txt = expand(config['datadirs']['bam'] + "/{file}/"+"{file}_reads_count_transcript.txt", file = SAMPLES),
            tab_zip =  expand(config['datadirs']['bam'] + "/{file}/" + "{file}_SJ.out.tab.gz", file = SAMPLES),
            ballgown = expand(config['datadirs']['ballgown'] + "/{file}/" + "{file}.gtf", file = SAMPLES),
            #boxplot_custom = expand(config['datadirs']['r'] + "ratio_transcript.pdf", file = SAMPLES)
elif workflow_type == "custom_and_salmon_modules":
    rule all:
        input:
            starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex",
            small_bam = expand(config['datadirs']['bam'] + "/{file}/" + "{file}_small_Aligned.sortedByCoord.out.bam", file = SAMPLES),
            gtf_stringtie = expand(config['datadirs']['stringtie'] + "/{file}/" + "{file}_transcriptome.gtf", file = SAMPLES),
            txt = expand(config['datadirs']['bam'] + "/{file}/"+"{file}_reads_count_transcript.txt", file = SAMPLES),
            tab_zip =  expand(config['datadirs']['bam'] + "/{file}/" + "{file}_SJ.out.tab.gz", file = SAMPLES),
            ballgown = expand(config['datadirs']['ballgown'] + "/{file}/" + "{file}.gtf", file = SAMPLES),
            salmon_index = config['refdir'] + "/salmon_index_v43/" + 'ctable.bin',
            salmon_sf =  expand(config['datadirs']['salmon'] + "/{file}/"  + 'quant.sf', file = SAMPLES),
            #boxplot_custom = expand(config['datadirs']['r'] + "ratio_transcript.pdf", file = SAMPLES),
            #boxplot_salmon = expand(config['datadirs']['r'] + "ratio_salmon.pdf", file = SAMPLES),
            #piecharts_salmon = expand(config['datadirs']['r'] + "pie_charts.pdf", file = SAMPLES),
            #total_salmon = expand(config['datadirs']['r'] + "total_salmon.pdf", file = SAMPLES)
elif workflow_type == "singlecells_module":
    rule all:
        input:
            qc_plot =expand(config['datadirs']['singlecell'] +"/{file}/"+ "qc_plots.png", file = SAMPLES),
            filtered_table = expand(config['datadirs']['singlecell'] +"/{file}/"+ "filtered_data.h5ad",file = SAMPLES),
            umap_plot = expand(config['datadirs']['singlecell'] +"/{file}/"+ "umap_plot.png", file = SAMPLES)
            
################### rules to create index files for salmon and star and single cell    ############################################################################################################################################
## create STAR genome index
rule STAR_index:
    input:
        fasta = config['reference']['fasta'][freeze],
        gtf   = config['reference']['gtf'][freeze]
    output:
        directory(config['reference']['stargenomedir'][freeze]),
        starindex = config['reference']['stargenomedir'][freeze] + "/" + "SAindex"
    params:
        star = config['tools']['star'],
    conda:
        config['conda']['STAR']
    threads: 2
    priority: 100
    shell:
        '{params.star} --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'
 

## create salmon gentrome
rule salmon_gentrome:
        input:
            fasta = config['reference']['fasta'][freeze],
            transcriptome_fasta = config['reference']['fasta_salmon'][freeze]
        output:
            temp(config['refdir'] + "/salmon_index_v43/" + 'gentrome.fa.gz')
        params:
            gentrome_sh = config['scripts']['gentrome_sh'],
            outdir=config['refdir'] + "/salmon_index_v43/"
        priority: 100
        shell:
            "sh {params.gentrome_sh} {input.fasta} {input.transcriptome_fasta} {params.outdir}"

## create salmon index
rule salmon_index:
        input:
            gentrome = config['refdir'] + "/salmon_index_v43/" + 'gentrome.fa.gz',
        output:
            salmon_index = config['refdir'] + "/salmon_index_v43/" + 'ctable.bin'
        params:
            outdir=config['refdir'] + "/salmon_index_v43/",
            salmon = config['tools']['salmon'],
            decoys = config['refdir'] + "/salmon_index_v43/" + 'decoys.txt'
        conda:
            config['conda']['salmon']
        priority: 80
        shell:
            """
            salmon index -t {input.gentrome} -i {params.outdir} \
            --decoys {params.decoys} {params.outdir} --gencode -p 2
            """

## create single cell index 
rule singlecell_index:
    input:
        fasta = config['reference']['fasta'][freeze],
        gtf = config['reference']['gtf'][freeze]
    output:
        idx = config['reference']['stargenomedir'][freeze] + "/" + "transcriptome.idx",
        t2g = config['reference']['stargenomedir'][freeze] + "/" + "transcripts_to_genes.txt",
        cdna = config['reference']['stargenomedir'][freeze] + "/" + "cdna.fa"
    params:
        kb = config['tools']['kb']
    conda:
        config['conda']['singlecell']
    threads: 2
    shell:
        "{params.kb} ref -i {output.idx} -g {output.t2g} -f1 {output.cdna} {input.fasta} {input.gtf}"

################### download and clean the read, if bam convert to fastq ##################################################################################################################################################
## Download the samples as fastq (PE-seq)
rule ffq_bash:
    input:
    output:
        f1 = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}_1.fastq.gz"),
        f2 = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}_2.fastq.gz")
    benchmark:
        config['datadirs']['benchmarks'] + "/{file}/" + "{file}_ffq.benchmark"
    params:
        outdir = config['datadirs']['fastq'] + "/{file}/",
        ffq = config['scripts']['bash']
    conda:
        config['conda']['ffq']
    threads: 2
    shell:
        'sh {params.ffq} {wildcards.file}'

## Clean and trim the reads (PE-seq)
rule fastp:
    input:
        f1 = config['datadirs']['fastq'] + "/{file}/" + "{file}_1.fastq.gz",
        f2 = config['datadirs']['fastq'] + "/{file}/" + "{file}_2.fastq.gz"
    output:
        f1_fp = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}_R1.fq"),
        f2_fp = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}_R2.fq")
    benchmark:
        config['datadirs']['benchmarks'] + "/{file}/" + "{file}_fastp.benchmark"
    params:
        fastp = config['tools']['fastp']
    conda:
        config['conda']['fastp']
    threads: 2
    shell:
        '{params.fastp} -i {input.f1} -I {input.f2} -o {output.f1_fp} -O {output.f2_fp}'


## convert bam files to fastq files (only if use bam as input or sliced-bam as input) (PE-seq)
"""
rule bam_to_fq:
    input:
        bam_start = config['datadirs']['bam'] +  "/{file}/"+ "{file}.bam",
    output:
        f1_fp_bam =  temp(config['datadirs']['fastq'] +  "/{file}/"+ "{file}_R1.fq"),
        f2_fp_bam =  temp(config['datadirs']['fastq'] +  "/{file}/"+ "{file}_R2.fq"),
        bam_start_sorted = temp(config['datadirs']['bam'] +  "/{file}/"+ "{file}_sorted.bam")
    benchmark:
        config['datadirs']['benchmarks'] +  "/{file}/"+ "{file}_bamtofq.benchmark"
    params:
        bedtools = config['tools']['bedtools'],
    conda:
        config['conda']['bedtools']
    threads: 2
    shell:
        "{params.samtools} sort -n -o {output.bam_start_sorted} {input.bam_start} | {params.bedtools} bamtofastq -i  {output.bam_start_sorted} -fq {output.f1_fp_bam} -fq2 {output.f2_fp_bam}"
"""

################### align reads with star and creating the sliced bam    ##################################################################################################################################################
## align the reads with star (PE-seq)
rule STAR_align:
    input:
        f1_fp = config['datadirs']['fastq'] +  "/{file}/"+ "{file}_R1.fq",
        f2_fp = config['datadirs']['fastq'] +  "/{file}/"+ "{file}_R2.fq"
    output:
        bam = temp(config['datadirs']['bam'] +  "/{file}/"+ "{file}_Aligned.sortedByCoord.out.bam"),
        tab = config['datadirs']['bam'] +  "/{file}/"+ "{file}_SJ.out.tab"
    benchmark:
        config['datadirs']['benchmarks'] +  "/{file}/"+ "{file}_staralign.benchmark"
    params:
        star = config['tools']['star'],
        genomedir = config['reference']['stargenomedir'][freeze],
        prefix = config['datadirs']['bam']  + "/{file}/"+ "{file}_"
    conda:
            config['conda']['STAR']
    threads: 2
    priority: 50
    shell:
        '{params.star} --runThreadN {threads} '
        '--genomeDir {params.genomedir} '
        '--readFilesIn {input.f1_fp} {input.f2_fp} '
        '--outFileNamePrefix {params.prefix} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outSAMunmapped Within '
        '--outFilterMultimapNmax 50 '
        '--outSAMattributes NH HI NM MD AS nM jM jI XS '
        '--outBAMsortingBinsN 100 '
        '--winAnchorMultimapNmax 100 '


## zip the tabs files to save space
rule zip_tabs:
        input:
            tab = config['datadirs']['bam'] + "/{file}/" + "{file}_SJ.out.tab"
        output:
            tab_zip = config['datadirs']['bam'] + "/{file}/" + "{file}_SJ.out.tab.gz"
        threads: 10
        shell:
            "gzip {input.tab}"

## creating the sliced bam for the transcript gene (PE-seq)
rule SAM_tools:
    input:
        bam = config['datadirs']['bam'] + "/{file}/" + "{file}_Aligned.sortedByCoord.out.bam"
    output:
        bam_idx = temp(config['datadirs']['bam'] + "/{file}/" + "{file}_Aligned.sortedByCoord.out.bam.bai"),
        small_bam = config['datadirs']['bam'] + "/{file}/" + "{file}_small_Aligned.sortedByCoord.out.bam",
        txt = config['datadirs']['bam'] + "/{file}/"+"{file}_reads_count_transcript.txt"
    params:
        samtools = config['tools']['samtools'],
        region = "chr7:140718327-140926929"
    conda:
        config['conda']['Samtools']
    threads: 2
    shell:
        """
        {params.samtools} index {input.bam}
        {params.samtools} view -bh {input.bam} {params.region} > {output.small_bam}
        {params.samtools} view -c -f 1 -F 12 {output.small_bam} >> {output.txt}
        """
################### quantifying the transcipts with stringtie and salmon  ##################################################################################################################################################
## Assemble transcripts for each sample with stringtie and a custom gtf file
rule stringtie_assembly:
        input:
            small_bam = config['datadirs']['bam'] + "/{file}/" + "{file}_small_Aligned.sortedByCoord.out.bam",
            gtf_personal = config['reference']['gtf_personal'][freeze]
        output:
            gtf_stringtie = config['datadirs']['stringtie'] + "/{file}/" + "{file}_transcriptome.gtf"
        benchmark:
            config['datadirs']['benchmarks'] + "/{file}/" + "{file}_stringtie_assembly.benchmark"
        params:
            name = "{file}",
            stringtie = config['tools']['stringtie'],
        conda:
            config['conda']['stringtie']
        threads: 2
        shell:
            "{params.stringtie} -p {threads} -e -c 1.5 -f 0.05 -G {input.gtf_personal} -o {output.gtf_stringtie} -l {params.name} {input.small_bam}"

## Estimate transcript abundances (gtf custom file):
rule stringtie_count:
        input:
         gtf_merge_stringtie = config['datadirs']['stringtie'] + "/{file}/" + "{file}_transcriptome.gtf",
         small_bam = config['datadirs']['bam'] + "/{file}/" + "{file}_small_Aligned.sortedByCoord.out.bam"
        output:
            ballgown = config['datadirs']['ballgown'] + "/{file}/"+ "{file}.gtf"
        benchmark:
            config['datadirs']['benchmarks'] + "/{file}/" + "{file}_stringtie_count.benchmark"
        params:
            stringtie = config['tools']['stringtie'],
        conda:
            config['conda']['stringtie']
        threads: 2
        shell:
            "{params.stringtie} -e -B -p {threads} -c 1.5 -f 0.05 -G {input.gtf_merge_stringtie} -o {output.ballgown} {input.small_bam}"

rule salmon_quant:
        input:
            f1 = config['datadirs']['fastq'] + "/{file}/" + "{file}_R1.fq",
            f2 = config['datadirs']['fastq'] + "/{file}/" + "{file}_R2.fq",
        output:
            salmon_sf = config['datadirs']['salmon'] + "/{file}/"  + "quant.sf"
        params:
            salmon_index = config['refdir'] + "/salmon_index_v43",
            salmon = config['tools']['salmon'],
            libtype = "A",
            zip_ext = "gz",
            extra = "--gcBias --seqBias --reduceGCMemory",
            outdir_all_salmon = directory(config['datadirs']['salmon'] + "/{file}/")
        benchmark:
            config['datadirs']['benchmarks'] + "/{file}/" + "{file}_salmon_quant.benchmark"
        threads: 2
        conda:
            config['conda']['salmon']
        priority: 60
        shell:
            """
            {params.salmon} quant --minAssignedFrags 9 -l {params.libtype} -i {params.salmon_index} -1 {input.f1} -2 {input.f2} -p {threads} -o {params.outdir_all_salmon}
            """


################### polyA module, only Lexogen Quant 3'UTR REV seq data   ################################################################################################################################################
## Download the samples as fastq (SE-seq)
rule ffq_bash_SE:
    input:
    output:
        f1_SE = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}.fastq.gz"),
    benchmark:
        config['datadirs']['benchmarks'] + "/{file}/" + "{file}_ffq.benchmark"
    params:
        outdir = config['datadirs']['fastq'] + "/{file}/",
        ffq = config['scripts']['bash']
    conda:
        config['conda']['ffq']
    threads: 2
    shell:
        'sh {params.ffq} {wildcards.file}'

## trimming and cleaning the reads (SE-seq)
rule fastp_SE:
    input:
        f1_SE = config['datadirs']['fastq'] + "/{file}/" + "{file}.fastq.gz",
    output:
        f1_SE_clean = temp(config['datadirs']['fastq'] + "/{file}/" + "{file}_fastq_1_clean.fq"),
    benchmark:
        config['datadirs']['benchmarks'] + "/{file}/" + "{file}_fastp.benchmark"
    params:
        fastp = config['tools']['fastp']
    conda:
        config['conda']['fastp']
    threads: 2
    shell:
        '{params.fastp} -y -Y 70 -g -x --poly_x_min_len 6 -i {input.f1_SE} -o {output.f1_SE_clean}'

## align using STAR (SE-seq)
rule STAR_align_SE:
    input:
        f1_SE_clean = config['datadirs']['fastq'] + "/{file}/" + "{file}_fastq_1_clean.fq",
    output:
        bam_SE= temp(config['datadirs']['bam'] + "/{file}/" + "{file}_SE_Aligned.sortedByCoord.out.bam"),
        tab = config['datadirs']['bam'] + "/{file}/" + "{file}_SE_SJ.out.tab"
    benchmark:
        config['datadirs']['benchmarks'] + "/{file}/" + "{file}_SE_staralign.benchmark"
    params:
        star = config['tools']['star'],
        genomedir = config['reference']['stargenomedir'][freeze],
        prefix = config['datadirs']['bam'] + "/{file}/" + "{file}_SE_"
    conda:
        config['conda']['STAR']
    threads: 2
    shell:
        "{params.star} --runThreadN {threads} "
        "--readFilesIn {input.f1_SE_clean} "
        "--genomeDir {params.genomedir} "
        "--outFilterType BySJout "
        "--outFileNamePrefix {params.prefix} "
        "--alignSJDBoverhangMin 1 "
        "--alignIntronMin 20 "
        "--alignIntronMax 1000000 "
        "--alignSJoverhangMin 8 "
        "--outFilterMultimapNmax 20 "
        "--outFilterMismatchNmax 14 "
        "--outFilterMismatchNoverLmax 0.2 "
        "--outSAMattributes NH HI NM MD XS "
        "--outSAMtype BAM SortedByCoordinate "


## zip the tabs files to save space
rule zip_tabs_SE:
    input:
        tab = config['datadirs']['bam'] + "/{file}/" + "{file}_SE_SJ.out.tab"
    output:
        tab_zip = config['datadirs']['bam'] + "/{file}/" + "{file}_SE_SJ.out.tab.gz"
    threads: 2
    shell:
        "gzip {input.tab}"

## creating the sliced bam for the transcript gene (SE-seq)
rule SAM_tools_SE:
    input:
        bam_SE= config['datadirs']['bam'] + "/{file}/" + "{file}_SE_Aligned.sortedByCoord.out.bam",
    output:
        bam_idx = temp(config['datadirs']['bam'] + "/{file}/" + "{file}_SE_Aligned.sortedByCoord.out.bam.bai"),
        small_bam_SE = config['datadirs']['bam'] + "/{file}/" + "{file}_SE_small_Aligned.sortedByCoord.out.bam",
        txt = config['datadirs']['bam'] + "/{file}/"+"{file}_SE_reads_count_transcript.txt"
    params:
        samtools = config['tools']['samtools'],
        region = "chr7:140718327-140926929"
    conda:
        config['conda']['Samtools']
    threads: 2
    shell:
        """
        {params.samtools} index {input.bam_SE}
        {params.samtools} view -bh {input.bam_SE} {params.region} > {output.small_bam_SE}
        {params.samtools} view -c -f 1 -F 12 {output.small_bam_SE} >> {output.txt}
        """
################### R code to generate plots and csv files for polyA   ################################################################################################################################################
## create the csv file for the polyA on transcript 3'UTR
rule polyA_r_script:
    input:
    params:
        r_polya = config['scripts']['r_polya']
    output:
        polyA_transcript_2 = config['datadirs']['r'] + "polyA_filtered_3UTR_transcript_2.csv",
        polyA_transcript_1 = config['datadirs']['r'] + "polyA_filtered_3UTR_transcript_1.csv"
    conda:
        config['conda']['r']
    shell:
        "Rscript {params.r_polya}"

## Plot the ratio boxplot
rule custom_r_script:
    input:
    params:
        r_custom = config['scripts']['r_custom']
    output:
        boxplot_custom = config['datadirs']['r'] + "ratio_transcript.pdf"
    conda:
        config['conda']['r']
    shell:
        "Rscript {params.r_custom}"

## Plot the pie charts and the ratio boxplot
rule salmon_r_script:
    input:
        salmon_sf =  expand(config['datadirs']['salmon'] + "/{file}/"  + 'quant.sf', file = SAMPLES)
    params:
        r_salmon = config['scripts']['r_salmon']
    output:
        boxplot_salmon = config['datadirs']['r'] + "ratio_salmon.pdf",
        piecharts_salmon = config['datadirs']['r'] + "pie_charts.pdf",
        total_salmon = config['datadirs']['r'] + "total_salmon.pdf"
    conda:
        config['conda']['r']
    priority: 10
    shell:
        "Rscript {params.r_salmon}"


################### single cell module, select the sequencing technology on the config file   ################################################################################################################################################
## Quality control and filtering
rule pseudoalignment:
    input:
        idx = config['reference']['stargenomedir'][freeze] + "/" + "transcriptome.idx",
        t2g = config['reference']['stargenomedir'][freeze] + "/" + "transcripts_to_genes.txt",
        f1 = config['datadirs']['fastq'] + "/{file}/" + "{file}_1.fastq.gz",
        f2 = config['datadirs']['fastq'] + "/{file}/" + "{file}_2.fastq.gz"
    output:
        h5ad = config['datadirs']['singlecell'] + "/{file}/" + "output.h5ad"
    params:
        kb = config['tools']['kb']
    conda:
        config['conda']['singlecell']
    shell:
        "{params.kb} count -i {input.idx} -g {input.t2g} -o {output.h5ad} -x 10xv3 --h5ad -t 2 {input.f1} {input.f2}"

rule quality_control:
    input:
        h5ad = config['datadirs']['singlecell'] + "/{file}/" + "output.h5ad"
    output:
        filtered=config['datadirs']['singlecell'] + "/{file}/" + "filtered_data.h5ad"
    conda:
        config['conda']['singlecell']
    script:
        "util/script_singlecell/qc_filtering.py {input.h5ad} {ouput.filtered}"

rule pca_plot:
    input:
        filtered=config['datadirs']['singlecell'] + "/{file}/" + "filtered_data.h5ad"
    output:
        pca_plot=config['datadirs']['singlecell'] + "/{file}/" +"pca_plot.png"
    conda:
        config['conda']['singlecell']
    script:
        "util/script_singlecell/pca_plot.py {input.h5ad} {output.pca_plot}"

rule umi_saturation:
    input:
         filtered=config['datadirs']['singlecell'] + "/{file}/" + "filtered_data.h5ad"
    output:
        qc_plots=config['datadirs']['singlecell'] + "/{file}/" +"qc_plots.png"
    conda:
        config['conda']['singlecell']
    script:
        "util/script_singlecell/qc_plots.py {input.h5ad} {output.qc_plots}"

rule clustering_and_umap:
    input:
        filtered=config['datadirs']['singlecell'] + "/{file}/" + "filtered_data.h5ad"
    output:
        umap_plot=config['datadirs']['singlecell'] + "/{file}/" +"umap_plot.png"
    conda:
        config['conda']['singlecell']
    script:
        "util/script_singlecell/umap_clustering.py {input.h5ad} {output.umap_plot}"
