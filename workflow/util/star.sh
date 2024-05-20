### fastp
fastp -g -y -Y 70 -x --poly_x_min_len 5 -i ../SRR21756204.fastq -o out.fq
### fastqc
fastqc -o ~/rnaseq/results/fastqc/ -t 6 *.fq
## STAR version=2.7.10b
## map single end reads to genome
STAR --runThreadN 2 \
    --readFilesIn out.fq \
    --genomeDir /CTGlab/db/GRCh38.primary_assembly.genome \
	--outFilterType BySJout \
	--alignSJDBoverhangMin 1 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignSJoverhangMin 8 \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 14 \
	--outFilterMismatchNoverLmax 0.2 \
	--outSAMattributes NH HI NM MD \
	--outSAMtype BAM SortedByCoordinate \
