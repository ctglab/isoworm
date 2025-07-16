cd results/lung_geo/fastq/$1 # Change only this line, the directory where the FASTQ files will be downloaded
# Download FASTQ files from the SRA using the provided accession number
ffq --ftp $1 | grep -Eo '"url": "[^"]*"' | grep -o ftp.sra.*$ | sed 's/.$//' | xargs -n 1 curl -O -C - --silent --show-error
