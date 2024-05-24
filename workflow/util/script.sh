cd results/lung/fastq/$1
ffq --ftp $1 | grep -Eo '"url": "[^"]*"' | grep -o ftp.sra.*$ | sed 's/.$//' | xargs -n 1 curl -O -C - --silent --show-error