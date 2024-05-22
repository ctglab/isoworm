ls
ls results/
ls results//fastq
cd results/fastq/$1
ffq --ftp $1 | grep -Eo '"url": "[^"]*"' | grep -o ftp.sra.*$ | sed 's/.$//' | xargs curl -O -C - -o nul
ls