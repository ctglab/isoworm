cd [../../results ]$1
ffq --ftp $1 | grep -Eo '"url": "[^"]*"' | grep -o ftp.sra.*$ | sed 's/.$//' | xargs wget -q -c
