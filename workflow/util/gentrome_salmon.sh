grep "^>" <$1 | cut -d " " -f 1 > $3/decoys.txt 
sed -i.bak -e 's/>//g' $3/decoys.txt 
cat $2 $1 > $3/gentrome.fa.gz