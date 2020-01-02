#accepts a file with the chromosome and clip position of the sequences to extract, plus species name

while read -r scaffold coord ; do 
  grep -A1 "scaffold=${scaffold};coordinate=${coord}" ${2}/unique_clips.fasta
done  < $1 | grep -v '^-' | paste - -  | sort -u | tr "\t" "\n" > ${1}.fasta 

cd-hit-est -i ${1}.fasta -o ${1}.cdhit -sc -sf -d 100

