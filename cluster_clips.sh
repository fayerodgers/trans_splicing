#accepts a file with the chromosome and clip position of the sequences to extract, plus species name

while read -r x coord y scaffold ; do 
  grep -A1 "scaffold=${scaffold};coordinate=${coord}" ${2}/clipped_reads.fa
done  < $1 | grep -v '^-' | paste - -  | sort -u | tr "\t" "\n" > ${2}/candidate_clips.fa 

cd-hit-est -i ${2}/candidate_clips.fa -o ${2}/candidate_clips.cdhit -sc -sf -d 100

