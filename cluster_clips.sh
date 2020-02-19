#accepts a file with the chromosome and clip position of the sequences to extract, plus species name

grep -A1 -v "gene=no_associated_gene" ${1}/clipped_reads.fa  | grep -v '^-' | paste - -  | sort -u | tr "\t" "\n" > ${1}/candidate_clips.fa 

cd-hit-est -i ${1}/candidate_clips.fa -o ${1}/candidate_clips.cdhit -sc -sf -d 100

