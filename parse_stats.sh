
#parse star stats and coordinates.txt

echo sample$'\t'total_reads$'\t'unique$'\t'multi$'\t'clipped > ${1}/stats.txt

total=$(grep 'Number of input reads' ${1}/../*Log.final.out | grep -o [0-9]\*)
unique=$(grep 'Uniquely mapped reads number' ${1}/../*Log.final.out | grep -o [0-9]\*)
multi=$(grep 'Number of reads mapped to multiple loci' ${1}/../*Log.final.out | grep -o [0-9]\*)
clipped=$(cut -f 1 ${1}/coordinates.txt | sort | uniq | wc -l)
#in_gene=$(awk '$6 !~/NA/{print $1}' ${1}/coordinates.txt  | sort | uniq | wc -l)
#in_gene_and_acceptor=$(awk '$8 ~ /acceptor/{print $1}' ${1}/coordinates.txt  | sort | uniq | wc -l)
#awk  '$6 !~/NA/ && $8 ~ /acceptor/{print $1,$5,$6,$4}' ${1}/coordinates.txt | sort | uniq | awk '{print $2,$3,$4}' | sort | uniq -c | sort -nr -k1,1 | awk '$1>=5{print $0}' > ${1}/transcripts.txt
#awk '$1>=5{print $3}' ${1}/transcripts.txt | sort | uniq -c | awk '{print $1}' | sort | uniq -c | sort -n -k2,2 > ${1}/hits_per_transcript.txt	
#candidate_transcripts=$(wc -l ${1}/transcripts.txt)
#awk '{print $3}' ${1}/coordinates.txt | sort | uniq -c | sort -n -k2,2 > ${1}/bases_clipped_distribution.txt
#awk '$7 !~ /NA/{print $7}' ${1}/coordinates.txt | sort | uniq -c | sort -n -k2,2 > ${1}/distance_distribution.txt	
echo ${1}$'\t'${total}$'\t'${unique}$'\t'${multi}$'\t'${clipped}$'\t'${in_gene}$'\t'${in_gene_and_acceptor}$'\t'${candidate_transcripts} >> ${1}/stats.txt	




	
	
