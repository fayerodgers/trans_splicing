
sed -n "/>Cluster ${2}$/,/^>/p" $1 | grep -o  '>.*;' | while read -r header; do
	grep -A1 $header $3 | grep -v '^-'
done > "cluster_${2}.fasta"

i=0

while read -r line; do 
	echo $line | sed -e "s/>.*/>${i}/" 
	i=$((i + 1))
done > "cluster_${2}_renamed.fasta"  < "cluster_${2}.fasta"

mafft --auto --clustalout cluster_${2}_renamed.fasta > cluster_${2}.aln
