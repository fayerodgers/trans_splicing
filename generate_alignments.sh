
sed -n "/>Cluster ${2}$/,/^>/p" $1 | grep -o  '>.*;' | while read -r header; do
	grep -A1 $header $3 | grep -v '^-'
done > "${4}/cluster_${2}.fasta"

i=0

while read -r line; do 
	echo $line | sed -e "s/>.*/>${i}/" 
	i=$((i + 1))
done > "${4}/cluster_${2}_renamed.fasta"  < "${4}/cluster_${2}.fasta"

#mafft --auto --clustalout ${4}/cluster_${2}_renamed.fasta > ${4}/cluster_${2}.aln
