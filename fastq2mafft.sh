cat $1  | paste - - - - | cut -f 1-2 |  sed -e 's/^@/>/' | tr '\t' '\n' > ${1}.fa

seqtk sample -s100 ${1}.fa 1000 > ${1}_sample.fa

mafft --auto --clustalout ${1}_sample.fa > ${1}_sample.clustal
