cat $1 | head -n 4000 | paste - - - - | cut -f 1-2 |  sed -e 's/^@/>/' | tr '\t' '\n' > 1000_reads_sample.fa

mafft --auto --clustalout 1000_reads_sample.fa > 1000_reads_sample.clustal
