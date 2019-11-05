# trans_splicing

Generate genome index:

```
bsub -o index.o -e index.e -R'select[mem>=25000] rusage[mem=25000] span[hosts=1]' -M 25000 -n 32 STAR --runThreadN 32 --runMode genomeGenerate --genomeDir . --genomeFastaFiles v8_final.fa
```

Map reads:

```
STAR --runThreadN 8 \
--genomeDir . \
--readFilesIn fastq_1.gz fastq_2.gz \
--readFilesCommand gunzip -c \
--alignSJoverhangMin 8 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--chimSegmentMin 10 \
--chimOutType WithinBAM \
--outSAMtype BAM SortedByCoordinate
```

Parse BAM file:

```
bsub -o test.o -e test.e -R'select[mem>=500] rusage[mem=500] span[hosts=1]' -M 500 \
'samtools view -h ./Aligned.sortedByCoord.out.bam | python ../parse_sam.py \
--gff GFF \
--nreads 10 \
--nbases 15 \
--upstream_bases 500'
```

Produces:

* ```unique_mappers_clipped.sam``` - SAM file of uniquely mapping reads with >```nbases``` soft clipped.
* ```multi_mappers_clipped.sam``` - SAM file of non-uniquely mapping reads with >```nbases``` soft clipped. 
* ```unique_clips.fasta``` - FASTA file of the clipped sequences from ```unique_mappers_clipped.sam```.
* ```multi_clips.fasta``` - FASTA file of the clipped sequences from ```multi_mappers_clipped.sam```.
* ```unique_mapper_candidates.txt``` - A file of transcripts with candidate splice leader acceptor sites. To appear in the file, a transcript must have a position in the exon containing the start codon, or within ```upstream_bases``` upstream of it, with at least ```nreads``` uniquely mapping reads being soft clipped at this position by at least ```nbases```. 
* ```multi_mapper_candidates.txt``` - As above, but for non-uniquely mapping reads.

Extract the clipped bases of reads that hit transcripts with a single candidate acceptor site:

```
cut -f 4 unique_mapper_candidates.txt | sort | uniq -c | awk '$1~/^1$/{print $2}' | while read -r id; do 
  grep $id unique_mapper_candidates.txt
done | cut -f 1,2 | while read -r scaffold coord ; do 
  grep -A1 "scaffold=${scaffold};coordinate=${coord}" unique_clips.fasta
done | grep -v '^-' | paste - -  | sort -u | tr "\t" "\n" > top.fasta
```

Cluster these clips with cd-hit:

```
cd-hit-est -i top.fasta -o top.cdhit -sc -sf -d 100
```

Summarise the clusters found by cd-hit (produces a table and plot of cluster size distribution):

```
python ../parse_cdhit.py --clusters top.cdhit.clstr --fasta top.cdhit --species species
```

Send interesting looking clusters to MAFFT to produce multiple alignments:

```
./generate_alignments.sh schistosoma_mansoni/top.cdhit.clstr cluster schistosoma_mansoni/top.fasta 
```


