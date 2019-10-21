# trans_splicing

Generate genome index:

```
bsub -o index.o -e index.e -R'select[mem>=25000] rusage[mem=25000] span[hosts=1]' -M 25000 -n 32 STAR --runThreadN 32 --runMode genomeGenerate --genomeDir . --genomeFastaFiles v8_final.fa
```

Map:

```
STAR --runThreadN 8 \
--genomeDir . \
--readFilesIn ERR1328128_1.fastq.gz ERR1328128_2.fastq.gz \
--readFilesCommand gunzip -c \
--alignSJoverhangMin 8 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--chimSegmentMin 10 \
--chimOutType WithinBAM \
--outSAMtype BAM SortedByCoordinate
```

Extract clips from uniquely mapping reads:

```
samtools view -h -q 10 hymenolepis_microstoma/Aligned.sortedByCoord.out.bam | python extract_clips.py --n 20
```
Map those clipped reads:

```
bsub -o map.o -e map.e -R'select[mem>=3000] rusage[mem=3000] span[hosts=1]' -M 3000 -n 8
STAR --runThreadN 8 \
--genomeDir . \
--readFilesIn clips.fasta \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 20 \
--alignIntronMax 1 \
--outSAMtype BAM SortedByCoordinate
```


