# trans_splicing

Generate BAMs:

```
SPECIES=<species>
mkdir -p $SPECIES/genome  #put genome.fa, wbps_annotation.gff3 and wbps_annotation.gtf here.
mkdir $SPECIES/fastq  #put FASTQs and metadata file here.
mkdir $SPECIES/bams

bsub -o $SPECIES/genome/index.o -e $SPECIES/genome/index.e -R'select[mem>=25000] rusage[mem=25000] span[hosts=1]' -M 25000 -n 32 STAR --runThreadN 32 --runMode genomeGenerate --genomeDir $SPECIES/genome --genomeFastaFiles $SPECIES/genome/genome.fa

cd $SPECIES/bams
perl $GIT_HOME/misc_useful_things/call_star.pl --index_dir  $SPECIES/genome --annotation $SPECIES/genome/wbps.gtf --files ../fastq/metadata.txt

for i in $(ls -d */ | sed -e 's/\///'); do bsub -o ${i}/o -e ${i}/e  -R'select[mem>=10000] rusage[mem=10000] span[hosts=1]' -M 10000 -n 8 $GIT_HOME/misc_useful_things/sort_and_index.sh ${i}/${i}Aligned.out.bam ; done

```

Recommended: use the wrapper to run for each BAM indiviudally, eg:

```
for i in $(ls -d */ | sed -e 's/\///'); do
  mkdir ${i}/out
  bsub -o ${i}/out/splicing.o -e ${i}/out/splicing.e -R'select[mem>=1000] rusage[mem=1000]' -M 1000 \ 
  $GIT_HOME/trans_splicing/trans_splicing.sh \
  -b ${i}/${i}Aligned.out.bam.sorted.bam  \
  -g ../genome/wbps.gff3  \
  -r 1  \         #report positions with at least r clipped reads
  -n 5 \          #minimum length of the clipped portion
  -u 500  \       #consider reads that align within this many bases upstream of a gene (useful if UTRs are not annotated)
  -s ${SPECIES}_${i} \
  -o ${i}/out \
  -t unique
done
```




