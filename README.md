# trans_splicing

Generate BAMs:

```
SPECIES=<species>
mkdir -p $SPECIES/genome  #put genome.fa, wbps_annotation.gff3 and wbps_annotation.gtf here.
mkdir $SPECIES/fastq  #put FASTQs and metadata file here.
mkdir $SPECIES/bams #BAMs will go here

#Trim adaptor sequences from FASTQs
module load trimmomatic/0.39--1 
perl ${GIT_HOME}/misc_useful_things/call_trimmomatic.pl --output_directory ${SPECIES}/fastq --files ${SPECIES}/fastq/metadata.txt --fasta <path/to/trimmomatic/adapter/file.fa>

#Generate genome index for mapping
bsub -o $SPECIES/genome/index.o -e $SPECIES/genome/index.e -R'select[mem>=25000] rusage[mem=25000] span[hosts=1]' -M 25000 -n 32 STAR --runThreadN 32 --runMode genomeGenerate --genomeDir $SPECIES/genome --genomeFastaFiles $SPECIES/genome/genome.fa --genomeSAindexNbases 12

#Map
cd $SPECIES/bams
perl $GIT_HOME/misc_useful_things/call_star.pl --index_dir  ../genome --annotation ../genome/wbps.gtf --files ../fastq/trimmed_metadata.txt

#Sort and index
for i in $(ls -d */ | sed -e 's/\///'); do bsub -o ${i}/sort.o -e ${i}/sort.e  -R'select[mem>=10000] rusage[mem=10000] span[hosts=1]' -M 10000 -n 8 $GIT_HOME/misc_useful_things/sort_and_index.sh ${i}/${i}Aligned.out.bam ; done

```

Parse BAMs to extract soft clips (clipped_reads.fa) and clipping coordinates (coordinates.txt):

```
for i in $(ls -d */ | sed -e 's/\///'); do
   mkdir ${i}/out
   bsub -o ${i}/out/splicing.o -e ${i}/out/splicing.e -R'select[mem>=1000] rusage[mem=1000]' -M 1000 \
   python ${GIT_HOME}/trans_splicing/parse_sam.py --bam ${i}/${i}Aligned.out.bam.sorted.bam  --gff ../genome/wbps.gff3 --nreads 1 --nbases 5 --upstream_bases 500 --out_dir ${i}/out --read_types unique --sample_id $i
done
```
Combine clips from all libraries and cluster:

```

for i in $(ls -d */ | sed -e 's/\///'); do
   grep -A1 "orientation=acceptor" ${i}/out/clipped_reads.fa  | grep -v '^-' | paste - -  | sort -u | tr "\t" "\n" >> all_clips.fa 
done

cd-hit-est -i all_clips.fa -o clips_cdhit -sc 1 -sf 1 -d 0
```

Summarise clusters:
```
bsub -o summarise_clusters.o -e summarise_clusters.e -R'select[mem>=1000] rusage[mem=1000]' -M 1000 \
"python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters clips_cdhit.clstr --fasta clips_cdhit | sort -nr -k2,2 > clusters_summary.txt"
```

Extract genes:
```
python $GIT_HOME/trans_splicing/extract_genes.py --clusters 0 1 2 --cdhit_clusters clips_cdhit.clstr --metadata ../fastq/metadata.txt
```




