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

Generate multiple alignments of interesting clusters:
```
module load mafft/7.407=1
clusters=(0 1 2)
for i in ${clusters[@]}; do
   bsub -o mafft.${i}.o -e mafft.${i}.e $GIT_HOME/trans_splicing/generate_alignments.sh clips_cdhit.clstr ${i} all_clips.fa
done
```


Extract genes associated with clusters of interest:
```
python $GIT_HOME/trans_splicing/extract_genes.py --clusters 0 1 2 --cdhit_clusters clips_cdhit.clstr --metadata ../fastq/trimmed_metadata.txt
```

To generate a combined table with three top clusters:
```
#sum all columns for each cluster, eg:
awk -v OFS="\t" '{ for(i=1; i<=NF;i++) j+=$i; print $0, j; j=0 }' SL2.txt | sed -e 's/gene://' | sort -k1,1 | sed -e 's/_\([123]\)/_\1_SL2/g' > temp.SL2.txt
#join all three tables
join -a 1 -a 2 -e 0 -1 1 -2 1 -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16 temp.SL1.txt temp.SL2.txt > temp
join -a 1 -a 2 -e 0 -1 1 -2 1 -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17 temp temp.SL3.txt > temp1
#calculate totals
awk -v OFS="\t" '{  j=$17+$33+$49; print $0, j; j=0 }' temp1 | awk '$NF>10{print $0}'> SLs.txt
```
Plotting
```
#stats on mapping 
for i in $(ls -d */ | sed -e 's/\///'); do $GIT_HOME/trans_splicing/parse_stats.sh $i/out; done
#ADD PLOT

#cluster sizes
Rscript $GIT_HOME/trans_splicing/plot_clusters.R --dir . --clusters clusters_summary.txt --title ${SPECIES}

#Gene numbers
Rscript $GIT_HOME/trans_splicing/plot_gene_numbers.R --dir . --trans_spliced_genes SLs.txt --library_counts SL1 SL2 SL3 --metadata ../fastq/trimmed_metadata.txt --title $SPECIES
```

