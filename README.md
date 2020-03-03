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

#first cluster only exact matches
cd-hit-est -i all_clips.fa -o cdhit_round1 -sc 1 -sf 1 -d 0 -aS 1 -aL 1

#take representative sequences of clusters with > 100 reads and recluster them
python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters cdhit_round1.clstr --fasta all_clips.fa --size 100 > cdhit_round1.cluster_sizes.txt

awk -v OFS="," '$2> 100{print $1,$3}' cdhit_round1.cluster_sizes.txt | sed -e 's/^/>/' | tr "," "\n" > round1_topreps.fa

cd-hit-est -i round1_topreps.fa -o cdhit_round2 -sc 1 -sf 1 -d 0 -aS 1 

#merge clusters of identical subsequences
python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters cdhit_round2.clstr --fasta all_clips.fa --members | sort -k1,1 > round2_megaclusters.txt

join -j 1 -e 'NA' -o 0,1.2,1.3,2.2 -a 1 -t $'\t'  cdhit_round1.cluster_sizes.txt round2_megaclusters.txt | sed -e 's/\([0-9]\)$/\1.mega/g' > clusters.txt

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
python ../../extract_genes.py --clusters 2.mega 3.mega 4.mega --cdhit_clusters cdhit_round1.clstr --metadata ../fastq/trimmed_metadata.txt --clusters_text_file clusters.txt
```

To generate a combined table with three top clusters:
```
#sum all columns for each cluster, eg:
awk -v OFS="\t" '{ for(i=1; i<=NF;i++) j+=$i; print $0, j; j=0 }' 2.mega.txt | sed -e 's/gene://' | sort -k1,1 | sed -e 's/_\([123]\)/_\1_SL2/g' > SL2.txt

#join all three tables
join -a 1 -a 2 -e 0 -1 1 -2 1 -o auto SL1.txt SL2.txt | join -a 1 -a 2 -e 0 -1 1 -2 1 -o auto - SL3.txt > temp

#calculate totals
#check using the correct fields!!
awk -v OFS="\t" '{  j=$17+$33+$49; print $0, j; j=0 }' temp > SLs.txt

awk '$NF>=10{print $0}' SLs.txt > SLs_10reads.txt
```
Plotting
```
#stats on mapping 
for i in $(ls -d */ | sed -e 's/\///'); do $GIT_HOME/trans_splicing/parse_stats.sh $i/out; done

#combine all stats file into one, and add a column for species, then:
Rscript $GIT_HOME/trans_splicing/plot_mapping_rate.R --dir . --stats stats.txt


#cluster sizes
Rscript $GIT_HOME/trans_splicing/plot_clusters.R --dir . --clusters clusters_summary.txt --title ${SPECIES}

#Gene numbers
Rscript $GIT_HOME/trans_splicing/plot_gene_numbers.R --dir . --trans_spliced_genes SLs.txt --library_counts SL1 SL2 SL3 --metadata ../fastq/trimmed_metadata.txt --title $SPECIES
```

