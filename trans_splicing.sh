USAGE="You need to provide the following arguments:
-b BAM file.
-g GFF file.
-r nreads.
-n nbases.
-u n upstream bases.
-s species name.
-o output directory.
-t alignments to analyse (unique/multi/all).
"

while getopts ":b:g:r:n:u:s:o:t:" args
do
	case $args in

		b)
			BAM=$OPTARG
			;;
		
		g)
			GFF=$OPTARG
			;;
		
		r)
			NREADS=$OPTARG
			;;

		n)
			NBASES=$OPTARG
			;;

		u)
			UPSTREAM_BASES=$OPTARG
			;;

		s)
			SPECIES=$OPTARG
			;;

		o)
			DIR=$OPTARG
			;;		
			
		t)
			TYPE=$OPTARG
			;;

		\?)
			printf "$USAGE"
			exit
			;;
	esac
done


echo "Parsing BAM file"	
python ${GIT_HOME}/trans_splicing/parse_sam.py --bam $BAM --gff $GFF --nreads $NREADS --nbases $NBASES --upstream_bases $UPSTREAM_BASES --out_dir $DIR --read_types $TYPE

echo "Generating summary statistics"
${GIT_HOME}/trans_splicing/parse_stats.sh $DIR

echo "Drawing plots"
source activate r_env
Rscript ${GIT_HOME}/trans_splicing/library_plots.R ${DIR}/hits_per_transcript.txt ${DIR}/distance_distribution.txt ${DIR}/bases_clipped_distribution.txt $DIR
conda deactivate 

echo "Clustering clips of candidate transcripts"
${GIT_HOME}/trans_splicing/cluster_clips.sh ${DIR}/transcripts.txt $DIR

echo "Plotting cluster sizes"
python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters ${DIR}/candidate_clips.cdhit.clstr  --fasta ${DIR}/candidate_clips.fa --species $SPECIES

echo "Generating alignments of top clusters"
TOP_ACCEPTOR=$(sort -nr -k2,2 candidate_clips.cdhit.clstr_summary.txt | head -n 1 | cut -f 1)
if [ -n "$TOP_ACCEPTOR" ]; then ${GIT_HOME}/trans_splicing/generate_alignments.sh ${DIR}/candidate_clips.cdhit.clstr $TOP_ACCEPTOR ${DIR}/candidate_clips.fa; fi

