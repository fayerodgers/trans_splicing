USAGE="You need to provide the following arguments:
-b BAM file.
-g GFF file.
-r nreads.
-n nbases.
-u n upstream bases.
-s species name.
-o output directory.
"

while getopts ":b:g:r:n:u:s:o:" args
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

		\?)
			printf "$USAGE"
			exit
			;;
	esac
done

export PYTHONPATH=${GIT_HOME}/isoseq_scripts

echo "Parsing BAM file"	
python ${GIT_HOME}/trans_splicing/parse_sam.py --bam $BAM --gff $GFF --nreads $NREADS --nbases $NBASES --upstream_bases $UPSTREAM_BASES --out_dir $DIR

#Make files of trans-splice sites for transcripts with 1) a single donor site 2) a single acceptor site
awk '$5~/donor/{print $4}' ${DIR}/unique_mapper_candidates.txt | sort | uniq -c | awk '$1~/^1$/{print $2}' | while read -r id; do grep $id ${DIR}/unique_mapper_candidates.txt | grep 'donor' | cut -f 1,2; done > ${DIR}/single_donor_sites

awk '$5~/acceptor/{print $4}' ${DIR}/unique_mapper_candidates.txt | sort | uniq -c | awk '$1~/^1$/{print $2}' | while read -r id; do grep $id ${DIR}/unique_mapper_candidates.txt | grep 'acceptor' | cut -f 1,2; done > ${DIR}/single_acceptor_sites

echo "Clustering clips"
${GIT_HOME}/trans_splicing/cluster_clips.sh ${DIR}/single_donor_sites $DIR

${GIT_HOME}/trans_splicing/cluster_clips.sh ${DIR}/single_acceptor_sites $DIR

echo "Plotting cluster sizes"
python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters ${DIR}/single_donor_sites.cdhit.clstr  --fasta ${DIR}/single_donor_sites.fasta --species $SPECIES

python ${GIT_HOME}/trans_splicing/parse_cdhit.py --clusters ${DIR}/single_acceptor_sites.cdhit.clstr  --fasta ${DIR}/single_acceptor_sites.fasta --species $SPECIES

echo "Generating alignments of top clusters"

TOP_DONOR=$(sort -nr -k2,2 single_donor_sites.cdhit.clstr_summary.txt | head -n 1 | cut -f 1)

if [ -n "$TOP_DONOR" ]; then ${GIT_HOME}/trans_splicing/generate_alignments.sh ${DIR}/single_donor_sites.cdhit.clstr $TOP_DONOR ${DIR}/single_donor_sites.fasta; fi

TOP_ACCEPTOR=$(sort -nr -k2,2 single_acceptor_sites.cdhit.clstr_summary.txt | head -n 1 | cut -f 1)

if [ -n "$TOP_ACCEPTOR" ]; then ${GIT_HOME}/trans_splicing/generate_alignments.sh ${DIR}/single_acceptor_sites.cdhit.clstr $TOP_ACCEPTOR ${DIR}/single_acceptor_sites.fasta; fi

