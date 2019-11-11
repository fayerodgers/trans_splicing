USAGE="You need to provide the following arguments:
-b BAM file.
-g GFF file.
-r nreads.
-n nbases.
-u n upstream bases.
-s species name.
"

while getopts ":b:g:r:n:u:s:" args
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

		\?)
			printf "$USAGE"
			exit
			;;
	esac
done

echo "Parsing BAM file"	
python ../parse_sam.py --bam $BAM --gff $GFF --nreads $NREADS --nbases $NBASES --upstream_bases $UPSTREAM_BASES

#Make files of trans-splice sites for transcripts with 1) a single donor site 2) a single acceptor site
awk '$5~/donor/{print $4}' unique_mapper_candidates.txt | sort | uniq -c | awk '$1~/^1$/{print $2}' | while read -r id; do    grep $id unique_mapper_candidates.txt | grep 'donor' | cut -f 1,2; done > single_donor_sites

awk '$5~/acceptor/{print $4}' unique_mapper_candidates.txt | sort | uniq -c | awk '$1~/^1$/{print $2}' | while read -r id; do    grep $id unique_mapper_candidates.txt | grep 'acceptor' | cut -f 1,2; done > single_acceptor_sites

echo "Clustering clips"
../cluster_clips.sh single_donor_sites $SPECIES

../cluster_clips.sh single_acceptor_sites $SPECIES

echo "Plotting cluster sizes"
python ../parse_cdhit.py --clusters single_donor_sites.cdhit.clstr  --fasta single_donor_sites.fasta --species $SPECIES

python ../parse_cdhit.py --clusters single_acceptor_sites.cdhit.clstr  --fasta single_acceptor_sites.fasta --species $SPECIES

echo "Generating alignments of top clusters"

TOP_DONOR=$(sort -nr -k2,2 single_donor_sites.cdhit.clstr_summary.txt | head -n 1 | cut -f 1)

TOP_ACCEPTOR=$(sort -nr -k2,2 single_acceptor_sites.cdhit.clstr_summary.txt | head -n 1 | cut -f 1)

../generate_alignments.sh single_donor_sites.cdhit.clstr $TOP_DONOR single_donor_sites.fasta

../generate_alignments.sh single_acceptor_sites.cdhit.clstr $TOP_ACCEPTOR single_acceptor_sites.fasta

