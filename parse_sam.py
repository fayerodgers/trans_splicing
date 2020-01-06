import sys
import json
import argparse
import isoseq
import trans_splicing as ts

parser=argparse.ArgumentParser(description='Find leading coding exons that overlap regions with high levels of soft clipping')
parser.add_argument('--bam', action='store',help='BAM')
parser.add_argument('--gff', action='store',help='GFF')
parser.add_argument('--nreads', action='store',type=int, help='number of reads that have to be clipped at a position to consider it of interest')
parser.add_argument('--nbases', action='store',type=int, help='only consider reads with at least this number of soft clipped bases')
parser.add_argument('--upstream_bases', action='store',type=int, help='consider a position if it is within this number of bases of the annotated start codon')
parser.add_argument('--out_dir', action='store', help='output directory')
parser.add_argument('--read_types', action='store', help='unique, multi or all')

args=parser.parse_args()


gff=open(args.gff,"r")
(genes,transcripts)=isoseq.parse_gff(gff)
clipped_reads = ts.parse_sam(args.bam, args.nbases, args.out_dir, args.read_types, genes, args.upstream_bases)
ts.print_files(clipped_reads,args.out_dir)


