import sys
import json
import re
import pysam
import matplotlib.pyplot as plt
import os


def parse_sam(bam, nbases, out_dir, alignments_to_analyse, genes, upstream_bases):
        clipped_reads=[]
        sam = pysam.AlignmentFile(bam, "rb")
	if alignments_to_analyse == 'unique':
        	out = pysam.AlignmentFile(os.path.join(out_dir,"unique.bam"),"wb",template=sam)
	elif alignments_to_analyse == 'multi':
		out = pysam.AlignmentFile(os.path.join(out_dir,"multi.bam"),"wb",template=sam)	
        for line in sam.fetch():
		if alignments_to_analyse == 'unique' and  line.mapping_quality != 255:
			continue
		elif alignments_to_analyse == 'multi' and line.mapping_quality == 255:
			continue
		clip = re.search('S|H',line.cigarstring)
		if not clip:
			continue
		lclip = get_left_clipping(line.cigarstring,line.query_sequence,nbases)
		if lclip is not None:
			alignment=initialise_alignment(line.query_name,'left',lclip,line.reference_name,line.reference_start)
			alignment_with_genes=find_overlapping_genes(alignment, genes, upstream_bases)
			clipped_reads = clipped_reads + alignment_with_genes
		rclip = get_right_clipping(line.cigarstring,line.query_sequence,line.reference_start,nbases)
		if rclip is not None:
			alignment=initialise_alignment(line.query_name,'right',rclip[1],line.reference_name,rclip[0])
			alignment_with_genes=find_overlapping_genes(alignment, genes, upstream_bases)
			clipped_reads = clipped_reads + alignment_with_genes
		if lclip is not None or rclip is not None:
			out.write(line)
#	print(json.dumps(clipped_reads,indent=4)) 
	return(clipped_reads)			

def initialise_alignment(read_name,left_or_right,clipped_seq,scaffold,coord):
	alignment={}
	alignment["read_name"] = read_name
	alignment["left_or_right"] = left_or_right
	alignment["nbases"] = len(clipped_seq)
	alignment["scaffold"] = scaffold
	alignment["coord"] = coord
	alignment["clipped_seq"] = clipped_seq
	return(alignment)

def get_left_clipping(cigar,seq,nbases):
        left=re.search('^([0-9]+)S',cigar)
        if left:
                bases=int(left.group(1))
		if bases < nbases:
                        return None
                clip=seq[0:(bases)]
                return(clip)
        else:
                return None

def get_right_clipping(cigar,seq,start_pos,nbases):
        right=re.search('([0-9]+)S$',cigar)
        if right:
                bases=int(right.group(1))
		if bases < nbases:
                        return None
                pattern=re.compile('\d+\D+')
                length=re.compile('\d+')
                operation=re.compile('\D+')
                blocks=pattern.findall(cigar)
                for block in blocks:
                        move=length.search(block)
                        op=operation.search(block)
                        if (op.group() == 'M') or (op.group() == 'D') or (op.group() == '=') or (op.group() == 'X') or (op.group() == 'N'):
                                pos=start_pos + int(move.group())
                clip=seq[-bases:]
                clipping=[pos,clip]
                return(clipping)
        else:
                return None

def check_overlap(this_clip,this_gene,upstream_bases,gene):
	if this_gene["strand"] == "+":
		if (this_gene["start"] - upstream_bases) <= this_clip["coord"] and this_gene["end"] >= this_clip["coord"]:
			this_clip["gene"] = gene
			this_clip["position"] = this_clip["coord"] - this_gene["start"]
			this_clip["orientation"] = determine_orientation("+",this_clip["left_or_right"])
			return(this_clip)
	elif this_gene["strand"] == "-":
		if this_gene["start"] <= this_clip["coord"] and (this_gene["end"] + upstream_bases ) >= this_clip["coord"]:
                        this_clip["gene"] = gene
                        this_clip["position"] = this_gene["end"] - this_clip["coord"]
			this_clip["orientation"] = determine_orientation("-",this_clip["left_or_right"])
			return(this_clip)
	else:
		return(None)	
	
def find_overlapping_genes(alignment,genes,upstream_bases):   
	alignment_list=[]  
	alignment_with_gene = None
	for gene in genes.keys():   
        	if genes[gene]["scaffold"] != alignment["scaffold"]:
			continue
		alignment_with_gene = check_overlap(alignment,genes[gene],upstream_bases,gene)
		if alignment_with_gene is not None:
			alignment_list.append(alignment_with_gene)
	if not alignment_list:
		alignment_list.append(alignment)
	return(alignment_list)			

def determine_orientation(strand,clip_type):
	if strand == '+' and clip_type == 'left':
        	return('acceptor')
        elif strand == '+' and clip_type == 'right':
        	return('donor')
        elif strand == '-' and clip_type == 'left':
        	return('donor')
        elif strand == '-' and clip_type == 'right':
        	return('acceptor')

def print_files(clipped_reads,out_dir,sample_id):
	coords_file = open(os.path.join(out_dir,"coordinates.txt"),"w")
	fasta = open(os.path.join(out_dir,"clipped_reads.fa"),"w")
	keys = ["read_name","left_or_right","nbases","scaffold","coord","gene","position","orientation"]
	for clip in clipped_reads:
		if "gene" in clip:
                    fasta.write(">" + clip["read_name"] + ";scaffold=" + clip["scaffold"] + ";coordinate=" + str(clip["coord"]) + ";gene=" + clip["gene"] + ";orientation=" + clip["orientation"] + ";sample=" + sample_id + ";\n" + clip["clipped_seq"] + "\n")
		to_print=[]
		for k in keys:	
			if k in clip:	
				to_print.append(str(clip[k]))
			else:
				to_print.append("NA")
                to_print.append(sample_id)
		coords_file.write("\t".join(to_print) + "\n")


##
