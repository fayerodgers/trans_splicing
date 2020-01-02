import sys
import json
import re
import pysam
import argparse
import matplotlib.pyplot as plt
import os
#sys.path.append('../isoseq_scripts')
import isoseq

parser=argparse.ArgumentParser(description='Find leading coding exons that overlap regions with high levels of soft clipping')
parser.add_argument('--bam', action='store',help='BAM')
parser.add_argument('--gff', action='store',help='GFF')
parser.add_argument('--nreads', action='store',type=int, help='number of reads that have to be clipped at a position to consider it of interest')
parser.add_argument('--nbases', action='store',type=int, help='only consider reads with at least this number of soft clipped bases')
parser.add_argument('--upstream_bases', action='store',type=int, help='consider a position if it is within this number of bases of the annotated start codon')
parser.add_argument('--out_dir', action='store', help='output directory')

args=parser.parse_args()


def add_to_coords(pos,scaffold,coords,clip_type):
	if scaffold not in coords:
		coords[scaffold]={}
	if pos not in coords[scaffold]:
               coords[scaffold][pos] = {}
	if clip_type not in coords[scaffold][pos]:
		coords[scaffold][pos][clip_type]=0
        coords[scaffold][pos][clip_type] += 1
	return(coords)

def find_overlapping_genes(genes,scaffold,pos,upstream_bases):			#for a given genomic coordinate, return all overlapping genes, and genes within upstream_bases 
	overlapping={}
	for gene in genes.keys():
		if genes[gene]["scaffold"] == scaffold and (genes[gene]["start"] - upstream_bases) <= pos and (genes[gene]["end"] + upstream_bases) >= pos:
			overlapping[gene]=genes[gene]["strand"]
	return(overlapping)

def get_transcripts(transcripts,genes):
	my_transcripts={}
	for gene, strand in genes.items():
		for transcript in transcripts.keys():
			if "parent" not in transcripts[transcript]:	#TODO: check that isoseq.parse_gff handles non-coding transcripts properly
				continue
			if transcripts[transcript]["parent"] == gene:
				my_transcripts[transcript] = strand
	return(my_transcripts) 

def get_first_cds(transcripts,my_transcripts):
	cds_coords={}
	for transcript,strand in my_transcripts.items():
		cds_coords[transcript] = []
		if strand == '+':
			first_cds_start = min(transcripts[transcript]["cds"].keys())
		elif strand == '-':
			first_cds_start = max(transcripts[transcript]["cds"].keys())	
		cds_coords[transcript].append(first_cds_start)
		cds_coords[transcript].append(transcripts[transcript]["cds"][first_cds_start]["end_coord"])
	return(cds_coords)

def get_exon(trancripts,cds_coords):
	exon_coords={}
	for transcript,coords in cds_coords.items():
		exon_coords[transcript]=[]
		for exon in transcripts[transcript]["exons"].keys():
			if exon <= coords[0] and transcripts[transcript]["exons"][exon] >= coords[0]:
				exon_coords[transcript]=[exon,transcripts[transcript]["exons"][exon]]
				break
	return(exon_coords)

def get_left_clipping(cigar,seq,start_pos,nbases):
	left=re.search('^([0-9]+)S',cigar)
	if left:
		bases=int(left.group(1))
		if bases < nbases:
			return None
		pos=start_pos
		clip=seq[0:(bases-1)]	
		clipping=[pos,clip]	 
		return(clipping)
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

def parse_sam(bam, nbases, unique_clips_file, multi_clips_file,out_dir):
	unique_mappers={}
	multi_mappers={}
	sam = pysam.AlignmentFile(bam, "rb")
	out_unique=pysam.AlignmentFile(os.path.join(out_dir,"unique.bam"),"wb",template=sam)
	out_multi=pysam.AlignmentFile(os.path.join(out_dir,"multi.bam"),"wb",template=sam)
	for line in sam.fetch():		
		clips=0
		lclip=get_left_clipping(line.cigarstring,line.query_sequence,line.reference_start,nbases)
		rclip=get_right_clipping(line.cigarstring,line.query_sequence,line.reference_start,nbases)
		if lclip is not None:
			write_fasta(lclip[1],line.query_name,line.mapping_quality, unique_clips_file, multi_clips_file,line.reference_name,lclip[0])
			if line.mapping_quality == 255:
				unique_mappers=add_to_coords(lclip[0],line.reference_name,unique_mappers,'left')
			else:
				multi_mappers=add_to_coords(lclip[0],line.reference_name,multi_mappers,'left')
			clips = clips + 1
		if rclip is not None:
			write_fasta(rclip[1],line.query_name,line.mapping_quality,unique_clips_file, multi_clips_file,line.reference_name,rclip[0])
                	if line.mapping_quality == 255:
	                	unique_mappers=add_to_coords(rclip[0],line.reference_name,unique_mappers,'right')
                	else:
                		multi_mappers=add_to_coords(rclip[0],line.reference_name,multi_mappers,'right')
			clips = clips + 1
		if clips == 0:
			continue
		elif line.mapping_quality == 255:
			out_unique.write(line)
		else:
			out_multi.write(line)
	return(unique_mappers,multi_mappers)

def write_fasta(sequence,read_name,mapq,unique_clips_file, multi_clips_file,scaffold,pos):
	if mapq == 255:
		unique_clips_file.write(">" + read_name + ";scaffold=" + scaffold + ";coordinate=" + str(pos) + ";\n" + sequence + "\n")
	else:
		multi_clips_file.write(">" + read_name + ";scaffold=" + scaffold + ";coordinate=" + str(pos) + ";\n" + sequence + "\n")
	
def determine_splice_type(my_transcripts,clip_type):
	splice_types={}
	for transcript, strand in my_transcripts.items():
		if strand == '+' and clip_type == 'left':
			splice_type = 'acceptor'
		elif strand == '+' and clip_type == 'right':
			splice_type = 'donor'
                elif strand == '-' and clip_type == 'left':
                        splice_type = 'donor'
                elif strand == '-' and clip_type == 'right':
                        splice_type = 'acceptor'
		splice_types[transcript]=splice_type
	return(splice_types)       

def determine_relative_position(exon,pos,my_transcripts,upstream_bases,cds):
	relative_positions={}
	for transcript,exon_coords in exon.items():
		if pos >= exon_coords[0] and pos <= exon_coords[1]: 
			rel_position = 'leading'
		elif my_transcripts[transcript] == '+' and (pos + upstream_bases) <= cds[transcript][0]:
			rel_position = 'leading'
		elif my_transcripts[transcript] == '-' and (pos - upstream_bases) >= cds[transcript][1]:
			rel_position = 'leading'
		else:
			rel_position = 'internal'
		relative_positions[transcript] = rel_position			
	return(relative_positions)		


def identify_ts_candidates(coords,nreads,genes,transcripts,out_file,upstream_bases):
	ts_candidates={}
	for scaffold,d1 in coords.items():
		for pos,d2 in d1.items():
			for clip_type,count in d2.items():
                		if count < nreads:
                        		continue
                		overlapping=find_overlapping_genes(genes,scaffold,pos,upstream_bases)
                		my_transcripts=get_transcripts(transcripts,overlapping)
				splice_types=determine_splice_type(my_transcripts,clip_type)
                		cds=get_first_cds(transcripts,my_transcripts)	#coordinates of the leading CDS
                		exon=get_exon(transcripts,cds)	
				relative_positions=determine_relative_position(exon,pos,my_transcripts,upstream_bases,cds)
				for transcript in my_transcripts:
					splice_type=splice_types[transcript]
					rel_position=relative_positions[transcript]
					splice_rel = str(splice_type + "_" + rel_position)
					out_file.write(scaffold + "\t" + str(pos) + "\t" + str(count) + "\t" + transcript + "\t" + splice_type + "\t"+ rel_position + "\n")
					if transcript not in ts_candidates:
						ts_candidates[transcript] = {}
					if splice_rel not in ts_candidates[transcript]:
						ts_candidates[transcript][splice_rel]=[]
					ts_candidates[transcript][splice_rel].append(pos)

	return(ts_candidates)

def count_sites(ts_candidates):
	counts={}
	splice_types=['donor_internal','donor_leading','acceptor_internal','acceptor_leading']
	for splice_type in splice_types:
		counts[splice_type]={}
		counts[splice_type]['single']=0
		counts[splice_type]['any']=0
	for transcript in ts_candidates.keys():
		for splice_rel in ts_candidates[transcript].keys():
			count=len(ts_candidates[transcript][splice_rel])
			if count == 1:
				counts[splice_rel]['single'] += 1
			counts[splice_rel]['any'] += 1
	return(counts)
							


def plot_counts(counts,transcript_type,out_dir):
	plt.style.use('ggplot')
	y_leading = [counts['acceptor_leading'][transcript_type],counts['donor_leading'][transcript_type]]
	y_internal = [counts['acceptor_internal'][transcript_type],counts['donor_internal'][transcript_type]]
	x=('Acceptor sites', 'Donor sites')
	x_pos = [1,2]
	p1 = plt.bar(x_pos,y_leading)
	p2 = plt.bar(x_pos,y_internal,bottom=y_leading)
#	plt.bar(x_pos,y, color='green')
	plt.ylabel("No. of transcripts")
	plt.title("Putative trans-spliced transcripts")
	plt.xticks(x_pos,x)
	plt.legend((p1[0], p2[0]), ('5prime sites', 'Internal sites'))
#	plt.show()
	plt.savefig(os.path.join(out_dir,(transcript_type + "_sites.png")))
	


##


gff=open(args.gff,"r")
unique_clips_file=open(os.path.join(args.out_dir,"unique_clips.fasta"),"w")
multi_clips_file=open(os.path.join(args.out_dir,"multi_clips.fasta"),"w")
unique_mapper_candidates=open(os.path.join(args.out_dir,"unique_mapper_candidates.txt"),"w")
multi_mapper_candidates=open(os.path.join(args.out_dir,"multi_mapper_candidates.txt"),"w")

(unique_mappers,multi_mappers)=parse_sam(args.bam, args.nbases, unique_clips_file, multi_clips_file,args.out_dir)
(genes,transcripts)=isoseq.parse_gff(gff)
ts_candidates=identify_ts_candidates(unique_mappers,args.nreads,genes,transcripts,unique_mapper_candidates,args.upstream_bases)
print(json.dumps(ts_candidates,indent=4))
counts=count_sites(ts_candidates)
print(json.dumps(counts,indent=4))
plot_counts(counts,'single',args.out_dir)
plot_counts(counts,'any',args.out_dir)

#identify_ts_candidates(multi_mappers,args.nreads,genes,transcripts,multi_mapper_candidates,args.upstream_bases)
#produce_summary_plots()


				

#print(json.dumps(coords,indent=4))
#print(json.dumps(genes,indent=4))


