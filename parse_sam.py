import sys
import json
import re
import argparse
sys.path.append('../isoseq_scripts')
import isoseq

parser=argparse.ArgumentParser(description='Find leading coding exons that overlap regions with high levels of soft clipping')
parser.add_argument('--gff', action='store',help='GFF')
parser.add_argument('--nreads', action='store',type=int, help='number of reads that have to be clipped at a position to consider it of interest')
parser.add_argument('--nbases', action='store',type=int, help='only consider reads with at least this number of soft clipped bases')
parser.add_argument('--upstream_bases', action='store',type=int, help='consider a position if it is within this number of bases of the annotated start codon')

args=parser.parse_args()


def add_to_coords(pos,scaffold,coords):
	if scaffold not in coords:
		coords[scaffold]={}
	if pos not in coords[scaffold]:
               coords[scaffold][pos] = 0
        coords[scaffold][pos] += 1
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

def parse_sam(sam, nbases, unique_sam_file, multi_sam_file, unique_clips_file, multi_clips_file):
	unique_mappers={}
	multi_mappers={}
	for line in sam:
		clips=0
        	if re.match('@',line): #headers
			unique_sam_file.write(line)
			multi_sam_file.write(line)
                	continue
        	fields=line.split("\t")
        	cigar=fields[5]
        	start_pos=int(fields[3])
        	scaffold=fields[2]
		read_name=fields[0]
		mapq=int(fields[4])
		seq=fields[9]
		lclip=get_left_clipping(cigar,seq,start_pos,nbases)
		rclip=get_right_clipping(cigar,seq,start_pos,nbases)
		if lclip is not None:
			write_fasta(lclip[1],read_name,mapq, unique_clips_file, multi_clips_file,scaffold,lclip[0])
			if mapq == 255:
				unique_mappers=add_to_coords(lclip[0],scaffold,unique_mappers)
			else:
				multi_mappers=add_to_coords(lclip[0],scaffold,multi_mappers)
			clips = clips + 1
		if rclip is not None:
			write_fasta(rclip[1],read_name,mapq,unique_clips_file, multi_clips_file,scaffold,rclip[0])
                        if mapq == 255:
                                unique_mappers=add_to_coords(rclip[0],scaffold,unique_mappers)
                        else:
                                multi_mappers=add_to_coords(rclip[0],scaffold,multi_mappers)
			clips = clips + 1
		if clips == 0:
			continue
		elif mapq == 255:
			unique_sam_file.write(line)
		else:
			multi_sam_file.write(line)
	return(unique_mappers,multi_mappers)

def write_fasta(sequence,read_name,mapq,unique_clips_file, multi_clips_file,scaffold,pos):
	if mapq == 255:
		unique_clips_file.write(">" + read_name + ";scaffold=" + scaffold + ";coordinate=" + str(pos) + ";\n" + sequence + "\n")
	else:
		multi_clips_file.write(">" + read_name + ";scaffold=" + scaffold + ";coordinate=" + str(pos) + ";\n" + sequence + "\n")
	
def identify_ts_candidates(coords,nreads,genes,transcripts,out_file,upstream_bases):
	for scaffold in coords.keys():
        	for pos in coords[scaffold].keys():
                	if coords[scaffold][pos] < nreads:
                        	continue
                	overlapping=find_overlapping_genes(genes,scaffold,pos,upstream_bases)
                	my_transcripts=get_transcripts(transcripts,overlapping)
                	cds=get_first_cds(transcripts,my_transcripts)
                	exon=get_exon(transcripts,cds)
                	for transcript,exon_coords in exon.items():
                        	if pos >= exon_coords[0] and pos <= exon_coords[1]:			#position falls in the first coding exon
					out_file.write(scaffold + "\t" + str(pos) + "\t" + str(coords[scaffold][pos]) + "\t" + transcript + "\t" + "in_leading_exon\n") #or within a speficied distance of the start codon
				elif my_transcripts[transcript] == '+' and (pos + upstream_bases) <= cds[transcript][0]:
					out_file.write(scaffold + "\t" + str(pos) + "\t" + str(coords[scaffold][pos]) + "\t" + transcript + "\t" + "upstream\n")
				elif my_transcripts[transcript] == '-' and (pos - upstream_bases) >= cds[transcript][1]:
					out_file.write(scaffold + "\t" + str(pos) + "\t" + str(coords[scaffold][pos]) + "\t" + transcript + "\t" + "upstream\n")		
			

##


gff=open(args.gff,"r")
sam=sys.stdin
unique_sam_file=open("./unique_mappers_clipped.sam","w")
multi_sam_file=open("./multi_mappers_clipped.sam","w")
unique_clips_file=open("./unique_clips.fasta","w")
multi_clips_file=open("./multi_clips.fasta","w")
unique_mapper_candidates=open("./unique_mapper_candidates.txt","w")
multi_mapper_candidates=open("./multi_mapper_candidates.txt","w")


(unique_mappers,multi_mappers)=parse_sam(sam, args.nbases, unique_sam_file, multi_sam_file, unique_clips_file, multi_clips_file)
(genes,transcripts)=isoseq.parse_gff(gff)
identify_ts_candidates(unique_mappers,args.nreads,genes,transcripts,unique_mapper_candidates,args.upstream_bases)
identify_ts_candidates(multi_mappers,args.nreads,genes,transcripts,multi_mapper_candidates,args.upstream_bases)



				

#print(json.dumps(coords,indent=4))
#print(json.dumps(genes,indent=4))


