import argparse
import re
import json


parser=argparse.ArgumentParser(description='extract genes associated with specified clusters')
parser.add_argument('--clusters',action='store',nargs='+', help = 'white space separated list of clusters to identify genes from"')
parser.add_argument('--cdhit_clusters',action='store',help='cdhit.clstr file')
parser.add_argument('--metadata', action='store',help='metadata.txt (sample names in first column)')
parser.add_argument('--clusters_text_file',action='store',help='clusters.txt')

args=parser.parse_args()

def parse_clusters(cluster,cluster_file,genes):
	clusters=open(cluster_file,"r")
	x=0
	for line in clusters:
            line=line.rstrip()
            header=re.match('>Cluster\s+(\d+)',line)
            if header:
                x=0
                this_cluster=int(header.group(1))
                if this_cluster == cluster:
                    x=1
                    continue
	    if x==1:
	        genes=parse_line(line,genes)
        clusters.close()
	return(genes)


def parse_line(line,genes):
	gene=re.search('gene=(.[^;]+)',line)
	try:
		this_gene=gene.group(1)
	except ValueError:
		"couldn't identify gene"
	if this_gene not in genes:
			genes[this_gene] = {};
	sample=re.search('sample=([^;]+)',line)
	try:
		this_sample=sample.group(1)
	except ValueError:
			"couldn't identify sample"
	if this_sample not in genes[this_gene]:
		genes[this_gene][this_sample] = 1
	else:
		genes[this_gene][this_sample] += 1
	return(genes)


def lookup_clusters(megacluster,clusters_text_file):
    clusters_file=open(clusters_text_file,"r")
    smaller_clusters=[]
    for line in clusters_file:
        line=line.rstrip()
        temp=line.split("\t")
        if temp[3] == megacluster:
            smaller_clusters.append(int(temp[0]))
    clusters_file.close()
    return(smaller_clusters)



metadata=open(args.metadata,"r")
samples=[]
for line in metadata:
    temp=line.split("\t")
    samples.append(temp[0])
samples=list(set(samples))
samples.sort()

for megacluster in args.clusters:
        smaller_clusters=lookup_clusters(megacluster,args.clusters_text_file)
        genes_in_this_cluster={}
        for cluster in smaller_clusters:
            genes_in_this_cluster=parse_clusters(cluster,args.cdhit_clusters,genes_in_this_cluster)
        FH=open(str(megacluster)+".txt","w")
        FH.write("\t"+"\t".join(samples)+"\n") 
        for gene in genes_in_this_cluster.keys():
            n=[]
            for sample in samples:
                if sample in genes_in_this_cluster[gene]:
                    n.append(str(genes_in_this_cluster[gene][sample]))
                else:
                    n.append("0")
            FH.write(gene+"\t"+"\t".join(n)+"\n")



