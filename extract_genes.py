import argparse
import re
import json


parser=argparse.ArgumentParser(description='extract genes associated with specified clusters')
parser.add_argument('--clusters',action='store',nargs='+', type=int, help = 'white space separated list of clusters to identify genes from')
parser.add_argument('--cdhit_clusters',action='store',help='cdhit.clstr file')
parser.add_argument('--metadata', action='store',help='metadata.txt (sample names in first column)')

args=parser.parse_args()

def parse_clusters(cluster,cluster_file):
	clusters=open(cluster_file,"r")
	genes={}
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


metadata=open(args.metadata,"r")
samples=[]
for line in metadata:
    temp=line.split("\t")
    samples.append(temp[0])
samples=list(set(samples))

for cluster in args.clusters:
        genes_in_this_cluster=parse_clusters(cluster,args.cdhit_clusters)
        print(json.dumps(genes_in_this_cluster,indent=4))
        FH=open(str(cluster)+".txt","w")
        FH.write("\t"+"\t".join(samples)+"\n") 
        for gene in genes_in_this_cluster.keys():
            n=[]
            for sample in samples:
                if sample in genes_in_this_cluster[gene]:
                    n.append(str(genes_in_this_cluster[gene][sample]))
                else:
                    n.append("0")
            FH.write(gene+"\t"+"\t".join(n)+"\n")



