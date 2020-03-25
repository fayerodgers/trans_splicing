import argparse
import re
import json
import sys
#import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='Parse CD-HIT output')
parser.add_argument('--clusters', action='store',help='cdhit.clstr file', required=True)
parser.add_argument('--fasta', action='store',help='all_clips.fa', required=True)
parser.add_argument('--size', type = int, action = 'store', help = 'print a representative sequence of clusters larger than this')
parser.add_argument('--members', action='store_true')
args=parser.parse_args()


def print_sizes(my_clusters, size, fasta):
    y=[]
    for cluster in my_clusters.keys():
        y.append(my_clusters[cluster]["size"])
        if my_clusters[cluster]["size"] < size:
            sys.stdout.write(str(cluster) + "\t" + str(my_clusters[cluster]["size"]) + "\t" + "NA" + "\n")
            continue
        regex = ">" + my_clusters[cluster]["longest"]
        for line in fasta:
            line=line.rstrip()
            m=re.match(regex,line)
            if m:
                longest_seq = next(fasta)
                longest_seq = longest_seq.rstrip()
                fasta.seek(0)
                break
        sys.stdout.write(str(cluster) + "\t" + str(my_clusters[cluster]["size"]) + "\t" + longest_seq + "\n")

def print_members(my_clusters):
    for megacluster in my_clusters.keys():
        for cluster in my_clusters[megacluster]["IDs"]:
            sys.stdout.write(str(cluster) + "\t" + megacluster + "\n")


clusters=open(args.clusters,"r")
fasta=open(args.fasta,"r") 

my_clusters={}

for line in clusters:
	line=line.rstrip()
	h=re.match('>Cluster\s([0-9]+)$',line)
	l=re.match('[0-9]+',line)
	s=re.search('\*',line)
	if h:	
		cluster=h.group(1)
		my_clusters[cluster] = {}
		my_clusters[cluster]["size"] = 0
		my_clusters[cluster]["longest"] = ""
                my_clusters[cluster]["IDs"] = []
        elif l:
		my_clusters[cluster]["size"] += 1
                b = re.search('>(.+)\.\.\.',line)
                if not b:
                    raise Exception("Couldn't parse sequence names")
                if s:
		    my_clusters[cluster]["longest"] = b.group(1)
                if args.members:
                    my_clusters[cluster]["IDs"].append(b.group(1))
                    
#print(json.dumps(my_clusters,indent=4))

if args.size:
    print_sizes(my_clusters, args.size, fasta)

if args.members:
    print_members(my_clusters)

#def print_members():

#x = range(len(y))
#print(x,y)
#plt.bar(x, y, width=10)
#plt.plot(x,y,'ro',markersize=1)
#plt.ylabel('Cluster size')
#plt.title(args.species)
#plt.show()
#plt.savefig(args.clusters + '_clusters.png')
