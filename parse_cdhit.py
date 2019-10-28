import argparse
import re
import json
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='Parse CD-HIT output')
parser.add_argument('--clusters', action='store',help='cdhit.clstr file')
parser.add_argument('--fasta', action='store',help='FASTA produced by CD-HIT')
parser.add_argument('--species', action='store',help='Species name for plot title')
args=parser.parse_args()

clusters=open(args.clusters,"r")
out_file=open("clusters_summary.txt","w")

my_clusters={}
my_sequences={}

for line in clusters:
	line=line.rstrip()
	h=re.match('>Cluster\s([0-9]+)$',line)
	l=re.match('[0-9]+',line)
	if h:	
		cluster=h.group(1)
		my_clusters[cluster] = 0
	elif l:
		my_clusters[cluster] += 1
	

for cluster, number in my_clusters.items():
	out_file.write(str(cluster) + "\t" + str(number) + "\n")		


y = my_clusters.values()
x = range(len(y))
#print(x,y)
#plt.bar(x, y, width=10)
plt.plot(x,y,'ro',markersize=1)
plt.ylabel('Cluster size')
plt.title(args.species)
#plt.show()
plt.savefig('clusters.png')
