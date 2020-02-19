import argparse
import re
import json
import sys
#import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='Parse CD-HIT output')
parser.add_argument('--clusters', action='store',help='cdhit.clstr file')
parser.add_argument('--fasta', action='store',help='FASTA produced by CD-HIT')
args=parser.parse_args()

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
	elif l:
		my_clusters[cluster]["size"] += 1
		if s:
			b = re.search('(>.+)\.\.\.',line)
			if b:
				my_clusters[cluster]["longest"] = b.group(1)
			
#print(json.dumps(my_clusters,indent=4))


y=[]
for cluster in my_clusters.keys():
	y.append(my_clusters[cluster]["size"])
	if my_clusters[cluster]["size"] < 10:
		sys.stdout.write(str(cluster) + "\t" + str(my_clusters[cluster]["size"]) + "\t" + "NA" + "\n")
		continue
	for line in fasta:
		line=line.rstrip()
		regex = my_clusters[cluster]["longest"]
		m=re.match(regex,line)
		if m:
			longest_seq = next(fasta)
			fasta.seek(0)
			break
	sys.stdout.write(str(cluster) + "\t" + str(my_clusters[cluster]["size"]) + "\t" + longest_seq + "\n")


#x = range(len(y))
#print(x,y)
#plt.bar(x, y, width=10)
#plt.plot(x,y,'ro',markersize=1)
#plt.ylabel('Cluster size')
#plt.title(args.species)
#plt.show()
#plt.savefig(args.clusters + '_clusters.png')
