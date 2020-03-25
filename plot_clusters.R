#!/usr/bin/env Rscript
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir",help="output directory")
parser$add_argument("--clusters",help="clusters.txt file")
parser$add_argument("--title", help="plot title")
args <- parser$parse_args()

args$dir<-'/Users/fr7/git_repos/trans_splicing/h_microstoma'
args$clusters<-'/Users/fr7/git_repos/trans_splicing/h_microstoma/clusters.txt'
args$title<-'Hymenolepis microstoma'

clusters<-read.table(args$clusters)
colnames(clusters)<-c('cluster','n','seq','megacluster')
temp<-clusters[which(!is.na(clusters$megacluster)),]
megaclusters<-clusters[which(is.na(clusters$megacluster)),c('cluster','n')]
for (megacluster in unique(temp$megacluster)){
  x<-temp[which(temp$megacluster == megacluster),'n']
  total<-sum(x)
  this_megacluster<-data.frame(cluster=megacluster,n=total)
  megaclusters<-rbind(megaclusters,this_megacluster)
}


#cluster sizes
plot_cluster_sizes<-function(clusters,dir,title){
  p2<-ggplot(clusters,aes(x=1:nrow(clusters),y=sort(n,decreasing=TRUE)))+
    geom_point(size=1)+
    xlab("Clusters (ordered by size)")+
    ylab("Cluster size (n reads)")+
    theme(text = element_text(size=20), plot.title = element_text(size=20, hjust=0.5), axis.text.x=element_blank())+
    labs(title = title)
  pdf(file.path(dir,paste0(title,".clusters.pdf")))
  print(p2)
  dev.off()
  return(p2)
}

p1<-plot_cluster_sizes(megaclusters,args$dir,args$title)
