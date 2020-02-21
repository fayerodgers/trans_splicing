#!/usr/bin/env Rscript
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir",help="output directory")
parser$add_argument("--clusters",help="clusters_summary.txt file")
parser$add_argument("--title", help="plot title")
args <- parser$parse_args()


#cluster sizes
plot_cluster_sizes<-function(clusters,dir,title){
  clusters<-read.table(clusters)
  colnames(clusters)<-c('cluster','n','seq')
  p2<-ggplot(clusters,aes(cluster,n,label=seq))+
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

p1<-plot_cluster_sizes(args$clusters,args$dir,args$title)
