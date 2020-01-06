#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

sites<-read.table(args[1])
colnames(sites)<-c('n','candidate_acceptor_sites')
p1<-ggplot(sites,aes(candidate_acceptor_sites,n))+
  geom_col()+
  xlab("Candidate acceptor sites per gene")+
  ylab("Number of genes")
pdf(paste0(args[4],"/hits_per_gene.pdf"),10,5)
print(p1)
dev.off()


distance<-read.table(args[2])
colnames(distance)<-c('n','bp_from_gene_start')
p2<-ggplot(distance,aes(bp_from_gene_start,n))+
  geom_col()+
  xlab("Distance (bp) of TS acceptor from start of gene")+
  ylab("Number of TS acceptor sites")
pdf(paste0(args[4],"/distance_distribution.pdf"),10,5)
print(p2)
dev.off()

clip_length<-read.table(args[3])
colnames(clip_length)<-c('n','nbases_clipped')
p3<-ggplot(clip_length,aes(nbases_clipped,n))+
  geom_col()+
  xlab("Number of clipped bases")+
  ylab("Number of reads")
pdf(paste0(args[4],"/bases_clipped.pdf"),10,5)
print(p3)
dev.off()