#!/usr/bin/env Rscript
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir",help="output directory")
parser$add_argument("--trans_spliced_genes",help="Table in which the first column should be gene IDs to include in counts")
parser$add_argument("--library_counts", nargs='+', help="Matrix of gene IDs by library IDs, one file per SL (white space separated)")
parser$add_argument("--metadata",help="metadata.txt, should have columns with 'sample' and 'meta'")
parser$add_argument("--title", help="plot title")
args <- parser$parse_args()

#args$dir<-"/Users/fr7/git_repos/trans_splicing/h_microstoma"
#args$trans_spliced_genes<-"/Users/fr7/git_repos/trans_splicing/h_microstoma/SLs_10reads.txt"
#args$metadata<-"/Users/fr7/git_repos/trans_splicing/h_microstoma/trimmed_metadata.txt"
#args$library_counts<-c('SL1', 'SL2','SL3')
#args$title<-'H microstoma'

ts_genes<-read.table(args$trans_spliced_genes,header=TRUE)
ts_genes<-as.character(ts_genes$gene)


files<-lapply(args$library_counts,function(x) file.path(args$dir,paste0(x,".txt")))
leaders<-lapply(files,function(x) read.table(x, header=TRUE))
names(leaders)<-args$library_counts

sls.df<-data.frame(sample=character(),count=integer(),SL=factor())
genes.df<-data.frame(sample=character(),gene=character())

samples<-colnames(leaders[[1]])
for (sample in samples){
  counts<-c()
  for (sl in names(leaders)){
    genes <- rownames(leaders[[sl]][which(leaders[[sl]][[sample]]>0),])
    genes <- intersect(genes,ts_genes)
    count <- length(genes)
    counts <- c(counts,count)
    if (count > 0){
      these_genes <- data.frame(gene=genes)
      these_genes$sample <- sample
      genes.df <- rbind(genes.df,these_genes)
    }
  }
  this_sample<-data.frame(SL=names(leaders),count=counts)
  this_sample$sample<-sample
  sls.df<-rbind(sls.df,this_sample)
  if (length(leaders)>1){
    any_sl_count<-data.frame(sample=sample,count=length(unique(genes.df[which(genes.df$sample == sample),"gene"])),SL="any")
    sls.df<-rbind(sls.df,any_sl_count)
  }
  
}

metadata<-read.table(args$metadata,header=TRUE)
metadata<-unique(data.frame(sample=metadata$sample,meta=metadata$meta))
sls.df$sample<-as.factor(sls.df$sample)
sls.df<-merge(x=sls.df,y=metadata,all=TRUE)

p2<-ggplot(sls.df,aes(meta,count,fill=SL))+
  geom_boxplot()+
  #    geom_text(aes(label=ifelse(n>5000,as.character(seq),'')),hjust=0,vjust=0,size=2)+
  xlab("Sample")+
  ylab("Genes with >=1 SL reads")+
  theme(text = element_text(size=20), plot.title = element_text(hjust=0.5), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  labs(title = args$plot_title)

pdf(file.path(args$dir,paste0(args$title,".genes.pdf")))
print(p2)
dev.off()
