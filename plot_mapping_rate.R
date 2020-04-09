#!/usr/bin/env Rscript
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir",help="output directory")
parser$add_argument("--stats",help="stats.txt file")
args <- parser$parse_args()


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#args$dir<-('/Users/fr7/git_repos/trans_splicing/h_microstoma/')
dir<-args$dir
#stats<-read.table('/Users/fr7/git_repos/trans_splicing/h_microstoma/stats.txt',header=T)
stats<-read.table(args$stats,header=T)
stats$unmapped<-stats$total-stats$unique-stats$multi
species<-levels(stats$species)
labels<-c(
  'h_microstoma' = 'Hymenolepis microstoma',
  's_mansoni' = 'Schistosoma mansoni',
  'e_multilocularis' = 'Echinococcus multilocularis'
)

y<-c('unmapped','multi','unique')
x<-'sample'
stats.1<-m$prepare_df_for_plotting(stats,species,x,y)
p1<-ggplot(stats.1,aes(sample,value))+
  geom_col(aes(fill=variable))+
  facet_grid(~species,scales="free", space="free_x", labeller=labeller(species = labels))+
  scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size=20),axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab('Number of reads (pairs)')+
  xlab('Sample')+
  labs(fill='Mapping')
pdf(paste0(dir,'Mapping.pdf'),30,5)
print(p1)
dev.off()