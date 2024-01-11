#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

f<-args[1];
print(f);
gr_obj =  rtracklayer::import.bed(f)
bed=read.delim(f,as.is=T,header=F)
ex_by_tx <- readRDS("ex_by_tx.RDS")

for(k in 1:length(gr_obj)){
    gnm_tx <- mapToTranscripts(gr_obj[k],ex_by_tx)
    gnm_df <- as.data.frame(gnm_tx)
    gnm_df <- data.frame(bed[rep(k,nrow(gnm_df)),],TxId=gnm_df$seqnames,cDNAStart=gnm_df$start,cDNAEnd=gnm_df$end) 
    write.table(gnm_df,file=gsub(".bed","_cDNA.bed",f),quote=F,sep="\t",append=T,  row.names=F, col.names=F)
}

