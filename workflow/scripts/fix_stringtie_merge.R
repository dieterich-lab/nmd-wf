ref <- rtracklayer::import(snakemake@input[[1]])

t <- rle(ref$transcript_id)
t$values <- make.unique(t$values)
ref$transcript_id <- inverse.rle(t) 

rtracklayer::export(ref, snakemake@output[[1]])
