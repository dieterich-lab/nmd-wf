q# Load required packages
suppressPackageStartupMessages({
  library(Biostrings)  # For manipulating DNA sequences
  library(stringr)     # For string operations
  library(ORFik)       # For identifying open reading frames in DNA sequences
  library(here)        # For file path handling
})

# Load utility functions
snakemake@source('utils.R')

# Define a function to find coding sequences in a set of DNA sequences
find_cds <- function(dna_seqs, transcripts){
  # Subset the DNA sequences by the list of transcripts
  dna_seqs <- dna_seqs[transcripts]
  
  # Use the ORFik package to identify open reading frames in the DNA sequences
  orfs <- ORFik:::orfs_as_List(
    fastaSeqs = as.character(dna_seqs, use.names = TRUE), 
    startCodon = "ATG", stopCodon = stopDefinition(1), 
    minimumLength = 0)
  
  # Add transcript IDs to the list of identified ORFs
  orfs$names <- names(dna_seqs)[orfs$index]
  
  # Split the identified ORFs by transcript ID
  orf_ranges <- split(IRanges(orfs$orf[[1]], orfs$orf[[2]]), orfs$names)
  
  # Combine all ORF ranges into a single GRanges object
  all_ranges <- unlist(orf_ranges, use.names = TRUE)
  all_cds <- GRanges(ranges=all_ranges, seqnames = names(all_ranges))
  
  # Rename the CDS IDs with the transcript ID and a rank number
  names(all_cds) <- paste(names(all_cds), make_ranks(names(all_cds)), sep='_')
  
  return(all_cds)
}

# Set the path to the input fasta file
fasta_path <- snakemake@input[[1]]

# Read the DNA sequences from the input fasta file
dna_seqs <- readDNAStringSet(fasta_path)
names(dna_seqs) <- str_split(names(dna_seqs), ' ', simplify = TRUE)[, 1]
# Load the transcript-to-gene annotation data
tx2gene <- readRDS(snakemake@input[[2]])

# Subset the transcript-to-gene annotation data to keep only valid transcripts
tx2gene <- tx2gene[tx2gene$keep, ]
transcripts <- tx2gene$transcript_id

# Identify all coding sequences in the input cDNA sequences
all_cds <- find_cds(
  dna_seqs,
  transcripts)

# Construct a data frame with information about the identified coding sequences
cds_info <- DataFrame(
  transcript = as.factor(seqnames(all_cds)),
  cds_id = names(all_cds),
  ranges = ranges(all_cds)
)

# Extract the DNA sequence of each coding sequence
cds_info$cds_seq <- subseq(dna_seqs[cds_info$transcript], start(cds_info$ranges), end(cds_info$ranges))

# Save the coding sequence data to a file
saveRDS(cds_info, snakemake@output[[1]])
