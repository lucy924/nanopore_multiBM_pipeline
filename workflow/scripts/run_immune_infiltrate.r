# run_immune_infiltrate.r
# NOTE - I don't think this is the way. Have abandoned for now.
# This code looks like it's for flow sorted (minfi)

library(minfi)
library(tibble)
library(dplyr)
library(GenomicRanges)

snakemake@source("shared_functions.r")

# ------------------------------------------------ #
sample_name <- snakemake@input[["sample_name"]]

dss_table <- read.table(file = snakemake@input[["epic_probe_results"]], sep = ",", header = TRUE)

# ------------------------------------------------ #
# Make GRanges object

# import chrom_sizes
chrom_sizes <- read.table(file = snakemake@input[["chrom_sizes"]], sep = "\t", header = FALSE)
named_integer_vector = chrom_sizes$V2
names(named_integer_vector) = chrom_sizes$V1
seq_info_hg38 = Seqinfo(genome = "hg38")

# https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRanges-class.html
GR_from_dss = makeGRangesFromDataFrame(dss_table,
    keep.extra.columns = FALSE,
    ignore.strand = TRUE,
    seqinfo = seq_info_hg38,
    seqnames.field = "chr",
    start.field = "pos",
    end.field = "endpos",
    starts.in.df.are.0based = TRUE
)

"
GRanges object with 96584 ranges and 0 metadata columns:
          seqnames    ranges strand
             <Rle> <IRanges>  <Rle>
      [1]     chr1   2653136      *
      [2]     chr1   2653375      *
      [3]     chr1   2654017      *
      [4]     chr1   2654071      *
      [5]     chr1   2654137      *
      ...      ...       ...    ...
  [96580]     chrM     16449      *
  [96581]     chrM     16454      *
  [96582]     chrM     16495      *
  [96583]     chrM     16542      *
  [96584]     chrM     16565      *
  -------
  seqinfo: 711 sequences (1 circular) from hg38 genome
  "

# Process 
rownames = paste0(dss_table$chr, ":", dss_table$pos, "-", dss_table$endpos)
meth_matrix = matrix(dss_table$X, ncol = 1)
rownames(meth_matrix) = rownames
unmeth_matrix = matrix((dss_table$N - dss_table$X), ncol = 1)
rownames(unmeth_matrix) = rownames


# Help with snakemake and R
# do_something <- function(data_path, out_path, threads, myparam) {
#     # R code
# }

# do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])
