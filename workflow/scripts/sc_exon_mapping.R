# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)

# input params -----------------------------------------------------------------

rs <- snakemake@config[['read_structure']]
inbam <- snakemake@input[['bam']]
strand <- as.logical(snakemake@config[['strand']])
gtf <- snakemake@config[['ref']][['gtf']]
threads <- as.integer(snakemake@threads)

# output files -----------------------------------------------------------------

outbam <- snakemake@output[[1]]

# process settings -------------------------------------------------------------

read_structure <- list(
    'bs1' = as.integer(rs$barcode_start_1),
    'bl1' = as.integer(rs$barcode_len_1),
    'bs2' = as.integer(rs$barcode_start_2),
    'bl2' = as.integer(rs$barcode_len_2),
    'us'  = as.integer(rs$umi_start),
    'ul'  = as.integer(rs$umi_len)
)

# perform exon mapping ---------------------------------------------------------

bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul

sc_exon_mapping(
    inbam = inbam,
    outbam = outbam,
    annofn = gtf,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = strand,
    nthreads = threads)
