# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)
library(SingleCellExperiment)

# input params -----------------------------------------------------------------

fs = snakemake@config[['filter_settings']]
rs = snakemake@config[['read_structure']]

swap_reads <- as.logical(rs$barcode_in_r1)

if (swap_reads) {
    fq_R1 <- snakemake@input[['r2']]
    fq_R2 <- snakemake@input[['r1']]
} else {
    fq_R1 <- snakemake@input[['r1']]
    fq_R2 <- snakemake@input[['r2']]
}

outfile <- snakemake@output[[1]]
outdir <- dirname(outfile)

# output files -----------------------------------------------------------------

combined_fq <- file.path(outdir, gsub('R[12]', 'combined', basename(fq_R1)))

# process settings -------------------------------------------------------------

filter_settings <- list(
    rmlow = as.logical(fs$remove_lowqual_reads),
    rmN   = as.logical(fs$remove_N_reads),
    minq  = as.integer(fs$min_base_q),
    numbq = as.integer(fs$max_bases_below_qual)
)

read_structure <- list(
    'bs1' = as.integer(rs$barcode_start_1),
    'bl1' = as.integer(rs$barcode_len_1),
    'bs2' = as.integer(rs$barcode_start_2),
    'bl2' = as.integer(rs$barcode_len_2),
    'us'  = as.integer(rs$umi_start),
    'ul'  = as.integer(rs$umi_len)
)

# run scpipe's trimming --------------------------------------------------------

sc_trim_barcode(combined_fq,
                fq_R1,
                fq_R2,
                read_structure = read_structure,
                filter_settings = filter_settings)
