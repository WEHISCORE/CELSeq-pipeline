# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)
library(tools)

# input params -----------------------------------------------------------------

rs <- snakemake@config[['read_structure']]
bam <- snakemake@input[['bam']]
bc_file <- snakemake@input[['barcodes']]

mito <- snakemake@config[['mito']]
max_mis <- as.integer(snakemake@config[['max_barcode_mismatches']])

threads <- as.integer(snakemake@threads)
outdir <- dirname(dirname(snakemake@output[[1]]))

# process settings -------------------------------------------------------------

read_structure <- list(
    'bs1' = as.integer(rs$barcode_start_1),
    'bl1' = as.integer(rs$barcode_len_1),
    'bs2' = as.integer(rs$barcode_start_2),
    'bl2' = as.integer(rs$barcode_len_2),
    'us'  = as.integer(rs$umi_start),
    'ul'  = as.integer(rs$umi_len)
)

# fix barcode length -----------------------------------------------------------

bc <- read.delim(bc_file, sep=',')
barcode_len <- max(read_structure$bl1,
                   read_structure$bl2)
file_bc_len <- nchar(bc[1, 2])

if (file_bc_len > barcode_len) {
    print('Writing trimmed barcode file.')

    bc[, 2] <- substr(bc$barcode, 1, barcode_len)

    bc_file <- paste0(file_path_sans_ext(bc_file), '_trimmed.csv')
    write.table(bc, file = bc_file, row.names = FALSE, quote = FALSE, sep = ',')
}

# run sc demultiplexing --------------------------------------------------------

has_UMI <- read_structure$ul > 0

sc_demultiplex(
    inbam = bam,
    outdir = outdir,
    bc_anno = bc_file,
    max_mis = max_mis,
    mito = mito,
    has_UMI = has_UMI,
    nthreads = threads)
