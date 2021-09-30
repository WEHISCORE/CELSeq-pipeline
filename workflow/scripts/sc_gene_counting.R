# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)

# input params -----------------------------------------------------------------


UMI_cor <- as.integer(snakemake@config[['correct_UMI_error']])
gene_fl <- as.logical(snakemake@config[['remove_low_abundance_genes']])

bc_file <- snakemake@config[['barcode_file']]
outdir <- dirname(dirname(snakemake@input[[1]]))

# use trimmed bc file if exists ------------------------------------------------

trimmed_barcodes_file <- file.path(dirname(bc_file), 'barcodes_trimmed.csv')

if (file.exists(trimmed_barcodes_file)) {
    bc_file <- trimmed_barcodes_file
}

# run sc gene counting --------------------------------------------------------

sc_gene_counting(
    outdir = outdir,
    bc_anno = bc_file,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
