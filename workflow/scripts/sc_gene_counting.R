# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)
library(tools)

# input params -----------------------------------------------------------------


UMI_cor <- as.integer(snakemake@config[['correct_UMI_error']])
gene_fl <- as.logical(snakemake@config[['remove_low_abundance_genes']])

outdir <- dirname(dirname(snakemake@input[['stats']]))
bc_file <- snakemake@input[['barcodes']]

# use trimmed bc file if exists ------------------------------------------------

trimmed_barcodes_file <- paste0(file_path_sans_ext(bc_file), '_trimmed.csv')
if (file.exists(trimmed_barcodes_file)) {
    bc_file <- trimmed_barcodes_file
}

# run sc gene counting --------------------------------------------------------

sc_gene_counting(
    outdir = outdir,
    bc_anno = bc_file,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
