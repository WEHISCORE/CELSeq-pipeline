# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)

# input params -----------------------------------------------------------------

datadir <- dirname(snakemake@input[[1]])
datadir <- file.path(getwd(), datadir) # need absolute path

organism <- snakemake@config[['organism']]
gene_id_type <- snakemake@config[['gene_id_type']]
pheno_data_file <- snakemake@config[['sample_sheet']]

outfile <- snakemake@output[[1]]

if (file.exists(pheno_data_file)) {
    pheno_data <- read.delim(pheno_data_file, sep=',', row.names = 1)
} else {
    pheno_data <- NULL
    print('WARNING: Sample sheet not found. Unable to add phenotype data.')
}

# create SCE objects from dir --------------------------------------------------
sce <- create_sce_by_dir(
        datadir = datadir,
        organism = organism,
        gene_id_type = gene_id_type,
        pheno_data = pheno_data,
        report = FALSE)

# save counts as sparse matrix
assay(sce, withDimnames = FALSE) <- as(
    assay(sce, withDimnames = FALSE),
    'dgCMatrix')

# write R output ---------------------------------------------------------------
saveRDS(
    sce,
    outfile,
    compress = 'xz'
)
