# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)

# input params -----------------------------------------------------------------

dirs <- lapply(snakemake@input, dirname)
print(dirs)

organism <- snakemake@config[['organism']]
gene_id_type <- snakemake@config[['gene_id_type']]

threads <- as.integer(snakemake@threads)
outfile <- snakemake@output[[1]]

# create SCE objects from dirs -------------------------------------------------
sce_list <- mclapply(dirs, function(dir) {
                create_sce_by_dir(
                    datadir = dir,
                    organism = organism,
                    pheno_data = NULL,
                    report = FALSE)
                },
                mc.cores = threads)

sce <- Reduce(function(x, y) .combine(x, y, rowData_by = NULL), sce_list)

# write R output ---------------------------------------------------------------
saveRDS(
    sce,
    outfile
)
