# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)
library(tools)

options(error = quote({dump.frames(); save.image(file = "last.dump.rda")}))

# input params -----------------------------------------------------------------

counts_file <- snakemake@input[['counts']]
datadir <- file.path(getwd(), dirname(counts_file)) # need absolute path

rs <- snakemake@config[['read_structure']]
fs = snakemake@config[['filter_settings']]

outfq <- snakemake@input[['trimfq']]
align_bam <- snakemake@input[['align_bam']]
map_bam <- snakemake@input[['map_bam']]

genome_index <- snakemake@config[['ref']][['star_index']]
gtf <- snakemake@config[['ref']][['gtf']]
ercc <- system.file("extdata", "ERCC92_anno.gff3", package = "scPipe")
strand <- as.logical(snakemake@config[['strand']])
bc_file <- snakemake@config[['barcode_file']]

max_mis <- as.integer(snakemake@config[['max_barcode_mismatches']])
UMI_cor <- as.integer(snakemake@config[['correct_UMI_error']])
gene_fl <- as.logical(snakemake@config[['remove_low_abundance_genes']])
organism <- snakemake@config[['organism']]
gene_id_type <- snakemake@config[['gene_id_type']]

fix_chr <- FALSE

# process params ---------------------------------------------------------------

swap_reads <- as.logical(rs$barcode_in_r1)
if (swap_reads) {
    fq_R1 <- snakemake@input[['r2']]
    fq_R2 <- snakemake@input[['r1']]
} else {
    fq_R1 <- snakemake@input[['r1']]
    fq_R2 <- snakemake@input[['r2']]
}

read_structure <- list(
    'bs1' = as.integer(rs$barcode_start_1),
    'bl1' = as.integer(rs$barcode_len_1),
    'bs2' = as.integer(rs$barcode_start_2),
    'bl2' = as.integer(rs$barcode_len_2),
    'us'  = as.integer(rs$umi_start),
    'ul'  = as.integer(rs$umi_len)
)

filter_settings <- list(
    rmlow = as.logical(fs$remove_lowqual_reads),
    rmN   = as.logical(fs$remove_N_reads),
    minq  = as.integer(fs$min_base_q),
    numbq = as.integer(fs$max_bases_below_qual)
)

# assume plate name is first part of sample name
plate <- strsplit(basename(outfq), '_')[[1]][1]

anno <- c(gtf, ercc)

# modify test data, otherwise report crashes -----------------------------------

test_dir <- '.test/results/sc_demultiplex/simu'
if (grepl(test_dir, datadir)) {
    cell_stat_file <- file.path(datadir, 'stat/cell_stat.csv')
    cell_stat <- read.delim(cell_stat_file, sep=',')

    # fill mapped_to_exon field with dummy records, to prevent report crashing
    cell_stat$mapped_to_exon <- cell_stat$mapped_to_ERCC

    write.csv(cell_stat, file = cell_stat_file, quote = FALSE, row.names = FALSE)

    # replace output data with scPipe's pregenerated data for QC reports
    data('sc_sample_data')
    sc_sample_data <- data.frame(sc_sample_data)
    cell_names <- colnames(sc_sample_data)
    sc_sample_data$gene_id <- rownames(sc_sample_data)
    counts <- sc_sample_data[,c('gene_id', cell_names)]

    write.csv(counts, file = counts_file, quote = FALSE, row.names = FALSE)
}

# create report from demultiplexed output --------------------------------------

# create the report Rmd (this will fail but we will fix it)
try(create_report(
    sample_name = plate,
    outdir = datadir,
    r1 = fq_R1,
    r2 = fq_R2,
    outfq = outfq,
    read_structure = read_structure,
    filter_settings = filter_settings,
    align_bam = align_bam,
    genome_index = genome_index,
    map_bam = map_bam,
    exon_anno = anno,
    stnd = strand,
    fix_chr = fix_chr,
    barcode_anno = bc_file,
    max_mis = max_mis,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl,
    organism = organism,
    gene_id_type = gene_id_type
))

# fix the report, based on code from Peter Hickey
# https://github.com/WEHISCORE/C122_Clucas/blob/14a13ec1b10b2d42976fd08784e5b47ea2429d59/code/scPipe.R
tmp <- readLines(file.path(datadir, "report.Rmd"))

# we want to quit after the mapping statistics have been rendered
last_line <- grep("plot_mapping", tmp, value = FALSE) %>% max()
tmp <- c(tmp[1:last_line],
        "knitr::knit_exit()",
        tmp[last_line:length(tmp)])

writeLines(tmp, file.path(datadir, "report.Rmd"))
knitr::wrap_rmd(
    file = file.path(datadir, "report.Rmd"),
    width = 120,
    backup = NULL)
rmarkdown::render(
    input = file.path(datadir, "report.Rmd"),
    output_file = file.path(datadir, "report.html"),
    knit_root_dir = ".")

# get rid of .nb
file.rename(
    from = file.path(datadir, "report.nb.html"),
    to = file.path(datadir, "report.html")
)
