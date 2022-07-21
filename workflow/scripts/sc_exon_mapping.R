# Init log file ----------------------------------------------------------------
log <- file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

library(scPipe)
library(tools)

# input params -----------------------------------------------------------------

rs <- snakemake@config[['read_structure']]
inbam <- snakemake@input[['bam']]
strand <- as.logical(snakemake@config[['strand']])
gtf <- snakemake@config[['ref']][['gtf']]
ercc <- system.file("extdata", "ERCC92_anno.gff3", package = "scPipe")
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

# if gene reference is in GFF format, create SAF, based on code by Peter Hickey
# https://github.com/WEHISCORE/C122_Clucas/blob/14a13ec1b10b2d42976fd08784e5b47ea2429d59/code/scPipe.R
anno <- NULL
if (file_ext(gtf) == "gff") {
    gr <- rtracklayer::import(gtf)
    exon_gr <- gr[gr$type == "exon"]
    saf <- data.frame(
        GeneID = sapply(exon_gr$Parent, "[[", 1),
        Chr = seqnames(exon_gr),
        Start = start(exon_gr),
        End = end(exon_gr),
        Strand = strand(exon_gr))

    ercc_gr <- rtracklayer::import(ercc)
    ercc_exon_gr <- ercc_gr[ercc_gr$type == "exon"]
    ercc_saf <- data.frame(
        GeneID = ercc_exon_gr$Name,
        Chr = seqnames(ercc_exon_gr),
        Start = start(ercc_exon_gr),
        End = end(ercc_exon_gr),
        Strand = strand(ercc_exon_gr))
    anno <- rbind(saf, ercc_saf)
} else {
    anno <- c(gtf, ercc)
}

# perform exon mapping ---------------------------------------------------------

bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul

sc_exon_mapping(
    inbam = inbam,
    outbam = outbam,
    annofn = anno,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = strand,
    nthreads = threads)
