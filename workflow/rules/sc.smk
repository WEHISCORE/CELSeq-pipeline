# single sample processing methods


rule trim_barcode:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
    output:
        "results/trimmed/{sample}_combined.fastq.gz",
    log:
        "logs/trim-barcode_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime="0-2:0:0",
    script:
        "../scripts/sc_trim_barcode.R"


rule exon_mapping:
    input:
        bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/exon_mapping/{sample}_exon_mapping.bam",
    log:
        "logs/exon_mapping_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: 16
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    script:
        "../scripts/sc_exon_mapping.R"


rule sc_demultiplex:
    input:
        bam="results/exon_mapping/{sample}_exon_mapping.bam",
    output:
        "results/sc_demultiplex/{sample}/stat/overall_stat.csv",
    log:
        "logs/sc_demultiplex_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime="0-12:0:0",
    script:
        "../scripts/sc_demultiplex.R"


rule sc_gene_counting:
    input:
        "results/sc_demultiplex/{sample}/stat/overall_stat.csv",
    output:
        "results/sc_demultiplex/{sample}/gene_count.csv",
    log:
        "logs/sc_gene_counting_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        runtime="0-6:0:0",
    script:
        "../scripts/sc_gene_counting.R"


rule create_sce_object:
    input:
        "results/sc_demultiplex/{sample}/gene_count.csv",
    output:
        "results/sc_demultiplex/{sample}/sce.rds",
    log:
        "logs/create_sce_object_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: 1
    resources:
        mem_mb=16384,
        runtime="0-6:0:0",
    script:
        "../scripts/create_sce_object.R"
