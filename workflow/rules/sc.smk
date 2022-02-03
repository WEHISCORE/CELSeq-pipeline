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
    threads: cluster["trimbarcode"]["threads"]
    resources:
        mem_mb=cluster["trimbarcode"]["mem_mb"],
        runtime=cluster["trimbarcode"]["runtime"],
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
    threads: cluster["exonmapping"]["threads"]
    resources:
        mem_mb=cluster["exonmapping"]["mem_mb"],
        runtime=cluster["exonmapping"]["runtime"],
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
    threads: cluster["scdemux"]["threads"]
    resources:
        mem_mb=cluster["scdemux"]["mem_mb"],
        runtime=cluster["scdemux"]["runtime"],
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
    threads: cluster["scgenecounting"]["threads"]
    resources:
        mem_mb=cluster["scgenecounting"]["mem_mb"],
        runtime=cluster["scgenecounting"]["runtime"],
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
    threads: cluster["createsceobject"]["threads"]
    resources:
        mem_mb=cluster["createsceobject"]["mem_mb"],
        runtime=cluster["createsceobject"]["runtime"],
    script:
        "../scripts/create_sce_object.R"


rule create_report:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        rds="results/sc_demultiplex/{sample}/sce.rds",
        trimfq="results/trimmed/{sample}_combined.fastq.gz",
        align_bam="results/aligned/{sample}/Aligned.sortedByCoord.out.bam",
        map_bam="results/exon_mapping/{sample}_exon_mapping.bam",
    output:
        "results/sc_demultiplex/{sample}/report.html",
    log:
        "logs/create_report_{sample}.log",
    conda:
        "../envs/scpipe.yaml"
    threads: cluster["createsceobject"]["threads"]
    resources:
        mem_mb=cluster["createsceobject"]["mem_mb"],
        runtime=cluster["createsceobject"]["runtime"],
    script:
        "../scripts/create_report.R"
