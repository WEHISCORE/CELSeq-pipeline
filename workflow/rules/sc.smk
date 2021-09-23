# single sample processing methods

rule trim_barcode:
    input:
        r1='fastq/{sample}_R1.fastq.gz',
        r2='fastq/{sample}_R2.fastq.gz',
    conda:
        '../envs/scpipe.yaml'
    threads:
        1
    resources:
        mem_mb=8192,
        runtime='0-2:0:0',
    log:
        'logs/trim-barcode_{sample}.log'
    output:
        'results/trimmed/{sample}_combined.fastq.gz'
    script:
        '../scripts/sc_trim_barcode.R'

rule exon_mapping:
    input:
        bam='results/aligned/{sample}/Aligned.sortedByCoord.out.bam'
    conda:
        '../envs/scpipe.yaml'
    resources:
        mem_mb=32768,
        runtime='0-12:0:0',
    log:
        'logs/exon_mapping_{sample}.log'
    output:
        'results/exon_mapping/{sample}_exon_mapping.bam'
    script:
        '../scripts/sc_exon_mapping.R'
