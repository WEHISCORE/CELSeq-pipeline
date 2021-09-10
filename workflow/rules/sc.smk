# single cell processing methods

rule trim_barcode:
    input:
        r1='fastq/{cell}_R1.fastq.gz',
        r2='fastq/{cell}_R2.fastq.gz',
    conda:
        '../envs/scpipe.yaml'
    threads:
        1
    resources:
        mem_mb=8192,
        runtime='0-2:0:0',
    log:
        'logs/trim-barcode_{cell}.log'
    output:
        'results/trimmed/{cell}_combined.fastq.gz'
    script:
        '../scripts/sc_trim_barcode.R'
