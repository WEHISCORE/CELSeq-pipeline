rule align:
    input:
        fq1='results/trimmed/{sample}_combined.fastq.gz',
    threads:
        24
    resources:
        mem_mb=65536,
        runtime='0-12:0:0',
    log:
        'logs/align_{sample}.log'
    output:
        'results/aligned/{sample}/Aligned.sortedByCoord.out.bam'
    params:
        index=config['ref']['star_index'],
        extra='--outSAMtype BAM SortedByCoordinate --sjdbGTFfile {gtf} {star_params}'.format(
            gtf=config['ref']['gtf'],
            star_params=config['params']['star']
        ),
    wrapper:
        '0.77.0/bio/star/align'
