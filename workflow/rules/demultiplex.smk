rule bcl2fastq:
    input:
        samplesheet=config['sample_sheet']
    output:
        fastq=expand(
            'results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz',
            sample=samples,
            lane=lanes,
            readend=READENDS
        )
    envmodules:
        'bcl2fastq/2.20.0'
    threads:
        16
    resources:
        mem_mb=98304,
        runtime='0-12:0:0'
    log:
        'logs/bcl2fastq.log'
    shell:
        '''
        bcl2fastq \
            --runfolder-dir {config[bcl_input]} \
            --output-dir results/bcl_output \
            --sample-sheet {input} \
            {config[params][bcl2fastq]}
        '''

rule mergelanes:
    input:
        fastq=[
            'results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz'.format(
                sample=sample,
                lane=lane,
                readend=readend
            ) for sample, lane, readend in zip(samples, lanes, READENDS)
        ]
    output:
        'fastq/{sample}_{readend}.fastq.gz',
    resources:
        mem_mb=12288,
        runtime='0-2:0:0'
    log:
        'logs/mergelanes_{sample}_{readend}.log'
    shell:
        'cat results/bcl_output/{wildcards.sample}_*_{wildcards.readend}_001.fastq.gz > {output}'
