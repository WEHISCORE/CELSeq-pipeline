rule bcl2fastq:
    input:
        samplesheet=config["bcl_sample_sheet"],
    output:
        fastq=expand(
            "results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz",
            sample=samples,
            lane=lanes,
            readend=READENDS,
        ),
    log:
        "logs/bcl2fastq.log",
    envmodules:
        "bcl2fastq/2.20.0",
    threads: cluster["bcl2fastq"]["threads"]
    resources:
        mem_mb=cluster["bcl2fastq"]["mem_mb"],
        runtime=cluster["bcl2fastq"]["runtime"],
    shell:
        """
        bcl2fastq \
            --runfolder-dir {config[bcl_input]} \
            --output-dir results/bcl_output \
            --sample-sheet {input} \
            {config[params][bcl2fastq]}
        """


rule mergelanes:
    input:
        fastq=[
            "results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz".format(
                sample=sample, lane=lane, readend=readend
            )
            for sample, lane, readend in zip(samples, lanes, READENDS)
        ],
    output:
        "fastq/{sample}_{readend}.fastq.gz",
    log:
        "logs/mergelanes_{sample}_{readend}.log",
    threads: cluster["mergelanes"]["threads"]
    resources:
        mem_mb=cluster["mergelanes"]["mem_mb"],
        runtime=cluster["mergelanes"]["runtime"],
    shell:
        "cat results/bcl_output/{wildcards.sample}_*_{wildcards.readend}_001.fastq.gz > {output}"
