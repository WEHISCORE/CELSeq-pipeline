# fastQC, fastqScreen and multiQC


rule fastQC:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
    output:
        r1="results/fastQC/{sample}_R1_fastqc.html",
        r2="results/fastQC/{sample}_R2_fastqc.html",
    log:
        "logs/fastQC_{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    resources:
        mem_mb=8192,
        runtime="0-12:0:0",
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -o results/fastQC
        """


rule fastqScreen:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        config_file=config["FastqScreen_config"],
    conda:
        "../envs/fastqscreen.yaml"
    threads: 8
    resources:
        mem_mb=8192,
        runtime="1-0:0:0",
    log:
        "logs/fastqScreen_{sample}.log",
    output:
        report3="results/fastqScreen/{sample}_R1_screen.html",
        report4="results/fastqScreen/{sample}_R2_screen.html",
    shell:
        """
        fastq_screen --conf {input.config_file} \
            --outdir results/fastqScreen \
            {input.r1} {input.r2}
        """


rule multiQC:
    input:
        report1=expand("results/fastQC/{sample}_R1_fastqc.html", sample=samples),
        report2=expand("results/fastQC/{sample}_R2_fastqc.html", sample=samples),
        report3=expand("results/fastqScreen/{sample}_R1_screen.html", sample=samples),
        report4=expand("results/fastqScreen/{sample}_R2_screen.html", sample=samples),
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        "MultiQC/1.10.1",
    resources:
        mem_mb=65536,
        runtime="0-12:0:0",
    log:
        "logs/multiQC.log",
    output:
        "results/qc_metrics/multiqc_report.html",
    shell:
        """
        multiqc results/fastQC/ results/fastqScreen/ \
            -o results/qc_metrics -f
        """
