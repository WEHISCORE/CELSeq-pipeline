# fastQC, fastqScreen and multiQC

rule fastQC:
    input:
        r1='fastq/{cell}_R1.fastq.gz',
        r2='fastq/{cell}_R2.fastq.gz',
    envmodules:
        'fastqc/0.11.8',
    threads:
        2
    resources:
        mem_mb=8192,
        runtime='0-12:0:0',
    log:
        'logs/fastQC_{cell}.log'
    output:
        report1='results/fastQC/{cell}_R1_fastqc.html',
        report2='results/fastQC/{cell}_R2_fastqc.html',
    shell:
        '''
        fastqc -t {threads} {input.r1} {input.r2} -o results/fastQC
        '''

rule fastqScreen:
    input:
        r2='fastq/{cell}_R2.fastq.gz',
    conda:
        '../envs/fastqscreen.yaml'
    envmodules:
        'bowtie2/2.3.4.1',
    threads:
        8
    resources:
        mem_mb=8192,
        runtime='1-0:0:0',
    log:
        'logs/fastqScreen_{cell}.log'
    output:
        report3='results/fastqScreen/{cell}_R2_screen.html'
    shell:
        '''
        fastq_screen --conf {config[FastqScreen_config]} \
            --outdir results/fastqScreen \
            {input.r2}
        '''

rule multiQC:
    input:
        report3=expand(
            'results/fastqScreen/{cell}_R2_screen.html',
            cell=cells
        ),
        report1=expand(
            'results/fastQC/{cell}_R1_fastqc.html',
            cell=cells
        ),
        report2=expand(
            'results/fastQC/{cell}_R2_fastqc.html',
            cell=cells
        ),
    envmodules:
        'MultiQC/1.10.1',
    resources:
        mem_mb=65536,
        runtime='0-12:0:0'
    log:
        'logs/multiQC.log'
    output:
        'results/qc_metrics/multiqc_report.html'
    shell:
        '''
        multiqc results/fastQC/ results/fastqScreen/ \
            -o results/qc_metrics -f
        '''
