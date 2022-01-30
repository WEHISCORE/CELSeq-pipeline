rule index:
    input:
        fasta=config["ref"]["fasta"],
        gtf=config["ref"]["gtf"],
    output:
        directory(config["ref"]["star_index"]),
    log:
        "logs/index.log",
    threads: 24
    resources:
        mem_mb=65536,
        runtime="0-12:0:0",
    wrapper:
        "0.77.0/bio/star/index"


rule align:
    input:
        fq1="results/trimmed/{sample}_combined.fastq.gz",
        index=config["ref"]["star_index"],
    output:
        "results/aligned/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/align_{sample}.log",
    threads: 24
    resources:
        mem_mb=65536,
        runtime="0-12:0:0",
    params:
        index=lambda w, input: "".join(os.path.splitext(input[1])),  # hack, otherwise linter complains
        extra="--outSAMtype BAM SortedByCoordinate --sjdbGTFfile {gtf} {star_params}".format(
            gtf=config["ref"]["gtf"], star_params=config["params"]["star"]
        ),
    wrapper:
        "0.77.0/bio/star/align"
