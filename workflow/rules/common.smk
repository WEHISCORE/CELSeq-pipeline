import pandas as pd
import os
from glob import iglob

#------------- globals ------------
READENDS = ['R1', 'R2']

#------------- set up samples ------------
process_from_bcl = bool(config['process_from_bcl'])

if process_from_bcl:
    # process sample sheet for bcl2fastq
    from sample_sheet import SampleSheet

    sample_sheet = SampleSheet(config['sample_sheet'])
    samples = ['%s_S%d' % (s.sample_name, idx + 1) for idx, s in enumerate(sample_sheet.samples)]

    lanes = int(config['lanes'])
    lanes = ['L%s' % str(lane).zfill(3) for lane in range(1, lanes + 1)]

else:
    # in this case, we expect the fastq files to already exist
    fqs = iglob('fastq/*_R1.fastq.gz')

    # extract basename of full file path
    base = [os.path.basename(i) for i in fqs]

    # extract sample name
    samples = []
    for f in base:
        samples.append(f.split('_R1')[0])

#------------- get cells from barcode files ------------
barcodes = pd.read_csv(config['barcode_file'])
cells = barcodes.iloc[:, 0].values

#------------- output functions ------------
def get_bcl2fastq_output():
    bcl2fastq_output = expand(
            'results/bcl_output/{sample}_{lane}_{readend}_001.fastq.gz',
            sample=samples,
            lane=lanes,
            readend=READENDS
    )
    return bcl2fastq_output

def get_mergelanes_output():
    mergelanes_output = expand(
            'fastq/{sample}_{readend}.fastq.gz',
            sample=samples,
            readend=READENDS
        )
    return mergelanes_output

def get_fastqc_output():
    fastqc_output = expand(
        'results/fastQC/{sample}_{readends}_fastqc.html',
        sample=samples,
        readends=READENDS
    )
    return fastqc_output

def get_fastqscreen_output():
    fastqscreen_output = expand(
        'results/fastqScreen/{sample}_{readends}_screen.html',
         sample=samples,
         readends=READENDS
    )
    return fastqscreen_output

def get_trim_barcode_output():
    trim_barcode_output = expand(
        'results/trimmed/{sample}_combined.fastq.gz',
         sample=samples
    )
    return trim_barcode_output

def get_align_output():
    align_output = expand(
        'results/aligned/{sample}/Aligned.sortedByCoord.out.bam',
         sample=samples
    )
    return align_output

def get_exon_mapping_output():
    exon_mapping_output = expand(
        'results/exon_mapping/{sample}_exon_mapping.bam',
         sample=samples
    )
    return exon_mapping_output

def get_sc_demultiplex_output():
    sc_demultiplex_output = expand(
        'results/sc_demultiplex/{sample}/stat/overall_stat.csv',
         sample=samples,
    )
    return sc_demultiplex_output

def get_sc_gene_counting_output():
    sc_gene_counting_output = expand(
        'results/sc_demultiplex/{sample}/gene_count.csv',
         sample=samples,
    )
    return sc_gene_counting_output

def get_create_sce_object_output():
    create_sce_object_output = expand(
        'results/sc_demultiplex/{sample}/sce.rds',
         sample=samples,
    )
    return create_sce_object_output
