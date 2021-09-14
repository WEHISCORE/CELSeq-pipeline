import pandas as pd
import os
from glob import iglob

from sample_sheet import SampleSheet

sample_sheet = SampleSheet(config['sample_sheet'])
samples = ['%s_S%d' % (s.sample_name, idx + 1) for idx, s in enumerate(sample_sheet.samples)]

lanes = int(config['lanes'])
lanes = ['L%s' % str(lane).zfill(3) for lane in range(1, lanes + 1)]

##------------- fastqs must be under directory called fastq ------------
#fqs = iglob('fastq/*_R1.fastq.gz')
#
## Extract basename of full file path
#base = [os.path.basename(i) for i in fqs]
#
##------------- extract fastq sample name ------------
#samples = []
#for f in base:
#    samples.append(f.split('_R1')[0])
#
#print(
#'''
#INPUTS:
##-------------------- SNAKEMAKE INPUTS --------------------
#%s
##----------------------------------------------------------
#''' % ', '.join(samples)
#)

#------------- globals ------------
READENDS = ['R1', 'R2']

#------------- output functions ------------
def get_bcl2fastq_output():
    bcl2fastq_output = expand(
            'bcl_output/{sample}_{lane}_{readend}_001.fastq.gz',
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
        'results/fastqScreen/{sample}_R2_screen.html',
         sample=samples
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
        'results/aligned/{sample}/AlignedSortedByCoord.out.bam',
         sample=samples
    )
    return align_output
