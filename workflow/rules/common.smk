import pandas as pd
import os
from glob import iglob

#------------- fastqs must be under directory called fastq ------------
fqs = iglob('fastq/*_R1.fastq.gz')

# Extract basename of full file path
base = [os.path.basename(i) for i in fqs]

#------------- extract fastq sample name ------------
cells = []
for f in base:
    cells.append(f.split('_R1')[0])

print(
'''
INPUTS:
#-------------------- SNAKEMAKE INPUTS --------------------
%s
#----------------------------------------------------------
''' % ', '.join(cells)
)

#------------- globals ------------
READENDS = ['R1', 'R2']

#------------- output functions ------------
def get_fastqc_output():
    fastqc_output = expand(
        'results/fastQC/{cell}_{readends}_fastqc.html',
        cell=cells,
        readends=READENDS
    )
    return fastqc_output


def get_fastqscreen_output():
    fastqscreen_output = expand(
        'results/fastqScreen/{cell}_R2_screen.html',
         cell=cells
    )
    return fastqscreen_output

def get_trim_barcode_output():
    trim_barcode_output = expand(
        'results/trimmed/{cell}_combined.fastq.gz',
         cell=cells
    )
    return trim_barcode_output
