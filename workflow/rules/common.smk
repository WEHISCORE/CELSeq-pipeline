import pandas as pd
import numpy as np
import os
import sys
import yaml
from glob import iglob

# ------------- globals ------------
READENDS = ["R1", "R2"]

# ------------- load cluster config ------------
with open("config/cluster.yaml", "r") as stream:
    try:
        cluster = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc, file=sys.stderr)

# ------------- set up samples ------------
process_from_bcl = bool(config["process_from_bcl"])
run_qc = bool(config["run_qc"])
demux_tool = config["demux_tool"].lower() if config["demux_tool"] else None

if process_from_bcl:

    if demux_tool == "bcl2fastq":
        # process sample sheet for bcl2fastq
        from sample_sheet import SampleSheet

        bcl_sample_sheet = SampleSheet(config["bcl_sample_sheet"])
        samples = [
            "%s_S%d" % (s.sample_name, idx + 1)
            for idx, s in enumerate(bcl_sample_sheet.samples)
        ]
        lanes = int(config["lanes"])
        lanes = ["_L%s" % str(lane).zfill(3) for lane in range(1, lanes + 1)]

    elif demux_tool == "bcl-convert":
        # we have to parse this by hand because bcl-convert uses
        # v2 sample sheets, which the python module does not support
        sample_info = []
        with open(config["bcl_sample_sheet"], "r") as ss:
            found_sample_info = False
            lines = ss.readlines()
            sample_info_idx = np.where(np.array(lines) == "[BCLConvert_Data]\n")[0]

            if len(sample_info_idx) == 0:
                print(
                    "Could not find [BCLConvert_Data] for demultiplexing. Check your sample sheet.",
                    file=sys.stderr,
                )
                sys.exit()

            for line in lines[sample_info_idx[0] + 1 :]:
                if line == "\n":
                    break
                sample_info.append(line.rstrip())

            sample_info = [s.split(",") for s in sample_info]
            sample_info = pd.DataFrame(sample_info[1:], columns=sample_info[0])
            samples = [
                "%s_S%d" % (s, idx + 1) for idx, s in enumerate(sample_info.Sample_ID)
            ]

            # make lane variable empty if there are no lanes
            cols = [col.lower() for col in sample_info.columns]
            if "lane" not in cols:
                lanes = [""]
            else:
                lanes = int(config["lanes"])
                lanes = ["_L%s" % str(lane).zfill(3) for lane in range(1, lanes + 1)]


else:
    # in this case, we expect the fastq files to already exist
    fqs = iglob("fastq/*_R1.fastq.gz")

    # extract basename of full file path
    base = [os.path.basename(i) for i in fqs]

    # extract sample name
    samples = []
    for f in base:
        samples.append(f.split("_R1")[0])

# ------------- set up barcode annotation file ------------
barcode_exists = os.path.exists(config["barcode_file"])
sample_sheet_exists = os.path.exists(config["sample_sheet"])

if not barcode_exists and sample_sheet_exists:
    # if barcode file doesn't exist, create it from sample sheet
    sample_sheet = pd.read_csv(config["sample_sheet"], sep=",")
    barcode_col = ["primer_name" in col for col in sample_sheet.columns]
    barcode_col = sample_sheet.columns[barcode_col].values

    if len(barcode_col) != 1:
        print("Multiple barcode fields detected.", file=sys.stderr)
        sys.exit()

    cell_ids = sample_sheet.iloc[:, 0].values
    barcodes = sample_sheet[barcode_col].iloc[:, 0].values

    if len(barcodes) != len(np.unique(barcodes)):
        print("Barcodes in sample sheet are not unique.", file=sys.stderr)
        sys.exit()

    barcode_df = pd.DataFrame.from_dict({"cell_id": cell_ids, "barcode": barcodes})
    barcode_df.to_csv(config["barcode_file"], sep=",", index=False)

elif not barcode_exists and not sample_sheet_exists:
    print("Neither barcode file or sample sheet was found!", file=sys.stderr)
    sys.exit()

# ------------- output functions ------------
def get_bcl2fastq_output():
    bcl2fastq_output = expand(
        "results/bcl_output/{sample}{lane}_{readend}_001.fastq.gz",
        sample=samples,
        lane=lanes,
        readend=READENDS,
    )
    return bcl2fastq_output


def get_mergelanes_output():
    mergelanes_output = expand(
        "fastq/{sample}_{readend}.fastq.gz", sample=samples, readend=READENDS
    )
    return mergelanes_output


def get_fastqc_output():
    fastqc_output = expand(
        "results/fastQC/{sample}_{readends}_fastqc.html",
        sample=samples,
        readends=READENDS,
    )
    return fastqc_output


def get_fastqscreen_output():
    fastqscreen_output = expand(
        "results/fastqScreen/{sample}_{readends}_screen.html",
        sample=samples,
        readends=READENDS,
    )
    return fastqscreen_output


def get_trim_barcode_output():
    trim_barcode_output = expand(
        "results/trimmed/{sample}_combined.fastq.gz", sample=samples
    )
    return trim_barcode_output


def get_align_output():
    align_output = expand(
        "results/aligned/{sample}/Aligned.sortedByCoord.out.bam", sample=samples
    )
    return align_output


def get_exon_mapping_output():
    exon_mapping_output = expand(
        "results/exon_mapping/{sample}_exon_mapping.bam", sample=samples
    )
    return exon_mapping_output


def get_sc_demultiplex_output():
    sc_demultiplex_output = expand(
        "results/sc_demultiplex/{sample}/stat/overall_stat.csv",
        sample=samples,
    )
    return sc_demultiplex_output


def get_sc_gene_counting_output():
    sc_gene_counting_output = expand(
        "results/sc_demultiplex/{sample}/gene_count.csv",
        sample=samples,
    )
    return sc_gene_counting_output


def get_create_sce_object_output():
    create_sce_object_output = expand(
        "results/sc_demultiplex/{sample}/sce.rds",
        sample=samples,
    )
    return create_sce_object_output


def get_create_report_output():
    create_report_output = expand(
        "results/sc_demultiplex/{sample}/report.html",
        sample=samples,
    )
    return create_report_output
