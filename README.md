# CELSeq Pipeline

A [snakemake](https://snakemake.readthedocs.io) pipeline for processing data generated using the [CEL-Seq](https://www.sciencedirect.com/science/article/pii/S2211124712002288) protocol. Takes BCL or fastq input files and generates a single-cell experiment object using [scpipe](https://github.com/LuyiTian/scPipe). Also performs QC using [Fastq Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), collated in a [MultiQC](https://multiqc.info/) report. In principle, the pipeline can be used for a range of singe-cell protocols.

### Installation ###

The only prerequisite is [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install snakemake, you will need to install a Conda-based Python3 distribution. For this, [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) is recommended. Once mamba is installed, snakemake can be installed like so:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Now activate the snakemake environment (you'll have to do this every time you want to run the pipeline):

```
conda activate snakemake
```

Now clone the repository:

```
git clone https://github.com/WEHISCORE/CELSeq-pipeline.git
cd CELSeq-pipeline
```

### Testing ###

If you would like to test the pipeline, first download the test data:

```
cd .test
./download_test_data.sh
```

Now run as follows:

```
snakemake --use-conda --conda-frontend mamba --cores 1
```

### Configuration ###

The configuration file is found under `config/config.yaml` and the config file for FastQ Screen is found under `config/fastq_screen.conf`. Please carefully go through these settings. The main settings to consider will be

- `process_from_bcl` -- set this to `True` only if converting from BCL files. If so, make sure the demultiplexing argument `bcl2fastq` is set properly (under `params`).
- `sample_sheet` -- this is the sample sheet for bcl2fastq conversion. Please check the [bcl2fastq documentation](https://sapac.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf) for more info. You can skip this if you're using fastq files.
- `barcode_file` -- contains a comma-separated file with an ID column (matching you well/cell IDs) and the corresponding barcode in the following format:
```
ID,Cell_Barcode
S1,ATATATAT
S2,GCGCGCGC
```
- `gtf` and `star_index` under `ref` -- make sure the chromosome names match for these and that you've generated an index for STAR-2.7.8, as this is the version used by the pipeline.
- `read_structure` -- ensure `barcode_in_r1` is set to `TRUE` if your barcodes are in R1 (which is standard for CEL-Seq). WEHI's modified CEL-Seq protocol uses a barcode size of 7 (`barcode_len_2` default), so set this to 8 if using a standard version of the protocol.

If you are running from BCL, make sure you put your BCL files under the `bcl_input` directory, and if running from fastqs, put them all under a `fastq` directory from where you run the pipeline.

### Running ###

Run the pipeline as follows:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 1
```

If you want to submit your jobs to the cluster using SLURM, use the following to run the pipeline:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --profile slurm --jobs 8 --cores 2
```

The pipeline will generate all results under a `results` directory. The final output will be under `results/sc_demultiplex/{sample}/sce.rds`. 
