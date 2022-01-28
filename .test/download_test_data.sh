#!/bin/sh

scpipe_extdata=https://raw.githubusercontent.com/LuyiTian/scPipe/master/inst/extdata

for file in barcode_anno.csv simu_R1.fastq.gz simu_R2.fastq.gz; do
	wget --no-check-certificate ${scpipe_extdata}/$file
done

mkdir -p metadata && mv barcode_anno.csv metadata
mkdir -p fastq && mv *fastq.gz fastq
