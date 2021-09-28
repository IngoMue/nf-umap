# nf-mapping

**My little DNA read mapping workflow which I wrote to use on sequences from museum specimen**

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04)
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10)
3. Download the pipeline, create a nextflow config profile that matches your cluster set-up ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) ) and start running your own analysis! If you want to use the existing `rackham.config` remember to specify your SNIC project ID (format: `snic20XX-XX-XXX`) as well as the path to `nf-polish/environment.yml`

    ```bash
    nextflow run nf-mapping/main.nf -profile rackham --reads 'READS'
    ```
4. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```
**Make sure** that there is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored).

## Overview

The pipeline was designed to follow up after read cleaning, specifically having the output from [nf-polish](https://github.com/MozesBlom/nf-polish) in mind, but it can of course be used on any paired-end reads (unpaired and/or merged). In its current state it runs through the following steps:

* Reference sequence indexing<sup>1</sup>
* Unpaired (UNP) read alignment<sup>1</sup>
* Merged (MRG) read alignment<sup>1</sup>
* Conversion to `.bam` format<sup>2</sup>
* Sorting and indexing of `.bam` files<sup>2</sup>
* Merging UNP and MRG files as well as different libraries for the same sample (if available)<sup>2</sup>
* Perform quality control and generate a html and raw txt report for each individual with bamqc<sup>3</sup>

The pipeline uses the following tools (Numbers show which step uses which tool):

<sup>1</sup>[`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) (version 2.2.1)

<sup>2</sup>[`samtools`](http://www.htslib.org/) (version 1.13)

<sup>3</sup>[`qualimap`](http://qualimap.conesalab.org/) (version 2.2.2d)

## Input

The pipeline uses `fastq(.gz)` reads as input. All reads need to be stored in a single directory.
Information that needs to be supplied either through flags or specified within `nextflow.config` includes:
* Reference sequence directory `refdir`
* Full name of the reference sequence `refname` (e.g. `GCF_000738735.5_ASM73873v5_genomic.fna.gz`)
* Prefix of the reference sequence `refprefix` (e.g. `GCF_000738735.5_ASM73873v5_genomic`)
* Absolute path where unpaired reads are stored `read_pairs` (`path/to/*_{R1,R2}.fastq.gz`)
* Asbolute path where merged reads a stored `merged_reads` (`path/to/*_U.fastq.gz`)
* Output directory `outdir`
* Specify whether you want to use merged and/or unpaired reads for mapping `inclMrgRds` `inclUnpRds` respectively (boolean, default true)
* Specify whether the server uses an X11 system `X11` (boolean), if it doesn't, qualimap will throw an error (default false)
* Decide whether the QC report should be a full .pdf report with figures or a simple .txt file `fullreport` (boolean, default true)

Reads should follow this naming convention `{SAMPLE_ID}_{LIBRARY}_{R1,R2,U}.{format}` where `{Library}` should be an "L" followed by any number of digits (e.g `L01` or `L001`). `{SAMPLE_ID}` can be anything. `{format}` can be specified in any way as long as it gets recognised by bwa-mem2 and samtools.

An example for using this workflow would look like this:
```bash
nextflow run main.nf -profile custom --refdir /some/path/RefDir/ --refname RefID_genomic.fa --refprefix RefID_genomic --read_pairs /some/path/'*_{R1,R2}.fastq.gz' --merged_reads /some/path/'*_U.fastq.gz' --outdir /some/path/results/ --X11 true --fullreport false
```
