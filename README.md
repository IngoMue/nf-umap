# nf-μmap

**nf-μmap [/ˈmjuːmæp/] is my little DNA read mapping workflow which I wrote to use on sequences obtained from natural history museum specimen**

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 22.10.2)
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10). By default nextflow will create a conda environment inside the work folder, this can take quite some time and ressources and is done without using a cluster's queuing system which can cause errors. To circumvent this, one option is to use [`Mamba`](https://github.com/mamba-org/mamba) instead which is much faster and efficient than conda. Alternatively, it is also possible to create the environment manually using the environment.yml file and change the `conda` flag within `profile.config` from the .yml file to the path in which the environment was created.
3. Download the pipeline, edit or create a config profile for the cluster you are using ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) ) and run the workflow. If you want to use the existing `profile.config` which is written for users of Uppmax' Rackham cluster, remember to specify your naiss project ID (format: `naiss20XX-XX-XXX`) as well as the path to `nf-umap/environment.yml`

    ```bash
    nextflow run main.nf -profile custom --refdir 'RefDir' --refname 'RefID.fa' --refprefix 'RefID' --read_pairs 'READS' --merged_reads 'MERGED_READS' --outdir 'OutputDir*
    ```
   (see below for a more detailed example)
   Other useful flags when running nextflow scripts:
   `-resume` will resume a workflow using cached intermediate files if a previous run was terminated early
   `-with-report` prints an execution report summarising input flags, directories, resource usage and information about each task
4. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```
**Make sure** that there is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored).

## Overview

The pipeline was designed to follow up after read cleaning, specifically having the output from [nf-polish](https://github.com/MozesBlom/nf-polish) in mind, but it can of course be used on any paired-end reads (unpaired and/or merged). In its current state it runs through the following steps:

* Reference sequence indexing <sup>1, 2</sup>
* Unpaired read pair (PRS) alignment <sup>2 or 3</sup>
* Merged (MRG) read alignment <sup>2 or 3</sup>
* Quickcheck whether any `.sam` files are truncated <sup>1</sup>
* Conversion to `.bam` format <sup>1</sup>
* Sorting and indexing of `.bam` files <sup>1</sup>
* Merging PRS and MRG files as well as different libraries for the same sample (if available) <sup>1</sup>
* Another quickcheck to see whether the final `.bam` files are truncated <sup>1</sup>
* Perform quality control and generate a html and raw txt report for each individual with bamqc <sup>4</sup>
* Investigate damage patterns typical for aDNA/hDNA <sup>5</sup>

Read group information (read group ID, Library, Sequencing Platform, Platform unit, sample ID) is added during the alignment steps using bwa mem's `-R` flag

The pipeline uses the following tools (Numbers show which step uses which tool):

<sup>1</sup>[`samtools`](http://www.htslib.org/) (version 1.13)

<sup>2</sup>[`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) (version 2.2.1)

<sup>3</sup>[`bwa`](http://bio-bwa.sourceforge.net/) (version 0.7.17)

<sup>4</sup>[`qualimap`](http://qualimap.conesalab.org/) (version 2.2.2d)

<sup>5</sup>[`DamageProfiler`](https://github.com/Integrative-Transcriptomics/DamageProfiler) (version 0.4.9)

## Input

The pipeline uses `fastq(.gz)` reads as input. All reads need to be stored in a single directory. The reference sequence should be in uncompressed FASTA (`.fasta`, `.fa`, `.fna`, ...) format.
Information that needs to be supplied either through flags or specified within `nextflow.config` includes:
* Asbolute path to the reference sequence directory `refseq` (`path/to/GCF_000738735.5_ASM73873v5_genomic.fa`)
* Prefix of the reference sequence `refprefix` (e.g. `GCF_000738735.5_ASM73873v5_genomic`)
* Absolute path where unpaired reads are stored `read_pairs` (`path/to/'*_{R1,R2}.fastq.gz'`). It's important to include the quotations as Nextflow won't recognise the expression and only use the first matching file otherwise.
* Asbolute path where merged reads are stored `merged_reads` (`path/to/'*_U.fastq.gz'`). Again remember to include the quotation marks.
* Output directory `outdir`
* Specify whether you want to use merged and/or unpaired reads for mapping `inclMrgRds` `inclUnpRds` respectively (boolean, default true)
* Decide whether bwa's newer version bwa-mem2 (identical output, faster, but larger index files, may not run on every cluster) or the original bwa-mem should be used for mapping `usebwamem2` (boolean, default true)
* Should merging be skipped? `skipmerge` (boolean, default false) If one library per individual and only one type of reads are used, merging becomes unnecessary
* Specify the read group platform unit `rgPU` (string, default `"X12345:123:ABCDEFGH1"`) this can be e.g. the flowcell ID, which can be found in fastq files received from sequencing deliveries, the first three values in the header are @(instrumentID):(runNumber):(flowcellID) (e.g. `A00621:496:HGLYGDSX2`)
* Should unmerged bam files be copied into the output directory? `publishInterBams` (boolean, default false) Note that these are the final bam files if merging is skipped.
* Specify whether the server uses an X11 system `X11` (boolean), if it doesn't, qualimap will throw an [error](http://qualimap.conesalab.org/doc_html/faq.html#x11problem) (default false)
* Specify whether you are using Mamba `UseMamba` (boolean, default true)
* Decide whether the QC report should be a full .pdf report with figures or a simple .txt file `fullreport` (boolean, default true)
* Decide whether you want to run DamageProfiler to asses damage patterns `runDMGprof` (boolean, default true)

Reads should follow this naming convention `{SAMPLE_ID}_{LIBRARY}_{R1,R2,U}.{format}` where `{Library}` should be an "L0" followed by any number of digits (e.g `L01` or `L001`, note that this is not the sequencing lane, this information can be instead added through the `rgPU` flag), for this it is also important that the filename does not contain a "_L0" anywhere else. `{SAMPLE_ID}` can be anything. `{format}` can be specified in any way as long as it gets recognised by bwa-mem2 and samtools.

To make sure neither `.sam` nor  `.bam` files are truncated, make sure the file lists generated by the quickcheck (`00_quickcheck/*.fofn`) are empty.

Note: The default value for the sequencing platform which gets added through the read group information is set for illumina, if you sequenced your sample on a different platform, you can change the `PL:ILLUMINA` flag accordingly in each `bwa_mem_*`/`bwa_mem2_*` process within `main.nf`

An example for using this workflow would look like this:
```bash
nextflow run main.nf -profile custom --refseq /some/path/RefDir/RefID_genomic.fa --refprefix RefID_genomic --read_pairs /some/path/'*_{R1,R2}.fastq.gz' --merged_reads /some/path/'*_U.fastq.gz' --outdir /some/path/results/ --usebwamem2 false --rgPU "A00621:496:HGLYGDSX2" --X11 true --fullreport false --runDMGprof true -with-report
```
