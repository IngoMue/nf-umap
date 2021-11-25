#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\

         ================================================
         NF-MAPPING - A nextflow mapping workflow
         https://github.com/IngoMue/nf-mapping     
         Author: Ingo A. Müller
         ================================================
         |refdir       : ${params.refdir}
         |refname      : ${params.refname}
         |refprefix    : ${params.refprefix}
         |readpairs    : ${params.read_pairs}
         |mergedreads  : ${params.merged_reads}
         |outdir       : ${params.outdir}
         |
         |
         |Include merged reads?.....${params.inclMrgRds}
         |Include unpaired reads?...${params.inclUnpRds}
         |
         |X11 unix server?..........${params.X11}
         |
         |Print full QC report?.....${params.fullreport}
         ================================================
         """
         .stripIndent()

refdir_ch = Channel.fromPath(params.refdir)

refseq_ch = Channel.fromPath("${params.refdir}/${params.refname}")

unp_reads = Channel.fromPath(params.read_pairs)
        .map { [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()

mrg_reads = Channel.fromPath(params.merged_reads)
        .map { [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()


process index_ref {
    label 'High_RAM'

    publishDir "${params.outdir}/01_RefIndex", mode:'copy'
    tag "$ref_file"

    input:
      file ref_file

    output:
      file '*'

    script:
    """
      bwa-mem2 index -p $params.refprefix $ref_file
    """
}


process bwa_mem2_UNP {
    label 'High_RAM'

    tag "$sample_id"

    input:
      tuple val(sample_id), file(sample_file), val(ref_ID), file(index)

    output:
      tuple val(sample_id), file("${sample_id}_UNP.sam")

    when:
      params.inclUnpRds == true

    script:
    """
      bwa-mem2 mem $params.refprefix $sample_file -t ${task.cpus} > ${sample_id}_UNP.sam
    """
}

process bwa_mem2_MRG {
    label 'High_RAM'

    tag "$mrg_id"

    input:
      tuple val(mrg_id), file(merged_file), val(ref_ID), file(index)

    output:
      tuple val(mrg_id), file("${mrg_id}_MRG.sam")

    when:
      params.inclMrgRds == true

    script:
    """
      bwa-mem2 mem $params.refprefix $merged_file -t ${task.cpus} > ${mrg_id}_MRG.sam
    """
}


process samtools_bam_srt_idx_UNP {
    publishDir "${params.outdir}/02_bams/UNPairedReads", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_UNP_sorted.bam")
      file("${sample_id}_UNP_sorted.bam.bai")

    when:
      params.inclUnpRds == true

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_UNP.bam

      samtools sort -o ${sample_id}_UNP_sorted.bam ${sample_id}_UNP.bam

      samtools index ${sample_id}_UNP_sorted.bam
    """
}

process samtools_bam_srt_idx_MRG {
    publishDir "${params.outdir}/02_bams/MergedReads", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_MRG_sorted.bam")
      file("${sample_id}_MRG_sorted.bam.bai")

    when:
      params.inclMrgRds == true

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_MRG.bam

      samtools sort -o ${sample_id}_MRG_sorted.bam ${sample_id}_MRG.bam

      samtools index ${sample_id}_MRG_sorted.bam
    """
}


process samtools_merge {
    publishDir "${params.outdir}/03_merged_bams/", pattern: '*{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(bam_file)

    output:
      tuple val(sample_id), file("${sample_id}.bam")
      file '*.bai'

    script:
    """
      samtools merge -o ${sample_id}.bam $bam_file

      samtools index ${sample_id}.bam
    """
}

process qualimap_bamqc {
    publishDir "${params.outdir}/04_mapping_stats/${sample_id}", pattern: '*.{txt,pdf}', mode: 'move'
    tag "$sample_id"
    label 'High_RAM'

    input:
      tuple val(sample_id), file(bam_file)

    output:
      file '*'

    script:
        if (params.X11 == true && params.fullreport == false)
            """
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -c
              cp ./${sample_id}/genome_results.txt ./${sample_id}_genome_results.txt
	    """
        else if (params.X11 == false && params.fullreport == false)
            """
              unset DISPLAY
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -c
              cp ./${sample_id}/genome_results.txt ./${sample_id}_genome_results.txt
            """
        else if (params.X11 == true && params.fullreport == true)
            """
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -outfile ${sample_id}_report.pdf -c
              cp ./${sample_id}/* .
            """
        else if (params.X11 == false && params.fullreport == true)
            """
              unset DISPLAY
              qualimap bamqc -bam ${bam_file} -outdir ${sample_id}/ -outfile ${sample_id}_report.pdf -c
              cp ./${sample_id}/* .
            """

}


workflow {
    index_ref(refseq_ch)
    refindex = refseq_ch.combine(index_ref.out).map { it -> tuple(it.simpleName, it) }

    if (params.inclUnpRds == true) {
        UNP_reads_idx = unp_reads.combine(refindex)
	bwa_mem2_UNP(UNP_reads_idx)
        samtools_bam_srt_idx_UNP(bwa_mem2_UNP.out)
    }

    if (params.inclMrgRds == true) {
        MRG_reads_idx = mrg_reads.combine(refindex)
        bwa_mem2_MRG(MRG_reads_idx)
        samtools_bam_srt_idx_MRG(bwa_mem2_MRG.out)

    }

    if (params.inclUnpRds == true && params.inclMrgRds == true) {
        samtools_merge(samtools_bam_srt_idx_UNP.out[0].mix(samtools_bam_srt_idx_MRG.out[0]).map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    } else if (params.inclUnpRds == true && params.inclMrgRds == false) {
        samtools_merge(samtools_bam_srt_idx_UNP.out[0].map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    } else if (params.inclUnpRds == false && params.inclMrgRds == true) {
        samtools_merge(samtools_bam_srt_idx_MRG.out[0].map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    }

    qualimap_bamqc(samtools_merge.out[0])
}
