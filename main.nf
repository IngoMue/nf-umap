#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\

         NF-MAPPING - A nextflow mapping workflow     
         ================================================
         |refdir       : ${params.refdir}
         |refname      : ${params.refname}
         |refprefix    : ${params.refprefix}
         |readpairs    : ${params.read_pairs}
         |mergedreads  : ${params.merged_reads}
         |outdir       : ${params.outdir}
         ================================================
         """
         .stripIndent()

refdir_ch = Channel.fromPath(params.refdir)

refseq_ch = Channel.fromPath("${params.refdir}/${params.refname}")

p_reads = Channel.fromPath(params.read_pairs)
        .map { [ it.name.tokenize("_")[0..-3].join("_"), it] }
        .groupTuple()

m_reads = Channel.fromPath(params.merged_reads)
        .map { [ it.name.tokenize("_")[0..-3].join("_"), it] }
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


process bwa_mem2_pairs {
    label 'High_RAM'

    publishDir "${params.outdir}/02_samfiles/pairs/", pattern: '*.sam', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sample_file), val(ref_ID), file(index)

    output:
      tuple val(sample_id), file("${sample_id}_P.sam")

    script:
    """
      bwa-mem2 mem $params.refprefix $sample_file > ${sample_id}_P.sam
    """
}

process bwa_mem2_merged {
    label 'High_RAM'

    publishDir "$params.outdir/02_samfiles/merged", pattern: '*.sam', mode: 'copy'
    tag "$mrg_id"

    input:
      tuple val(mrg_id), file(merged_file), val(ref_ID), file(index)

    output:
      tuple val(mrg_id), file("${mrg_id}_M.sam")

    script:
    """
      bwa-mem2 mem $params.refprefix $merged_file > ${mrg_id}_M.sam
    """
}


process samtools_bam_srt_idx_p {
    publishDir "${params.outdir}/03_bams/", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_P_sorted.bam")
      file("${sample_id}_P_sorted.bam.bai")

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_P.bam

      samtools sort -o ${sample_id}_P_sorted.bam ${sample_id}_P.bam

      samtools index ${sample_id}_P_sorted.bam
    """
}

process samtools_bam_srt_idx_m {
    publishDir "${params.outdir}/03_bams/", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_M_sorted.bam")
      file("${sample_id}_M_sorted.bam.bai")

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_M.bam

      samtools sort -o ${sample_id}_M_sorted.bam ${sample_id}_M.bam

      samtools index ${sample_id}_M_sorted.bam
    """
}


process samtools_merge {
    publishDir "${params.outdir}/04_merged/", pattern: '*{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(bams)

    output:
      file '*'

    script:
    """
      samtools merge -o ${sample_id}.bam $bams

      samtools index ${sample_id}.bam
    """
}

workflow {
    index_ref(refseq_ch)
    refindex = refseq_ch.combine(index_ref.out).map { it -> tuple(it.simpleName, it) }
    p_reads_idx = p_reads.combine(refindex)
    m_reads_idx = m_reads.combine(refindex)
    bwa_mem2_pairs(p_reads_idx)
    bwa_mem2_merged(m_reads_idx)
    samtools_bam_srt_idx_p(bwa_mem2_pairs.out)
    samtools_bam_srt_idx_m(bwa_mem2_merged.out)
    samtools_merge(samtools_bam_srt_idx_p.out[0].mix(samtools_bam_srt_idx_m.out[0]).map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
}
