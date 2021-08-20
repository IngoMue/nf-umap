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
        .map { [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()

m_reads = Channel.fromPath(params.merged_reads)
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


process bwa_mem2_PE {
    label 'High_RAM'

    tag "$sample_id"

    input:
      tuple val(sample_id), file(sample_file), val(ref_ID), file(index)

    output:
      tuple val(sample_id), file("${sample_id}_PE.sam")

    script:
    """
      bwa-mem2 mem $params.refprefix $sample_file > ${sample_id}_PE.sam
    """
}

process bwa_mem2_SE {
    label 'High_RAM'

    tag "$mrg_id"

    input:
      tuple val(mrg_id), file(merged_file), val(ref_ID), file(index)

    output:
      tuple val(mrg_id), file("${mrg_id}_SE.sam")

    script:
    """
      bwa-mem2 mem $params.refprefix $merged_file > ${mrg_id}_SE.sam
    """
}


process samtools_bam_srt_idx_PE {
    publishDir "${params.outdir}/02_bams/PairedEnd", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_PE_sorted.bam")
      file("${sample_id}_PE_sorted.bam.bai")

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_PE.bam

      samtools sort -o ${sample_id}_PE_sorted.bam ${sample_id}_PE.bam

      samtools index ${sample_id}_PE_sorted.bam
    """
}

process samtools_bam_srt_idx_SE {
    publishDir "${params.outdir}/02_bams/SingleEnd", pattern: '*_sorted{.bam,.bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sam_file)

    output:
      file("${sample_id}_SE_sorted.bam")
      file("${sample_id}_SE_sorted.bam.bai")

    script:
    """
      samtools view -bS $sam_file > ${sample_id}_SE.bam

      samtools sort -o ${sample_id}_SE_sorted.bam ${sample_id}_SE.bam

      samtools index ${sample_id}_SE_sorted.bam
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

process samtools_coverage {
    publishDir "${params.outdir}/04_mapping_stats/", pattern: '*{.stats,.hist}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(bam_file)

    output:
      file '*'

    script:
    """
      samtools coverage -o ${sample_id}.stats ${bam_file}

      samtools coverage -m -o ${sample_id}.hist ${bam_file}
    """
}


workflow {
    index_ref(refseq_ch)
    refindex = refseq_ch.combine(index_ref.out).map { it -> tuple(it.simpleName, it) }
    PE_reads_idx = p_reads.combine(refindex)
    SE_reads_idx = m_reads.combine(refindex)
    bwa_mem2_PE(PE_reads_idx)
    bwa_mem2_SE(SE_reads_idx)
    samtools_bam_srt_idx_PE(bwa_mem2_PE.out)
    samtools_bam_srt_idx_SE(bwa_mem2_SE.out)
    samtools_merge(samtools_bam_srt_idx_PE.out[0].mix(samtools_bam_srt_idx_SE.out[0]).map { [ it.name.split(/_L0\d+/)[0], it] }.groupTuple())
    samtools_coverage(samtools_merge.out[0])
}
