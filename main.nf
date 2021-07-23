#!/usr/bin/env nextflow

log.info """\
         nf - mapping - A simple nextflow mapping workflow     
         =================================================
         refdir        : ${params.refdir}
         refname       : ${params.refname}
         reads         : ${params.reads}
         outdir        : ${params.outdir}
         """
         .stripIndent()

refdir_ch = Channel.fromPath(params.refdir)

refseq_ch = Channel.fromPath("${params.refdir}/${params.refname}")

reads_ch = Channel.fromPath(params.reads)
        .map{ [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()

process index_ref {
    label 'High_RAM'

    publishDir "$params.refdir", pattern: '*.fai', mode: 'copy'
    tag "$ref_file.baseName"

    input:
      file(ref_file) from refseq_ch

    output:
      file 'index.ref'

    script:
    """
      bwa-mem2 index $ref_file
    """
}

refseq_ch = Channel.fromPath("${params.refdir}/${params.refname}")

process bwa_mem2 {
    label 'High_RAM'

    publishDir "$params.outdir/sam/", pattern: '*.sam', mode: 'copy'
    tag "$sample_id"

    input:
      each file(ref_file) from refseq_ch
      tuple val(sample_id), file(sample_file) from reads_ch

    output:
      tuple val(sample_id), file("${sample_id}.sam") into sams_ch

    script:
    """
      bwa-mem2 mem $ref_file $sample_file > ${sample_id}.sam
    """
}


process convert_to_bam {
    publishDir "$params.outdir/bam/", pattern: '*.bam', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sample_file) from sams_ch

    output:
      tuple val(sample_id), file("${sample_id}.bam") into bams_ch

    script:
    """
      samtools view -bS $sample_file > ${sample_id}.bam
    """
}


process sort_index_bam {
    publishDir "$params.outdir/sorted/", pattern: '*_sorted.{bam,bam.bai}', mode: 'copy'
    tag "$sample_id"

    input:
      tuple val(sample_id), file(sample_file) from bams_ch

    output:
      file("${sample_id}_sorted.bam")
      file("${sample_id}_sorted.bam.bai")

    script:
    """
      samtools sort -o ${sample_id}_sorted.bam $sample_file

      samtools index ${sample_id}_sorted.bam
    """
}
