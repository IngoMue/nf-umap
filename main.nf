refseq_ch = Channel.fromPath(params.refdir)


reads_ch = Channel.fromPath(params.reads)
        .map{ [ it.name.tokenize("_")[0..-2].join("_"), it] }
        .groupTuple()

process index_ref {
    publishDir "$params.outdir/indexref/", pattern: '*.ref', mode: 'copy'
    tag "$ref_file.baseName"
    input:
      file(ref_file) from refseq_ch
    output:
      file 'index.ref' into idx_ref
    script:
    """
      bwa-mem2 index $ref_file > index.ref
    """
}



process bwa_mem2 {
    publishDir "$params.outdir/sam/", pattern: '*.sam', mode: 'copy'
    tag "$sample_id"
    input:
      each file(ref_file) from idx_ref
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


process sort_bam {
    publishDir "$params.outdir/sorted/", pattern: '*_sorted.bam', mode: 'copy'
    tag "$sample_id"
    input:
      tuple val(sample_id), file(sample_file) from bams_ch
    output:
      tuple val(sample_id), file("${sample_id}_sorted.bam") into sorted_ch
    script:
    """
      samtools sort -o ${sample_id}_sorted.bam $sample_file
    """
}

process index_bam {
    publishDir "$params.outdir/sorted/", pattern: '*.bai', mode: 'copy'
    tag "$sample_id"
    input:
      tuple val(sample_id), file(sample_file) from sorted_ch
    output:
      file("${sample_id}.bai")
    script:
    """
      samtools index $sample_file
    """
}

