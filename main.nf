// Running fastqc
process fastqc {

        publishDir("${params.outdir}/fqc_report", mode: 'copy')
        conda "$projectDir/conda_envs/env-preprocessing"

        input:
        tuple val(sample_id), path(fastq)

        output:
        tuple val(sample_id), path("${sample_id}_R1_fastqc.html"), path("${sample_id}_R1_fastqc.zip"), \
        path("${sample_id}_R2_fastqc.html"), path("${sample_id}_R2_fastqc.zip")

        script:
        """
        fastqc ${fastq[0]} ${fastq[1]} -o .
        """

}

// Trimming Reads
process fastpTrim {

        publishDir("${params.outdir}/trimmed", mode: 'copy')
        conda "$projectDir/conda_envs/env-preprocessing"

        input:
        tuple val(sample_id), path(fastq)

        output:
        tuple val(sample_id), path("${sample_id}_R1_fastpTrimmed.fq.gz"), \
        path("${sample_id}_R2_fastpTrimmed.fq.gz"), emit: trimmed_ch

        script:
        """
        fastp \
        -i ${fastq[0]} \
        -I ${fastq[1]} \
        -o ${sample_id}_R1_fastpTrimmed.fq.gz \
        -O ${sample_id}_R2_fastpTrimmed.fq.gz \
        --detect_adapter_for_pe
        """
}


// Indexing with bowtie2
process bowtie2Build {

        publishDir("${params.refdir}", mode: 'copy')
        conda "$projectDir/conda_envs/env-alignment"

        input:
        path ref

        output:
        path '*', emit: bt2idx_ch

        script:
        """
        bowtie2-build ${ref} ${ref.baseName}
        """
}


// Alingning Reads using bowtie2
process bowtie2Aln {

        publishDir("${params.outdir}/bams", mode: 'copy')
        conda "$projectDir/conda_envs/env-alignment"

        input:
        tuple val(sample_id), path(read1), path(read2)
        path ref
        path bt2idx

        output:
        tuple path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam_ch


        script:
        """
        bowtie2 -x ${ref.baseName} -1 ${read1} -2 ${read2} -S /dev/stdout | \
        samtools view -h -b - | samtools sort -o ${sample_id}_sorted.bam -
        samtools index ${sample_id}_sorted.bam
        """
}