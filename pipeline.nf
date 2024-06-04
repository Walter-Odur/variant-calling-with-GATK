nextflow.enable.dsl=2

/* GATK somatic variant calling pipeline
 * Usage: nextflow run /mnt/c/Users/Walterrickman/Desktop/cancer/nextflow/GatkVariantCallingPipeline.nf
 * Author: Cancer Genomics working Group
 * African Centers Of Excellence in Bioinformatics and Data Intensive Sciences (ACE-Uganda)
 * Infectious Diseases Institute (IDI)
 * Makerere University, Kampala Uganda
 */

// Define some parameters here
params.fastq    = "/mnt/c/Users/Walterrickman/Desktop/cancer/nextflow/TestDir/reads/*_R{1,2}.fastq.gz"
params.outdir   = "/mnt/c/Users/Walterrickman/Desktop/cancer/nextflow/TestDir/results"
params.ref      = "/mnt/c/Users/Walterrickman/Desktop/cancer/nextflow/TestDir/reference/Agy99.fasta"
params.refdir   = "/mnt/c/Users/Walterrickman/Desktop/cancer/nextflow/TestDir/reference/"

include {fastqc}        from './bowtie2Main.nf'
include {fastpTrim}     from './bowtie2Main.nf'
include {bowtie2Build}  from './bowtie2Main.nf'
include {bowtie2Aln}    from './bowtie2Main.nf'

// Running the Workflow

workflow {
        fastq_ch        = Channel.fromFilePairs(params.fastq)
        fastqc(fastq_ch).view()
        fastpTrim(fastq_ch)

        fastpTrim.out.view()

        bt2build_ch     = Channel.fromPath(params.ref)
        bowtie2Build(bt2build_ch)
        bowtie2Build.out.view()

        bowtie2Aln(fastpTrim.out.trimmed_ch, bt2build_ch, bowtie2Build.out.bt2idx_ch)

}