# <font color="green"> **SOMATIC VARIANT DETECTION IN CANCER** </font> 
---
## **STEP2:** <font color="green">Alignment To The Reference And Preparing Bam Files For Variant            calling</font>

---
### Tools Involved:
>>- <font color="red">BWA
>>- Samtools
>>- GATK
</font>

### <font color="blue">(a): Align Reads To The Reference</font>
---
>This process involves mapping sequencing reads to their corresponding positions on a reference genome purposely to determine where they originate from within the genome
>
>A number of tools have been developed to perform reads alignment to the reference genome. 
>Some of these include:
>>- BWA (Burrows-Wheeler Aligner)
>>- Bowtie
>>- HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts)
>>- STAR (Spliced Transcripts Alignment to a Reference)
>>- Minimap2
>>- Novoalign
>>- Subread

>Here we shall use <font color="blue">BWA</font> due to its speed, accuracy and efficiency in handling short read sequencing errors and small indels.
>
><font color="blue">**BWA**</font> is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
>It consists of three algorithms: <font color="blue">*BWA-backtrack*</font>, <font color="blue">*BWA-SW*</font> and <font color="blue">*BWA-MEM*</font>. 
>The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
>BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but <mark> <font color="blue">***BWA-MEM***</font>, which is the latest, is generally recommended for highquality queries as it is faster and more accurate</mark>.
> BWA-MEM also has better performance than BWAbacktrack for 70-100bp Illumina reads.
:::success
We shall therefore use **BWA MEM** to perform alignment
:::
>Before performing the alignment, the reference <font color="red">**must**</font> first be <font color="red">**Indexed**</font>.
>- This enables fast, efficient, and accurate alignment of sequencing reads by creating a searchable data structure that optimizes memory usage and computational load.
>This is done using <font color="blue">**bwa index**</font> command as below.
```bash=
bwa index Homo_sapiens_assembly38.fasta
```
>However, you can download the human reference genome build 38 <mark>(Homo_sapiens_assembly38.fasta)</mark> and it's fasta index <mark>(Homo_sapiens_assembly38.fasta.fai)</mark> from <font color="green">**google cloud**</font> following this [**LINK**](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)
>
>Alternatively, you can download them by running the following ***commands:***
:::success
Download the Human Refernce Genome Build 38 from google cloud
:::
```bash=
wget \
    -P ./reads \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
```
:::success
Download the corresponding fasta index file
:::
```bash=
wget \
    -P ./reads \
    https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
```
>The <font color="red">***-P***</font> option in the <font color="red">***wget***</font> command is used to specify the directory to which the file should be downloaded, in this we are downloading the reference into the directory called reads which must be present in the our working directory before running the wget command. When the option is not specified, the reference is downloaded in the current working directory.

:::success
:green_book: Rember that in **STEP ONE** we did preprocessing where we trimmed reads to remove poor quality ones as well as sequencing adapters. We saved our trimmed reads in a directory called <mark>trim</mark>
:desktop_computer: Here you are to align those trimmed reads to the reference genome which you have just obtained.
::: 
:::warning
:no_entry_sign: It is better and much easier to work with <font color="purple">**sample IDs/basename**</font> other than the full <font color="purple">**sample name**</font> 
:::
>Run the following commands to extract sample IDs
```bash=
ls ./reads/*fastq |sed 's/_[^_]*.fastq$//g' |sort -u |xargs -i{} basename {} > ./samples.txt
```
:::success
Check whether sample IDs have been extracted correctly.
:::
```bash=
cat samples.txt
```
>The output should look like this, without <font color="red">**(_*fastq)**</font> extension
```bash=
tumor
normal
```
:::success
Run the following commands to perform the <font color="blue">**alignment**</font>.
:::

```bash=
# ensure you have a directory called bams to store alignment results
mkdir ./bams
```
```bash=
# Align paired end reads to the reference
for id in $(cat samples.txt]; do
    bwa mem \
        -t 4 -aM ./ref/Homo_sapiens_assembly38.fasta \
        -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA" \
        ./trim/${id}_1_trim.fastq \
        ./trim/${id}_2_trim.fastq |\
        samtools view -bS - > ./bams/${id}.bam
done
```
:::info
Th <mark>for loop</mark> enables you to perform alignment on multiple samples as above
:::
>Breakdown of the <font color="red">**options**</font> used
>>- <font color="red">-t:</font> This specifies the number of threads(CPU) to use, <font color="green">*4 CPUs*</font> in this case. This speeds up the process. You can adjust this value depending on your need and the capability of your machine
>>- <font color="red">-a:</font> Outputs all alignments for SE or unpaired PE reads, not just the best one
>>- <font color="red">-M:</font> Marks shorter split hits as secondary alignment, for <font color="green">*Picard compatibility*</font>
>>- <font color="red">-R:</font> Adds read group information to the output BAM file.
>>- <font color="red">-samtools view:</font> This command converts SAM format which is the alignment output into BAM format.
>>- <font color="red">-b:</font> specifies that the output should be in BAM format
>>- <font color="red">-S:</font> Specifies that the Input to samtools is in SAM format, which is the output from bwa mem
>>- <font color="red">**-**:</font>  This (-) at the end of bS option specifies that the input for samtools is the output of the previous command (bwa mem, in this case)

:::info
>The next step is to <font color="orange">index the output bam</font> However, this requires <font color="orange">sorted bam files</font>; 
:::
:::success
:computer: Sort bam files
><font color="orange">Sorting</font> is to organize the aligned reads either by their genomic coordinates or by read names, which is essential for downstream analyses such as duplicate marking, variant calling, and efficient data retrieval.
:::
```bash=
for id in $(cat samples.txt]; do
    samtools sort ./bams/${id}.bam -o ./bams/${id}_sorted.bam 
done
```
:::success
:computer: Index the  sorted bam files
>This enables quick retrieval of specific regions of the genome, significantly speeding up data access and enabling efficient querying of aligned sequencing data.
>You can as well generate the alignment summary statistics for the bam files using <mark>*samtools flagstat*</mark> command
:::
```bash=
for id in $(cat samples.txt]; do
    samtools index ./bams/${id}_sorted.bam
    samtools flagstat ./bams/${id}_sorted.bam > ./bams/${id}.stats.txt || exit 1
done
```
:::info
>The next step is to mark sequence duplicates. This is done using GATK (Genome Analysis Toolkit)
>GATK is a software suite developed by the Broad Institute for analyzing high-throughput sequencing data. It offers a wide range of utilities that perform a number of tasks such as quality control, data processing, duplicate marking, base recalibration, variant calling among others
:::
   
   
   
   
   
   
```
#### <font color="blue">(b): Mark Duplicates</font>
```bash=
# Marh Duplicates if not yet done

if ls ./gatk/*_mkdp.bam 1> /dev/null 2>&1; then
    echo "Duplicates already marked 👌👌👌👌"
else
    echo "Marking duplicates in aligned samples 😊😊😊😊"
    for id in ${IDs[@]}; do
        gatk MarkDuplicates \
            --INPUT ./bams/${id}_sorted.bam \
            --OUTPUT ./gatk/${id}_mkdp.bam \
            --METRICS_FILE ./gatk/${id}_metrics.txt \
            --REMOVE_DUPLICATES false \
            --CREATE_INDEX true || exit 1
        samtools flagstat ./gatk/${id}_mkdp.bam > ./gatk/${id}_mkdp.bam.flagstat || exit 1
    done
    echo "Marking Duplicates Done 👌👌👌👌"
fi
```
#### <font color="blue">(c.): Base Quality Score Recalibration (BQSR)</font>
```bash=+
if ls ./recal/*_recal.bam 1> /dev/null 2>&1; then
    echo "BQSR already done 👌👌👌👌"
else
    echo "Performing BQSR 😊😊😊😊"
    for id in ${IDs[@]}; do
    
        # Generate recalibration table
        gatk BaseRecalibrator \
            -I ./gatk/${id}_mkdp.bam \
            -R ${ref} \
            --known-sites ./ref/Mills_and_1000G_gold_standard.indels.hg38.vcf \
            --known-sites ./ref/Homo_sapiens_assembly38.known_indels.vcf \
            --known-sites ./ref/1000G_phase1.snps.high_confidence.hg38.vcf \
            -O ./recal/${id}_recal_data.table || exit 1

        # Apply the recalibration data
        gatk ApplyBQSR \
            -R ${ref} \
            -I ./gatk/${id}_mkdp.bam \
            --bqsr-recal-file ./recal/${id}_recal_data.table \
            -O ./recal/${id}_recal.bam || exit 1
    done
    echo "Base Recalibration Done 👌👌👌👌"
fi

# This step outputs recalibrated bam file with the extension _recal.bam
```

:::success
### **:congratulations: <font color="purple">Great Job!! We Are Now Heading To Variant Calling!!**</font>
:::
> Step Three (3) will be our final step for somatic variant calling, and we have a <font color="blue">**filtered vcf file**</font> containing somatic sequence variants.

> Let's proceed to the final analysis step.

:::info
### **Step :three::** <font color="green">Somatic Variant Calling and Annotation</font> 
:::

#### Tools Involved:
>>- <font color="red">Mutect2
>>- GetPileupSummaries
>>- CalculateContamination
>>- FilterMutectCalls
>>- SnpEff
</font>

#### <font color="purple">**Define Tumor-Normal Sample Pair**</font>
```bash=+
declare -A samples=(
    ["SRR25434464"]="SRR25434464_recal.bam"
    ["SRR25434465"]="SRR25434465_recal.bam"
    ["SRR25434466"]="SRR25434466_recal.bam"
    ["SRR25434467"]="SRR25434467_recal.bam"
)
```
> We shall be using functions to simplify and fasten the process

#### <font color="blue">(a): Variant Calling</font>
>- Create a function to call variants

```bash=+
call_variants() {
        local tumor_sample=$1
        local normal_sample=$2
        java \
                -Dsamjdk.use_async_io_read_samtools=false \
                -Dsamjdk.use_async_io_write_samtools=true \
                -Dsamjdk.use_async_io_write_tribble=false \
                -Dsamjdk.compression_level=2 \
                -Xmx16g \
                -Djava.io.tmpdir=tempDir \
                -jar /opt/ohpc/admin/spack/0.17.0/opt/spack/linux-rocky8-sandybridge/gcc-8.5.0/gatk-4.2.2.0-7fb72tl4fvlec6vn7yw46o4zdq2cspnz/bin/gatk-package-4.2.2.0-local.jar Mutect2 \
                -I ./recal/${tumor_sample}_recal.bam \
                --tumor-sample ${tumor_sample} \
                -I ./recal/${normal_sample}_recal.bam \
                --normal-sample ${normal_sample} \
                -O ./mut2/${tumor_sample}_vs_${normal_sample}.Mut2.raw.vcf \
                --germline-resource ./ref/af-only-gnomad.hg38.vcf \
                --tmp-dir ./tempDir
}
```
#### <font color="blue">(b): Variant Filtering </font>
>- Create a function to filter variants

```bash=+
filter_variants() {
        local tumor_sample=$1
        local normal_sample=$2
        local mut2_vcf=mut2/${tumor_sample}_vs_${normal_sample}.Mut2.raw.vcf

        if ls ./mut2/${tumor_sample}_vs_${normal_sample}.Mut2.Filtered.vcf 1> /dev/null 2>&1; then
                echo "Variant filtering already done👌👌👌👌"
        else
                echo "performing Variant filtering 😊😊😊😊"
                for sample in ${tumor_sample} ${normal_sample}; do
                        if [ -e ./recal/${sample}_recal.bam ]; then
                                echo "Found file for ${sample}_recal.bam"
                                echo "Getting pileup information for the samples at sites of known mutations"
                #---------
                # [i]:  Get pileup summaries
                #---------
                                gatk GetPileupSummaries \
                                        -I ./recal/${sample}_recal.bam \
                                        -O ./pileups/${sample}.pileups.table \
                                        --variant ./ref/af-only-gnomad.hg38.vcf \
                                        --intervals ./ref/Illumina_Exome_TargetedRegions_v1.2.hg38.bed \
                                        --TMP_DIR tempDir \
                                        --java-options '-Xmx16g -Djava.io.tmpdir=tempDir'
                        else
                                echo "No corresponding file found for ${sample}_recal.bam"
                        fi
                done

                if [ -e pileups/${tumor_sample}.pileups.table ] && [ -e pileups/${normal_sample}.pileups.table ]; then
                        echo "Found pileup files for ${tumor_sample} and ${normal_sample}"

                        echo "Estimating contamination"

                #---------
                # [i]:  Calculate contamination
                #---------
                        gatk CalculateContamination \
                                -I ./pileups/${tumor_sample}.pileups.table \
                                -O ./pileups/Tumor_Normal.contamination.table \
                                --matched-normal ./pileups/${normal_sample}.pileups.table \
                                --TMP_DIR ./tempDir \
                                --java-options '-Xmx16g -Djava.io.tmpdir=tempDir'

                        echo "Applying filter to variant calls"

                #---------
                # [ii]: Apply the Filter to variants
                #---------
                        gatk FilterMutectCalls \
                                --variant ${mut2_vcf} \
                                -O ./mut2/${tumor_sample}_vs_${normal_sample}.Mut2.filtered.vcf \
                                --contamination-table ./pileups/Tumor_Normal.contamination.table \
                                --reference ./ref/Homo_sapiens_assembly38.fasta \
                                --TMP_DIR ./tempDir \
                                --java-options '-Xmx16g -Djava.io.tmpdir=tempDir'

                        echo "Removing variants that failed to pass the filter check"
                
                #---------
                # [iii]: Extract Variants that passed the filter
                #---------
                        awk \
                                '/^#/ {print $0; next} $7=="PASS" {print $0}' \
                                ./mut2/${tumor_sample}_vs_${normal_sample}.Mut2.filtered.vcf > \
                                ./mut2/${tumor_sample}_vs_${normal_sample}.Mut2.Filtered.vcf
                else
                        echo "No corresponding pileup files found for ${tumor_sample} and ${normal_sample}"
                fi
                echo "Variant filtering done👌👌👌👌"
        fi
}
```

#### <font color="blue">(c.): Variant Annotation</font>
>- Create a function to annotate variants

```bash=+
annotate_variants() {
        local filtered_vcf=mut2/${tumor_sample}_vs_${normal_sample}.Mut2.Filtered.vcf
        if ls ./annotations/*_annsnpEff.vcf 1> /dev/null 2>&1; then
                echo "variants already annotated👌👌👌👌"
        else
                echo "Annotating variants 😊😊😊😊"
                java -jar \
                        ./snpEff/snpEff.jar \
                        hg38 \
                        ${filtered_vcf} > \
                        ./annotations/$(basename ${filtered_vcf} .vcf)_annsnpEff.vcf
        fi
}
```

:::success
### **:congratulations: <font color="purple">Now Lets Execute Our Functions Here!!**</font>
:::

```bash=+
call_variants "SRR25434464" "SRR25434465" &
call_variants "SRR25434466" "SRR25434467" &
wait

filter_variants "SRR25434464" "SRR25434465" &
filter_variants "SRR25434466" "SRR25434467" &
wait

annotate_variants \
        mut2/SRR25434464_vs_SRR25434465.Mut2.Filtered.vcf &
annotate_variants \
        mut2/SRR25434466_vs_SRR25434467.Mut2.Filtered.vcf &

echo "All processes completed successfully!"

echo "Time taken: $(($duration / 3600))h: $((($duration % 3600) / 60))m: $(($duration % 60))s"
```

:::success
**A brief description of some few commands used*** 
:::
> <font color="red">**||exit 1**</font>
>>- To <font color="blue">**exit the script**</font> once that particular script block has not run successfully.
>>
> <font color="red">**&**</font> on at the ends of the functions:
>>- To allow <font color="blue">**parallel run**</font> (Processing multiple sample pairs simultaneosly)

> <font color="red">**wait**</font> after calling the function:
>>- To <font color="blue">**inhibit**</font> two adjascent processes from running in parallel as the next process must wait for the output of the running process and takes in as its input. 

>**<font color="blue">Congratulations! You have just learnt a Comprehensive pipeline for detecting somatic mutation. All the Best In Your Cancer Genomics Journey!
</font>**
