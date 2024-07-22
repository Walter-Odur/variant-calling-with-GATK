# **GATK PIPELINE**
:::info
## :green_book: **A brief Description of the pipeline**
:::

> This pipeline is developed to detect [**somatic mutations**](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/somatic-mutation) in tumor samples. It takes in as input <font color="blue">*tumor*</font> samples and their <font color="blue">*matched normal*</font> samples and outputs <font color="blue">*VCF*</font> files containing somatic variants.
> It provides step-wise procedures for <font color="#1936C9">**GATK**</font> analysis starting from <font color="blue">Quality control, Alignment, base quality score recalibration, variant calling and annotation.</font> The pipeline leverages several essential bioinformatics tools and follows a structured workflow as outlined below.

:::info
## :inbox_tray: **Prerequisites**
:::
Ensure you have the following directories, files, and tools before running the script:
> #### Directories & Files
>>1. reads/
>>> - This directory should contains all the NORMAL and TUMOR samples to be analyzed.
>>2. ref/
>>> - This directory should contain all reference files required for the analysis, including:
>>>>>(i). Human reference genome assembly build 38 (Homo_sapiens_assembly38.fasta).
>>>>>
>>>>>(ii). Known-sites files 
>>>>>
>>>>>(iii).TargetedRegions/interval file
>>>>>

> All reference files and their coresponding indexes are found [**Here**](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)
>> #### Example command to download the references:

```bash=
wget -P ref/ <download_link_for_ref_file>
```
:::warning 
> <font color="red">Ensure all reference files are **unzipped**.</font>
:::

> #### Conda Environment
>> Create and activate new Conda environment to install the necessary tools.

> #### Tools
>> Ensure the following tools are installed:
>> 
>> *Click on the <font color="blue">links</font> below to download <font color="blue">gatk</font> and <font color="blue">snpEff</font>*
>>>- [GATK v4.2.2.0](https://github.com/broadinstitute/gatk/releases?page=2)
>>>- [snpEff 4.1](https://pcingola.github.io/SnpEff/snpeff/running/)

>> *Use mamba/conda to install the following:*
>>>- Samtools 1.2.0
>>>- FastQC v0.12.1
>>>- Fastp 0.23.4
>>>- BWA 0.7.18

```bash=+
# using mamba
mamba install \
    -c bioconda \
    -c conda-forge samtools fastqc fastp bwa
    
# using conda
conda install \
    -c bioconda \
    -c conda-forge samtools fastqc fastp bwa
```

:::info
## :computer: Analysis Steps
:::
#### Initial Setup
```bash=+
#!/bin/bash

# You can time your analysis here
SECONDS=0

# Extract sample IDs
ls ./reads/*fastq | sed 's/_[^_]*.fastq$//g' | sort -u | xargs -i{} basename {} > ./samples.txt || exit 1

# Define variable for sample IDs and reference genome
IDs=$(cat ./samples.txt)
ref=ref/Homo_sapiens_assembly38.fasta

# Create directories for downstream analyses
mkdir -p fastQC/ trim/ bams/ gatk/ recal/ pileups/ mut2/ annotations/ tempDir/
```
:::info
### **Step :one::** <font color="green">Quality Control</font> 
:::

#### Tools Involved:
>>- <font color="red">fastqc
>>- fastp
</font>

#### <font color="blue">(a): Run FastQC</font>
```bash=+
# Check if FastQC is already done then don't run: else run FastQC

if ls fastQC/*fastqc* 1> /dev/null 2>&1; then
    echo "FastQC already done 👌👌👌👌"
else
    echo "Performing quality checks 😊😊😊😊"
    for id in ${IDs[@]}; do
        fastqc \
            ./reads/${id}*.fastq \
            -o ./fastQC || exit 1
    done
    echo "QC Done 👌👌👌👌"
fi
```
#### <font color="blue">(b): Trim Poor Quality Reads</font>
```bash=+
# Check if trimming is already done then don't trim: else trim

if ls ./trim/*trim.fastq 1> /dev/null 2>&1; then
    echo "Trimming already done 👌👌👌👌"
else
    echo "Trimming poor quality reads 😊😊😊😊"
    for id in ${IDs[@]}; do
        fastp \
            -i ./reads/${id}*1.fastq \
            -I ./reads/${id}*2.fastq \
            -o ./trim/${id}_1_trim.fastq \
            -O ./trim/${id}_2_trim.fastq || exit 1
    done
    echo "Trimming Done 👌👌👌👌"
fi
```
:::info
### **Step :two::** <font color="green">Alignment To The Reference And Preparing Bam Files For Variant calling</font>
:::

#### Tools Involved:
>>- <font color="red">BWA
>>- Samtools
>>- GATK
</font>

#### <font color="blue">(a): Align Reads To The Reference</font>
```bash=+
# Perform alignment if not yet done
# The outputs of this script block are indexed bam files

if ls ./bams/*.bam 1> /dev/null 2>&1; then
    echo "Alignment already done 👌👌👌👌"
else
    echo "Performing alignment & generating alignment metrics 😊😊😊😊"
    for id in ${IDs[@]}; do
        bwa mem \
            -t 16 -aM ${ref} \
            -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA" \
            ./trim/${id}_1_trim.fastq \
            ./trim/${id}_2_trim.fastq |\
            samtools view -bS - > ./bams/${id}.bam || exit 1
        samtools sort ./bams/${id}.bam -o ./bams/${id}_sorted.bam || exit 1
        samtools index ./bams/${id}_sorted.bam || exit 1
        samtools flagstat ./bams/${id}_sorted.bam > ./bams/${id}.stats.txt || exit 1
    done
    echo "Alignment Done 👌👌👌👌"
fi
```
#### <font color="blue">(b): Mark Duplicates</font>
```bash=+
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
