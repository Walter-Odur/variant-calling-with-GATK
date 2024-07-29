# <font color="green"> **SOMATIC VARIANT DETECTION IN CANCER GENOME** </font> 
---
## **STEP1:** <font color="violet">PRE-PROCESSING</font>
---
## **STEP2:** <font color="violet">ALIGNMENT TO THE REFERENCE AND PREPARING BAM FILES FOR VARIANT CALLING</font>

---
### Tools Involved:
>>- <font color="red">BWA
>>- Samtools
>>- GATK
</font>

### <font color="blue">2(a): Align Reads To The Reference</font>
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
mkdir -p ./bams
```
```bash=
# Align paired end reads to the reference
for id in $(cat samples.txt); do
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
>Breakdown of the <font color="red">**options**</font> used
:::
>
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
for id in $(cat samples.txt); do
    samtools sort ./bams/${id}.bam -o ./bams/${id}_sorted.bam 
done
```
:::success
:computer: Index the  sorted bam files
>This enables quick retrieval of specific regions of the genome, significantly speeding up data access and enabling efficient querying of aligned sequencing data.
>You can as well generate the alignment summary statistics for the bam files using <mark>*samtools flagstat*</mark> command
:::
```bash=
for id in $(cat samples.txt); do
    samtools index ./bams/${id}_sorted.bam
    samtools flagstat ./bams/${id}_sorted.bam > ./bams/${id}.stats.txt || exit 1
done
```
:::info
>The next step is to mark sequence duplicates. This is done using GATK (Genome Analysis Toolkit)
>GATK is a software suite developed by the Broad Institute for analyzing high-throughput sequencing data. It offers a wide range of utilities that perform a number of tasks such as quality control, data processing, duplicate marking, base recalibration, variant calling among others
:::
### <font color="blue">2(a): Mark Sequence Duplicates </font>
> The <font color="red">*GATK MarkDuplicates*</font> tool identifies and tags duplicate reads (reads originating from a single DNA fragment) in a BAM or SAM file.
> These duplicates arise during sample preparation or from optical sensor errors.
> 
> The tool compares sequences at the 5' positions, uses molecular barcodes if available, and marks duplicates in the SAM flags field with a hexadecimal value of 0x0400 (decimal 1024), producing a new SAM or BAM file with these annotations
> 
:::success
:computer: Run the following commands to mark duplicates
:::
```bash=
# Make a directory to store the output of gatk MarkDuplicates
mkdir -p ./gatk
```
```bash=
for id in $(cat samples.txt); do
    gatk MarkDuplicates \
        --INPUT ./bams/${id}_sorted.bam \
        --OUTPUT ./gatk/${id}_mkdp.bam \
        --METRICS_FILE ./gatk/${id}_metrics.txt \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true
done
```
:::info
>Breakdown of the <font color="red">**options**</font> used
:::
>>- <font color="red">- -REMOVE_DUPLICATES false:</font> Indicates that duplicates should be marked but not removed from the BAM file.
>>If <font color="red"> true</font> do not write duplicates to the output file
instead of writing them with appropriate flags set.
>>- <font color="red">- -INPUT:</font> One or more input SAM or BAM files to analyze. Must be coordinate sorted.
>>- <font color="red">- -OUTPUT:</font> The output file to write marked records to.
>>- <font color="red">- -METRICS_FILE:</font> File to write duplication metrics to. This file contains statistics about the duplicates found in the input BAM.
>>- <font color="red">- -CREATE_INDEX true:</font> Creates an index for the output BAM file, enabling quick access to specific regions.

### <font color="blue">2(c.): Base Quality Score Recalibration (BQSR)</font>
>The goal of this step is to correct for systematic bias that affect the assignment of base
quality scores by the sequencer.
>
>BQSR is done in two steps:
>> The first pass <mark>*(Generation of recalibration table)*</mark> consists of calculating error empirically and finding patterns in how error varies with basecall features over all bases. The relevant observations are written to a recalibration table.
>> The second pass <mark>*(Applyng generated recalibration data)*</mark> consists of applying numerical corrections to each individual basecall based on the patterns identified in the first step (recorded in the recalibration table) and write out the recalibrated data to a new BAM or CRAM file.
:::success
:computer: Generate Recalibration Table
:::
>Tool: <font color="red"> *gatk BaseRecalibrator* </font>
```bash=
for id in $(cat samples.txt); do
    gatk BaseRecalibrator \
        -I ./gatk/${id}_mkdp.bam \
        -R ${ref} \
        --known-sites ./ref/Mills_and_1000G_gold_standard.indels.hg38.vcf \
        --known-sites ./ref/Homo_sapiens_assembly38.known_indels.vcf \
        --known-sites ./ref/1000G_phase1.snps.high_confidence.hg38.vcf \
        -O ./recal/${id}_recal_data.table
done
```
:::success
:computer: Apply Recalibration data to bam file
:::
> Tool: <font color="red"> *gatk ApplyBQSR*</font>
```bash=
for id in $(cat samples.txt); do
    gatk ApplyBQSR \
        -R ${ref} \
        -I ./gatk/${id}_mkdp.bam \
        --bqsr-recal-file ./recal/${id}_recal_data.table \
        -O ./recal/${id}_recal.bam
done
```
:::info
>Breakdown of the <font color="red">**options**</font> used
:::
>>- <font color="red">--bqsr-recal-file:</font> Input recalibration table for BQSR
>>- <font color="red">--known-sites:</font> One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis. You can obtain the files from [**GOOGLE CLOUD**](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false). After downloading, make sure you index the sites files using <font color="red"> *gatk IndexFeatureFile*</font> 
>>- <font color="red">-I:</font> BAM/SAM/CRAM file containing reads
>>- <font color="red">-O:</font> The file to which the output should be written
>>- <font color="red">-R:</font> Reference sequence


## **STEP3:** <font color="violet">VARIANT CALLING AND ANNOTATION</font>
Downloading population variant file and its index (used to filter out potential germline mutations from the final list of somatic mutations). These are available at [link](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) 
```{bash}
wget -P reference https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -P reference https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
```
Unzipping the vcf file downloaded above.
```{bash}
gunzip $dir3/af-only-gnomad.hg38.vcf.gz
```

```{bash}
k5="af-only-gnomad.hg38.vcf.gz"

# Creating variables with tumor and normal sample names respectively
tsample=$(cat tsamples)
nsample=$(cat nsamples)

for i in $tsample
do
  # Assuming names are separated by space in the respective files
  name1=$(echo "$i" | cut -d ' ' -f 1)
  j=$(echo "$nsample" | cut -d ' ' -f 1)
  echo "$name1,$j"
  START="$( date +%s )"
  gatk  Mutect2 -R $dir3/$ref -I $dir5/${name1}_tumor_recal.bam --tumor-sample ${name1} -I $dir7/${j}_normal_recal.bam --normal-sample ${j} -O $dir6/${name1}_tumor_${j}_normal.Mutect2.raw.vcf --germline-resource $dir3/$k5 --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
  END="$( date +%s )"
  echo "==================Variant calling ${i,j} : $[ $END - $START ] seconds===============" >> time.txt
# Shift remaining names in both variables
  tsample=$(echo "$tsample" | cut -d ' ' -sf 2-)
  nsample=$(echo "$nsample" | cut -d ' ' -sf 2-)
done
```
:::success
### :feet: Step 3b: Variant filtering
:::



```{bash}
for i in $tsample
do
  # Assuming names are separated by space
  name1=$(echo "$i" | cut -d ' ' -f 1)
  j=$(echo "$nsample" | cut -d ' ' -f 1)
  echo "$name1,$j"
  if [ -e $dir5/${name1}_tumor_recal.bam ]
	then
		echo "Found file for ${i}_tumrecal_sorted.bam"
	       	echo "Getting pileup information for the samples at sites of known mutations"
		START="$( date +%s )"
	       	gatk GetPileupSummaries -I $dir5/${name1}_tumor_recal.bam -O $dir6/${name1}_tumor.pileups.table --variant $dir3/$k5 --intervals Illumina_Exome_TargetedRegions_v1.2.hg38.bed --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
		END="$( date +%s )"
		echo "==================Generating pileup summaries for ${name1} : $[ $END - $START ] seconds===============" >> time.txt
	else
		echo "No corresponding file found for ${i}_tumrecal_sorted.bam"
	fi
	if [ -e $dir7/${j}_normal_recal.bam ]
	then
		echo "Found file for ${i}_normrecal_sorted.bam"
		START="$( date +%s )"
	       	gatk GetPileupSummaries -I $dir7/${j}_normal_recal.bam -O $dir6/${j}_normal.pileups.table --variant $dir3/$k5 --intervals Illumina_Exome_TargetedRegions_v1.2.hg38.bed --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
		END="$( date +%s )"
		echo "==================Generating pileup summaries for ${j} : $[ $END - $START ] seconds===============" >> time.txt
	else
		echo "No corresponding file found for ${i}_normrecal_sorted.bam"
	fi
	if [ -e $dir6/${name1}_tumor.pileups.table ] || [ -e $dir6/${j}_normal.pileups.table ]
	then
		echo "Found files for ${name1}_tumor.pileups.table and ${j}_normal.pileups.table"
	       	echo "Estimation of contamination"
		START="$( date +%s )"
   	   	gatk CalculateContamination -I $dir6/${name1}_tumor.pileups.table -O $dir6/${name1}_Tumor_${j}_Normal.contamination.table --matched-normal $dir6/${j}_normal.pileups.table --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
		END="$( date +%s )"
		echo "==================Calculating contaminations for ${name1,j} : $[ $END - $START ] seconds===============" >> time.txt
	       	echo "Application of a first filter to variant calls"
		START="$( date +%s )"
	       	gatk FilterMutectCalls --variant $dir6/${name1}_tumor_${j}_normal.Mutect2.raw.vcf -O $dir6/${name1}_Tumor_${j}_Normal.Mutect2.oncefiltered.vcf --contamination-table $dir6/${name1}_Tumor_${j}_Normal.contamination.table --reference $dir3/$ref --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
		END="$( date +%s )"
		echo "==================Applying first filtering parameters ${name1,j} : $[ $END - $START ] seconds===============" >> time.txt
#		echo "Applying the second pass filter"
#		START="$( date +%s )"
#	       	gatk FilterByOrientationBias --variant $dir6/${name1}_Tumor_${j}_Normal.Mutect2.oncefiltered.vcf  -O $dir6/${name1}_Tumor_${j}_Normal.Mutect2.twicefiltered.vcf
#		END="$( date +%s )"
#		echo "==================Applying second filtering parameters ${name1,j} : $[ $END - $START ] seconds===============" >> time.txt
     		echo "Removing variants that failed to pass the filter check"
		START="$( date +%s )"
	       	awk '/^#/ {print $0; next} $7=="PASS" {print $0}' $dir6/${name1}_Tumor_${j}_Normal.Mutect2.oncefiltered.vcf > $dir6/${name1}_T_${j}_N.vcf
		END="$( date +%s )"
		echo "==================Removing variants that failed to pass the filters for ${name1,j} : $[ $END - $START ] seconds===============" >> time.txt
	else
		echo "No corresponding file found for SRR25434460_tumor.pileups.table and SRR25434461_normal.pileups.table"
	fi
# Shift remaining names in both variables
  tsample=$(echo "$tsample" | cut -d ' ' -sf 2-)
  nsample=$(echo "$nsample" | cut -d ' ' -sf 2-)
done
```
:::success
### :feet: Step 3c: variant annotation
:::


```{bash}
for i in $tsample
do
  # Assuming names are separated by space
  name1=$(echo "$i" | cut -d ' ' -f 1)
  j=$(echo "$nsample" | cut -d ' ' -f 1)
  echo "$name1,$j"
  echo "Variant Annotation of $name1 and $j using snpeff"
  START="$( date +%s )"
  snpEff download -v hg38
  END="$( date +%s )"
  echo "================== Downloading the hg38 reference database : $[ $END - $START ] seconds===============" >> time.txt
  START="$( date +%s )"
  snpEff hg38 $dir6/${name1}_T_${j}_N.vcf > $dir6/${name1}_T_${j}_N_annsnpEff.vcf
  snpEff hg38 -s ${name1}_T_${j}_N_summary.html -stats ${name1}_T_${j}_N_stats.txt $dir6/${name1}_T_${j}_N.vcf > $dir6/${name1}_T_${j}_N_annsnpEff.vcf
  END="$( date +%s )"
  echo "==================Applying annotations for ${name1,j} : $[ $END - $START ] seconds===============" >> time.txt
 # Shift remaining names in both variables
  tsample=$(echo "$tsample" | cut -d ' ' -sf 2-)
  nsample=$(echo "$nsample" | cut -d ' ' -sf 2-)
done
```


1. Create the repository on GitHub(website); if it already exists, skip it.
2. Go to terminals
	1. git clone https://github.com/repository.git
	2. cd repository
	3. git status
	4. git add .
	5. git status
	6. git commit -m "[what are you do]"
	7. git push
	8. Enter: the username
	9. Enter: token (ps: password had been canceled)
        - Settings
        - Developer settings
        - Personal access tokens