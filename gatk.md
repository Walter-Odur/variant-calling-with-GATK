![aceLogo](https://hackmd.io/_uploads/Bk0UcfrtC.png)


# <font color="darkblue"> **SOMATIC VARIANT DETECTION IN CANCER GENOME** </font> 
### Background
On 29th November, 2023 members of the ACE Cancer Genomics Working Group agreed on building a custom variant calling pipeline which was to be curated, documented and finally adopted by the members of the group for use. We tested two popular tools (Vascan and GATK) used for variant calling as documented in this [PAPER](https://doi.org/10.1007/978-1-4939-8876-1_21).
The main objective was to build a scableble, automated and well documented variant calling pipeline for our daily use.

<font color="red">This text is red.</font>


## :feet: **STEP1A:** Sample Acquisition
Samples.txt file was created with all the accession numbers of the samples in BioProject [PRJNA851929](https://www.ebi.ac.uk/ena/browser/view/PRJNA851929). The samples of tumor tissue were obtained from needle biopsy while the non tumor tissue were obtained from mononuclear cells isolated by Ficoll-gradient from EDTA-anticoagulated whole blood extracted at diagnosis and preserved at -80C.

```bash
for id in $(cat Samples.txt)
do
    fasterq-dump $sample
done
```
## :feet: **STEP1B:** PRE-PROCESSING
### a). Quality Control using FASTQC

FASTQC is a Java-based quality control tool providing per-base and per-read quality profiling features. It is a program that will read through the raw sequence files (fastq) and perform a quality check on the data. It produces a report that summaries the information. We run a quality assessment on both the normal paired end reads and the tumor paired end reads.

A summary report is generated and it is divided into 10 analysis modules each giving details on the results obtained. An evaluation of the results obtained in each module is denoted by a green tick (for passed), a yellow exclamation mark (for any warnings) and a red cross (for failed). Go to [LINK](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details about fastqc

This step is important for;
* **Identifying Quality Issues:** FastQC provides valuable insights into the quality of sequencing reads, such as low-quality bases, adapter contamination, and nucleotide bias.
* **Targeted Trimming:** Based on the FastQC report, you can tailor your trimming strategy to remove specific regions of low quality or contamination.

```bash
mkdir -p fastqcs
fastqc reads/*.fastq -o fastqcs
```
MultiQC which aggregate results from Fastqc analyses across many samples into a single report can be employed.

MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools. Go to [LINK](https://multiqc.info/) for details about MultiQC
```bash
multiqc fastqcs/.
```

### b). Trimming reads using fastp
Trimming sequence data is often performed after running FastQC to improve read quality. High-quality reads are essential for accurate downstream analyses like alignment, assembly, and variant calling. Trimming can significantly improve the reliability and sensitivity of these processes.

**Common Trimming Methods include:**

* **Quality-based trimming:** Removes bases that fall below a certain quality threshold.
* **Adapter trimming:** Identifies and removes adapter sequences that may be present at the ends of reads.
* **Length-based trimming:** Removes reads that are shorter or longer than a specified length.

**Popular Trimming Tools include:**

* **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**
* **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)**
* **[BBMap](https://sourceforge.net/projects/bbmap/)**
* **[TrimGalore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)**
* **[Fastp](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)**

NB: Click the tool for more details

Fastp was chosen as the tool to be used for trimming. Fastp is an ultra-fast tool to perform quality control, read filtering and base correction with a single scan FASTQ files. It includes most features of FASTQC + Cutadapt + Trimmomatic + AfterQC while running 2–5 times faster than any of them alone. 

Default fastp parameters were used i.e. adapter trimming is enabled by default, quality filtering (phred quality score of 15 as default), length filtering, global trimming which means trim all reads in the front or tail.

By carefully analyzing the FastQC report and applying appropriate trimming techniques, you can significantly enhance the quality and usability of your sequencing data.

```bash
cat samples.txt
```
The output should look like this, without **(_fastq)** extension
```bash
tumor
normal
```

```bash
mkdir -p trims

for id in $(cat samples.txt)
do
    fastp -i reads/${id}_1.fastq -I reads/${id}_2.fastq.gz -o trims/${id}_1_trim.fastq -O trims/${id}_2_trim.fastq
done
```
## :feet: **STEP2:** REFERENCE BASED ALIGNMENT

### a). Alignment of Reads To The Reference

Reference-based alignment methods utilize the sequence for each read (and it’s mate in paired-end data) to find potential mapping locations by exact match or scoring sequence similarity. Mapping locations indicate transcripts of origin, and the number of reads originating from a given transcript inform of how much of the transcript was present in the sample.

Popular tools for **DNA** alignment include;
- BWA (Burrows-Wheeler Aligner)
- Bowtie

Popular tools for **RNA** alignment include;
- HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts)
- STAR (Spliced Transcripts Alignment to a Reference)
- Kalisito


We used **BWA** for this pipeline since were dealing with DNA sequences and also the added advantages of BWA over other tools such as its speed, accuracy and efficiency in handling short reads.

**BWA** is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
It consists of three algorithms: **BWA-backtrack**, **BWA-SW** and **BWA-MEM**. 
The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for highquality queries as it is faster and more accurate.
BWA-MEM also has better performance than BWAbacktrack for 70-100bp Illumina reads.

We shall therefore use **BWA MEM** to perform alignment.

Before performing the alignment, the reference must first be **Indexed**. This enables fast, efficient, and accurate alignment of sequencing reads by creating a searchable data structure that optimizes memory usage and computational load. This is done using **bwa index** command below.
```bash
# Downloading the reference genome
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

# Indexing 
bwa index Homo_sapiens_assembly38.fasta
```
Follow this [**LINK**](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) if you would prefer to use a different reference sequence

```bash
# Creating bams directory to store alignment results
mkdir -p ./bams

# Align paired end reads to the reference
for id in $(cat samples.txt) 
do
    bwa mem \
        -t 4 -aM ./ref/Homo_sapiens_assembly38.fasta \
        -R "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA" \
        ./trim/${id}_1_trim.fastq \
        ./trim/${id}_2_trim.fastq |\
        samtools view -bS - > ./bams/${id}.bam
done
# -t for number of threads
# -R for specifying read groups
# -aM for outputting all alignments and marking the short ones
# -bS for specifying the input sam file and converting it to a bam file
```
The output are unsorted bam files as below
```bash
samtools view -h ./bams/normal.bam | head -n4
```
![image](https://hackmd.io/_uploads/SkentfBKC.png)

The next step is to index the output bam after sorting it. 
Sorting a bam file is to organize the aligned reads either by their genomic coordinates or by read names, which is essential for downstream analyses such as duplicate marking, variant calling, and efficient data retrieval.
```bash
for id in $(cat samples.txt)
do
    samtools sort ./bams/${id}.bam -o ./bams/${id}_sorted.bam 
done
```
After sorting the bam files we then index them.
Indexing sorted bam files enables quick retrieval of specific regions of the genome, significantly speeding up data access and enabling efficient querying of aligned sequencing data. We also generated the alignment summary statistics for the bam files using **Samtools FlagStat** command
```bash
for id in $(cat samples.txt)
do
    samtools index ./bams/${id}_sorted.bam
    samtools flagstat ./bams/${id}_sorted.bam > ./bams/${id}.stats.txt || exit 1
done
```

The next step is to mark sequence duplicates. This is done using GATK (Genome Analysis Toolkit).
GATK is a software suite developed by the Broad Institute for analyzing high-throughput sequencing data. It offers a wide range of utilities that perform a number of tasks such as quality control, data processing, duplicate marking, base recalibration, variant calling among others

### b). Mark Sequence Duplicates 
The **GATK MarkDuplicates** tool identifies and tags duplicate reads (reads originating from a single DNA fragment) in a BAM or SAM file. These duplicates arise during sample preparation or from optical sensor errors. This tool works by comparing sequences at the 5' positions, uses molecular barcodes if available, and marks duplicates in the SAM flags field with a hexadecimal value of 0x0400 (decimal 1024), producing a new SAM or BAM file with these annotations
```bash
# Making directory to store the output of gatk MarkDuplicates
mkdir -p ./gatk

for id in $(cat samples.txt)
do
    gatk MarkDuplicates \
        --INPUT ./bams/${id}_sorted.bam \
        --OUTPUT ./gatk/${id}_mkdp.bam \
        --METRICS_FILE ./gatk/${id}_metrics.txt \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true
done
```
### c). Base Quality Score Recalibration (BQSR)
The **gatk BaseRecalibrator** tool aims at correcting the systematic bias that affect the assignment of base quality scores by the sequencer. BQSR is done in two steps:

1. First step is the generation of recalibration table which consist of calculating error empirically and finding patterns in how error varies with basecall features over all bases. This step requires knownsites which can be downloaded from [LINK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false))
```bash
# Making directory to store the recalibration table
mkdir -p ./recal

# Downloading known sites
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

# Unzipping them
gunzip 1000G_phase1.snps.high_confidence.hg38.vcf.gz
gunzip Homo_sapiens_assembly38.known_indels.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Variable
known1="1000G_phase1.snps.high_confidence.hg38.vcf"
known2="Homo_sapiens_assembly38.known_indels.vcf"
known3="Mills_and_1000G_gold_standard.indels.hg38.vcf"

# Indexing known sites
gatk IndexFeatureFile --input $known1
gatk IndexFeatureFile --input $known2
gatk IndexFeatureFile --input $known3
```

Generating recalibration tables
```bash
for id in $(cat samples.txt); do
    gatk BaseRecalibrator \
        -I ./gatk/${id}_mkdp.bam \
        -R ${ref} \
        --known-sites ./ref/$known1 \
        --known-sites ./ref/$known2 \
        --known-sites ./ref/$known3 \
        -O ./recal/${id}_recal_data.table
done
```


2. Then applying generated recalibration data to each individual basecall based on the patterns identified in the recalibration table by writing out a new BAM or CRAM file.

```bash
for id in $(cat samples.txt); do
    gatk ApplyBQSR \
        -R ${ref} \
        -I ./gatk/${id}_mkdp.bam \
        --bqsr-recal-file ./recal/${id}_recal_data.table \
        -O ./recal/${id}_recal.bam
done
```


## :feet: **STEP3:** VARIANT CALLING AND ANNOTATION
### a). Variant calling
This is a computational process used to identify genetic variations between individuals or populations. These variations can include single nucleotide polymorphisms (SNPs), insertions, deletions, and larger structural variants.

*Applications of variant calling:*
* Disease genetics: Identifying genetic variants associated with specific diseases.
* Population genetics: Studying genetic variation across populations and understanding evolutionary history.
* Personalized medicine: Tailoring medical treatments based on an individual's genetic makeup.
* Agricultural genomics: Improving crop yields and disease resistance through genetic modification.

*Tools used for variant calling:*

* [GATK](https://gatk.broadinstitute.org/hc/en-us) (Genome Analysis Toolkit): A popular tool for variant calling and analysis.
* [FreeBayes](https://github.com/freebayes/freebayes): A probabilistic variant caller.
* [VarScan](https://varscan.sourceforge.net/): A tool for detecting somatic variants in cancer genomes.
* [Samtools](https://github.com/samtools/samtools): A toolkit for manipulating and analyzing sequencing data, including variant calling.

We used GATK and VarScan on the same dataset for this analysis. By accurately identifying genetic variants, variant calling plays a crucial role in advancing our understanding of human genetics, disease biology, and population evolution.

We need to download a population variant file  and its index. This file is used to filter out potential germline mutations from the final list of somatic mutations. These are available at [LINK](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

NB: We create tsample.txt and nsamples.txt for tumor and non-tumor samples respectively
```bash
wget -P ref https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget -P ref https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

# Unzipping the vcf file downloaded above.
gunzip ./ref/af-only-gnomad.hg38.vcf.gz

pop="af-only-gnomad.hg38.vcf.gz"

# Creating variables with tumor and normal sample names respectively
tsample=$(cat tsamples)
nsample=$(cat nsamples)
mkdir -p VCFs

for id in $(cat sample)
do
  # Assuming names are separated by space in the respective files
  tumor=$(echo "$id" | cut -d ' ' -f 1)
  normal=$(echo "$nsample" | cut -d ' ' -f 1)
  #echo "$tumor,$normal"
  gatk  Mutect2 -R ref/$ref -I recal/${tumor}_tumor_recal.bam --tumor-sample ${tumor} -I $recal/${normal}_normal_recal.bam --normal-sample ${normal} -O VCFs/${tumor}_tumor_${normal}_normal.Mutect2.raw.vcf --germline-resource ref/$pop --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
done
```

### b). Variant filtering
**Variant filtering** is a crucial step in genomic analysis that involves removing low-quality or artifact variants from a set of identified variants. This process is essential for ensuring the accuracy and reliability of downstream analyses, such as disease association studies or population genetics.

*Common filtering criteria:*
* Quality score: Variants with low quality scores, often based on the base quality scores of the underlying reads, are typically filtered out.
* Depth of coverage: Variants with insufficient depth of coverage (i.e., the number of reads supporting the variant) may be filtered to reduce the likelihood of false positives.
* Allele frequency: Variants with very low or very high allele frequencies may be filtered, as they could be artifacts or common polymorphisms.
* Genotype quality: Variants with low genotype quality, indicating uncertainty in the assigned genotype, can be filtered.
* Hardy-Weinberg equilibrium (HWE): Variants that deviate significantly from HWE expectations may be filtered, as this could indicate genotyping errors or population stratification.
* Strand bias: Variants that show a strong bias towards one strand of DNA may be filtered, as this could be indicative of sequencing artifacts.
* Codon bias: Variants that disrupt the codon usage bias of a gene may be filtered, as they could be non-coding variants or sequencing errors.
* Known variants: Variants that have been previously reported in databases (e.g., dbSNP) can be filtered if they are not of interest.

Best filtering practices such as using multiple filtering criteria, evaluating the impact of filtering were considered. By carefully selecting and applying appropriate filtering criteria, we significantly improved the quality and reliability of our vgenomic analysis results.

```bash
# Making directory for filtered results
mkdir -p filtered

for id in $(cat sample)
do
    # Generating pileup files for Tumor and Matched Normals 
  	gatk GetPileupSummaries -I recal/${id}_Tumor_recal.bam -O filtered/${id}_Tumor.pileups.table --variant ./$k5 --intervals Illumina_Exome_TargetedRegions_v1.2.hg38.bed --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
	# Estimation of contamination
	gatk CalculateContamination -I filtered/${id}_Tumor.pileups.table -O filtered/${id}_Tumor_Normal.contamination.table --matched-normal filtered/${id}_Normal.pileups.table --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
	# Application of a first filter to variant calls"
	gatk FilterMutectCalls --variant VCFs/${tumor}_tumor_${normal}_normal.Mutect2.raw.vcf -O VCFs/${tumor}_Tumor_${normal}_Normal.Mutect2.oncefiltered.vcf --contamination-table filtered/${id}_${tumor}_Tumor_${normal}_Normal.contamination.table --reference $ref --tmp-dir temp_dir --java-options '-Xmx16g -Djava.io.tmpdir=tmp'
	# Removing variants that failed to pass the filter check
	awk '/^#/ {print $0; next} $7=="PASS" {print $0}' VCFs/${id}_${tumor}_Tumor_${normal}_Normal.Mutect2.oncefiltered.vcf > VCFs/${id}_Filtered.vcf
done
```

### c). Variant Annotation

**Variant annotation** is the process of adding biological context to identified genetic variants. This involves associating variants with genes, regulatory regions, and other relevant genomic features, as well as providing information about their potential functional consequences.

*Key steps in variant annotation:*

1. Gene mapping: Identify the genes that are located near or within the regions containing the variants.
2. Exon/intron prediction: Determine whether the variants are located within coding exons or non-coding introns.
3. Consequence prediction: Predict the potential functional consequences of the variants, such as missense mutations, nonsense mutations, frameshifts, in-frame deletions/insertions, and regulatory variants.
4. Population frequency: Determine the frequency of the variants in different populations.
5. Disease association: Identify any known associations between the variants and specific diseases or traits.
6. Functional annotation: Provide additional information about the functional context of the variants, such as their role in protein domains or regulatory motifs.

*Tools used for annotation:*

* [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/): A popular tool for annotating variants with information from multiple databases.
* [SnpEff](http://pcingola.github.io/SnpEff/): Another widely used tool for variant annotation.

*Importance of variant annotation:*

* Prioritization of variants: Annotation can help prioritize variants for further investigation based on their potential functional consequences.
* Interpretation of variants: Annotation provides context for understanding the significance of variants and their potential impact on health or disease.
* Functional validation: Annotation can guide the design of experiments to validate the functional consequences of variants.

By adding biological context to identified genetic variants, variant annotation plays a crucial role in advancing our understanding of these variants.


```bash

# Downloading the hg38 database
snpEff download -v hg38

for id in $(cat sample)
do
  # Annotation
  snpEff hg38 VCFs/${id}_Filtered.vcf > VCFs/${id}_Filtered_annotated_snpEff.vcf 
done
```
## Acknowledgements
We appreciate ACE-Uganda for letting us use their server, and we thank the Cancer Working Group and the GATK team for their support and helpful feedback. We’re especially grateful to Mr. Fredrick Kakembo for his leadership and guidance throughout this project. 


## Authors
N. Gloria,
O. Walter, 
S. Syrus

## GATK Team
![GATK_TEAM_2024](https://hackmd.io/_uploads/SJT8nuY6C.jpg)





## References
1. Ulintz, P. J., Wu, W., & Gates, C. M. (2019) Bioinformatics Analysis of Whole Exome Sequencing Data. Methods in molecular biology (Clifton, N.J.), 1881, 277–318. https://doi.org/10.1007/978-1-4939-8876-1_21
2. Poplin R, Ruano-Rubio V, DePristo MA, Fennell TJ, Carneiro MO, Van der Auwera GA, Kling DE, Gauthier LD, Levy-Moonshine A, Roazen D, Shakir K, Thibault J, Chandran S, Whelan C, Lek M, Gabriel S, Daly MJ, Neale B, MacArthur DG, Banks E. (2017) Scaling accurate genetic variant discovery to tens of thousands of samples bioRxiv, 201178. DOI: 10.1101/201178
3. Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. (2021) Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2, giab008, https://doi.org/10.1093/gigascience/giab008
4. Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). (2012) A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3; 6(2):80-92. PMID: 22728672
5. Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].
6. Shifu Chen. (2023) Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107
