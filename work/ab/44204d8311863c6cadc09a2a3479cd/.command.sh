#!/bin/bash -ue
bowtie2 -x Agy99 -1 P7741_R1_fastpTrimmed.fq.gz -2 P7741_R2_fastpTrimmed.fq.gz -S /dev/stdout |         samtools view -h -b - | samtools sort -o P7741_sorted.bam -
samtools index P7741_sorted.bam
