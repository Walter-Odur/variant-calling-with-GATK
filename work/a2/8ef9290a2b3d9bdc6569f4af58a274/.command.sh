#!/bin/bash -ue
fastp         -i P7741_R1.fastq.gz         -I P7741_R2.fastq.gz         -o P7741_R1_fastpTrimmed.fq.gz         -O P7741_R2_fastpTrimmed.fq.gz         --detect_adapter_for_pe
