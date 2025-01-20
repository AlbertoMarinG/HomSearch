#!/bin/bash
# Example bash file to align ChIP-Seq fastq files to reference genome, then filter, sort and index
genfile="PATH-TO-GENOME-FILE"
genname="NAME-OF-GENOME-ASSEMBLY (e.g., hg38 or mm10)"
inR1="PATH-TO-READ1-FASTQ-FILE (e.g. RAD51_NoEP_R1.fastq.gz"
inR2="PATH-TO-READ2-FASTQ-FILE (e.g. RAD51_NoEP_R2.fastq.gz"
out="PATH-TO-OUTPUT-FILE (no extension) (e.g. RAD51_NoEP)"
bowtie2 -p 8 -q --local -X 1000 -x $genfile -1 $inR1 -2 $inR2 -S "$out"_"$genname".sam
samtools view -h -S -b -F 0x08 -q 25 "$out"_"$genname".sam > "$out"_"$genname"_unsorted.bam
samtools fixmate -m "$out"_"$genname"_unsorted.bam "$out"_"$genname"_fixmate.bam
samtools sort "$out"_"$genname"_fixmate.bam > "$out"_"$genname"_sorted.bam
samtools markdup -r "$out"_"$genname"_sorted.bam "$out"_"$genname"_final.bam
samtools index "$out"_"$genname"_final.bam
samtools flagstat "$out"_"$genname"_final.bam > "$out"_"$genname"_flagstats.txt
