#!/usr/bin/env bash
set -euo pipefail

# fq1_file路径中不要出现空格！都是*.cleaned_1.fastq.gz格式，不再进行trim

fq1_file="EMS_4.28_fq1.txt"
fq2_file="EMS_4.28_fq2.txt"
threads=16
genome_fasta="/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bwa_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

echo -e "Start time: " $(date "+%H:%M:%S %Y-%m-%d")

echo "--------------------------------------------"
echo "1. Perform fastqc"
set +eu && source activate seq && set -eu
[ -d 1_fastqc ] || mkdir 1_fastqc
cat $fq1_file $fq2_file | xargs fastqc -t $threads -o 1_fastqc
multiqc -f --fullnames -o 1_fastqc 1_fastqc
set +eu && conda deactivate && set -eu
echo -e "Complete time: " $(date "+%H:%M:%S %Y-%m-%d")

echo -e "\n\n--------------------------------------------"
echo "2. Mapping to genome with bwa"
set +eu && source activate gatk3 && set -eu
[ -d 2_mapping ] || mkdir 2_mapping
for i in $(cat $fq1_file); do
	prefix=${i##*/}
	prefix=${prefix%.cleaned_1.fastq.gz}
	fq2=${i/%_1.fastq.gz/_2.fastq.gz} ##################
	bwa mem -t $threads -R '@RG\tID:'$prefix'\tSM:'$prefix'\tPL:illumina' $genome_fasta \
	"$i" "$fq2" | samtools view -@ $threads -bhS -o "$prefix".bam -
	samtools sort -@ $threads -o 2_mapping/"$prefix"_sorted.bam "$prefix".bam
	samtools index -@ $threads 2_mapping/"$prefix"_sorted.bam
	rm "$prefix".bam
done
set +eu && conda deactivate && set -eu
echo -e "Complete time: " $(date "+%H:%M:%S %Y-%m-%d")

echo -e "\n\n--------------------------------------------"
echo "3. Deduplication and realignment"
set +eu && source activate gatk3 && set -eu
[ -d 3_dedup_realign ] || mkdir 3_dedup_realign
for i in $(cat $fq1_file); do
	prefix=${i##*/}
	prefix=${prefix%.cleaned_1.fastq.gz}
	sortbam="2_mapping/${prefix}_sorted.bam"
	picard MarkDuplicates -I $sortbam -O 3_dedup_realign/"$prefix"_dedup.bam \
	-M 3_dedup_realign/"$prefix"_metrics.txt --REMOVE_DUPLICATES true
	samtools index -@ $threads 3_dedup_realign/"$prefix"_dedup.bam
	java -Xmx10g -jar /home/lining/miniconda2/envs/gatk3/opt/gatk-3.8/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator -R $genome_fasta -I 3_dedup_realign/"$prefix"_dedup.bam -o 3_dedup_realign/"$prefix".intervals
	java -Xmx10g -jar /home/lining/miniconda2/envs/gatk3/opt/gatk-3.8/GenomeAnalysisTK.jar \
	-T IndelRealigner -R $genome_fasta -I 3_dedup_realign/"$prefix"_dedup.bam -targetIntervals 3_dedup_realign/"$prefix".intervals \
	-o 3_dedup_realign/"$prefix"_dedup_realign.bam
	samtools index -@ $threads 3_dedup_realign/"$prefix"_dedup_realign.bam
	rm 3_dedup_realign/"$prefix"_metrics.txt
	rm 3_dedup_realign/"$prefix"_dedup.bam
	rm 3_dedup_realign/"$prefix"_dedup.bam.bai
	rm 3_dedup_realign/"$prefix".intervals
done
set +eu && conda deactivate && set -eu
echo -e "Complete time: " $(date "+%H:%M:%S %Y-%m-%d")

echo -e "\n\n--------------------------------------------"
echo "4. Call variant and filtering"
[ -d 4_variant ] || mkdir 4_variant
for i in $(cat $fq1_file); do
	prefix=${i##*/}
	prefix=${prefix%.cleaned_1.fastq.gz}
	bamfile="3_dedup_realign/${prefix}_dedup_realign.bam"
	samtools mpileup --adjust-MQ 50 --fasta-ref $genome_fasta \
	--BCF --uncompressed --output-tags DP,AD,INFO/AD --ignore-RG \
	/bios-store1/LN/Project/20240329_YHR/hen1_8_dedup_realign.bam $bamfile | \
	bcftools call --output-type v --format-fields GQ --variants-only --multiallelic-caller \
	--output 4_variant/"$prefix".vcf
	bcftools filter --threads $threads -g 5 -G 10 \
	-e '%QUAL < 20 || INFO/MQ < 20 || %MIN(FORMAT/DP) < 4' -s LowQual -O z \
	-o 4_variant/"$prefix"_filtered.vcf.gz 4_variant/"$prefix".vcf
	bcftools index 4_variant/"$prefix"_filtered.vcf.gz
done
echo -e "Complete time: " $(date "+%H:%M:%S %Y-%m-%d")
