set -exo pipefail

bwa index tmp/contigs.fa

ls input/rnaseq | grep -E '\.(R1|R2)' | cut -d'.' -f1 | sort | uniq > tmp/carlos_paired
for i in $(cat tmp/carlos_paired); do
  bwa mem -t 36 -M tmp/contigs.fa input/rnaseq/${i}.R1.fastq.gz input/rnaseq/${i}.R2.fastq.gz | samtools view -b > tmp/$i.bam
done

ls input/rnaseq | grep -E '_1|_2' | cut -d'_' -f1 | sort | uniq > tmp/sra_paired
for i in $(cat tmp/sra_paired); do
  bwa mem -t 36 -M tmp/contigs.fa input/rnaseq/${i}_1.fastq.gz input/rnaseq/${i}_2.fastq.gz | samtools view -b > tmp/$i.bam
done

ls input/rnaseq | grep -vE '_1|_2|\.R1|\.R2' | cut -d'.' -f1 > tmp/sra_interleaved
for i in $(cat tmp/sra_interleaved); do
  bwa mem -t 36 -M -p tmp/contigs.fa input/rnaseq/${i}.fastq.gz | samtools view -b > tmp/$i.bam
done

samtools merge -@36 -rn - tmp/*.bam > tmp/merging
mv tmp/merging tmp/merged.bam