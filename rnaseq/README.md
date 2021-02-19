# Make RNAseq map to scaffold contigs

## Map RNAseq reads to the PacBio assembly and create a map

```bash
ln -s $(readlink -f input/asm/purged.1.fasta) tmp/contigs.fa
# Map the RNAseq reads using bwa-mem
./map_rnaseq.sh

# Eliminate multi-mapping reads.
samtools view -H tmp/merged.bam > tmp/filtered.sam
samtools view -F 4 -q 3 tmp/merged.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' \
  >> tmp/filtered.sam

# Now we will run AUGUSTUS. Annotations are required by agouti.
cp -r /share/apps/sbcs/augustus/3.2.2/config tmp/augustus_3.2.2_config

# Run AUGUSTUS using fly as the model species.
augustus --AUGUSTUS_CONFIG_PATH=tmp/augustus_3.2.2_config --gff3=on \
  --species=fly tmp/contigs.fa > tmp/augustus.gff

# Run agouti
python ext/AGOUTI/agouti.py scaffold -assembly tmp/contigs.fa \
  -bam tmp/filtered.sam -gff tmp/augustus.gff \
  -outdir tmp/rnaseq_merged_filtered_agouti

# Convert to BED
samtools faidx tmp/contigs.fa
./make_rnaseq_map.rb tmp/contigs.fa.fai \
  tmp/rnaseq_merged_filtered_agouti/*.scaffolding_paths.txt > tmp/agouti.bed
```
