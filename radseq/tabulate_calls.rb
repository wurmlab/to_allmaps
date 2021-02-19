#!/usr/bin/env ruby
# Copyright 2020 Anurag Priyam, Queen Mary U London
#
# Reads a VCF file and outputs for each sample the number of called sites,
# and the number homozygous and heterzygous calls.

require 'matrix'

# The input VCF file.
input_vcf = ARGV.shift

# Check the file exists and is not empty.
if !File.exist?(input_vcf) || File.zero?(input_vcf)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Read each line of genotype matrix from the VCF file and store it in an
# array. The array (of arrays) is later converted into a matrix.
genotypes = []
IO.foreach(input_vcf) do |line|
    # Skip all but the last line of the header.
    next if line[0..1] == '##'

    # Subsequent lines will either contain sample names or genotypes.
    # Information relevant to us is present from 10th column onwards.
    cols = line.split
    genotypes << cols[9..-1]
end

# Turn genotypes array (of arrays) into a matrix.
genotypes = Matrix[*genotypes]

# For each sample (column) output the number of called sites, the number
# of homzygous sites, and the number of heterozygous sites.
genotypes.column_vectors.each do |column|
    num_called = 0
    num_homozygous = 0
    column[1..-1].each do |genotype|
        gt = genotype.split(':').first
        next if gt.include? '.'
        num_called += 1
        num_homozygous += 1 if gt.split('/').uniq.length == 1
    end
    total_sites = column.size - 1
    perc_called = (num_called / total_sites.to_f).round(2)
    perc_homozygous = (num_homozygous / total_sites.to_f).round(2)
    
    puts [column[0], perc_called, perc_homozygous].join("\t")
end