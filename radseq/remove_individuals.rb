#!/usr/bin/env ruby
# Copyright 2020 Anurag Priyam, Queen Mary U London
#
# Removes the given individuals from a VCF files.

require 'matrix'

# The input VCF file and the file containing list of individuals to remove.
input_vcf = ARGV.shift
indiv_lst = ARGV.shift

# Check the files exist and are not empty.
if [input_vcf, indiv_lst].any? { |f| f.nil? || !File.exist?(f) || File.zero?(f) }
    puts 'Usage: remove_individuals.rb vcf_file ids_file'
    exit!
end

# Read in the file containing id of individuals to remove. The ids must
# match the column names in VCF
individuals_to_remove = File.readlines(indiv_lst).map!(&:chomp!)

# Go over each line of the VCF file, print the header lines and store
# the genotype lines in an array. The array is later converted into a
# matrix, so that we can easily retrieve the columns.
lines = []
IO.foreach(input_vcf) do |line|
    # Print all but the last line of the header.
    if line[0..1] == '##'
        puts line
        next
    end
    # Subsequent lines will either contain sample names or genotypes.
    # Information relevant to us is present from 10th column onwards.
    lines << line.split
end
matrix = Matrix.rows(lines)

# Remove given columns and recreate the matrix.
retained_columns = []
matrix.column_vectors.each_with_index do |column, i|
    next if i >=9 && individuals_to_remove.include?(column[0])
    retained_columns << column
end
matrix = Matrix.columns(retained_columns)

# Print
matrix.to_a.each do |row|
    puts row.join("\t")
end