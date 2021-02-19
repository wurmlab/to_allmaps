#!/usr/bin/env ruby
# Copyright 2020 Anurag Priyam, Queen Mary U London

# The input VCF file.
input_file = ARGV.shift
cut_off_p_value = ARGV.shift
missing_threshold = ARGV.shift

# Check the file exists and is not empty.
if !File.exist?(input_file) || File.zero?(input_file)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Helpers.
def hamming(line1, line2)
    line1 = line1[1..-1]
    line2 = line2[1..-1]
    h = 0
    line1.zip(line2).each do |g1, g2|
        next if g1 == 'U' || g2 == 'U'
        h = h + 1 if g1 != g2
    end
    h
end
def flip(line)
    line.map do |elem|
        case elem
        when 'A' then 'B'
        when 'B' then 'A'
        else elem
        end
    end
end

# Parse input.
input = Hash.new { |h,k| h[k] = [] }
num_sites = -1
IO.foreach(input_file) do |line|
    num_sites += 1
    line = line.split
    line[0] = 'locus' if line[0] == '#'
    input[line[0].split('_').first] << line
end
num_indviduals = input['locus'][0].length - 1

# Sort positions per contig.
input.each do |contig, sites|
    sites.sort_by! do |site|
        site.first.split('_').last.to_i
    end
end

# Phase.
input.each do |contig, lines|
    # Skip if it is the header row or the contig has only one marker (we can't
    # phase a single marker).
    next if contig == 'locus' || lines.length == 1

    # Flip.
    i = 0
    j = 1
    loop do
        l1 = lines[i]
        l2 = lines[j]

        if hamming(l1, flip(l2)) < hamming(l1, l2)
            lines[j] = flip(l2)
        end

        i = i + 1
        j = j + 1
        break if j == lines.length
    end
end

# Flip and double the sites.
input.each do |contig, lines|
    next if contig == 'locus'
    lines.concat(lines.map {|l| flip(l)})
end

# Print MSTmap header.
puts <<-HEADER
population_type DH
population_name LG
distance_function kosambi
cut_off_p_value #{cut_off_p_value}
no_map_dist 15.0
no_map_size 2
missing_threshold #{missing_threshold}
estimation_before_clustering no
detect_bad_data yes
objective_function COUNT
number_of_loci #{num_sites * 2}
number_of_individual #{num_indviduals}
HEADER

# Print marker id and genotypes table.
input.each do |contig, lines|
    lines.each do |line|
        puts line.join(' ')
    end
end
