#!/usr/bin/env ruby
# Copyright 2020 Anurag Priyam, Queen Mary U London

# Input is output of MSTmap
input_file = ARGV.shift

# Check the file exists and is not empty.
if !File.exist?(input_file) || File.zero?(input_file)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Parse input.
groups = []
IO.foreach(input_file) do |line|
    next if line[0] == ';' || line[0] == "\n"

    line = line.split
    if line[0] == 'group'
        groups << []
        next
    end

    *ref, pos = line[0].split('_')
    ref = ref.join('|')
    pos = pos.to_i
    cM = line[1].to_f
    groups.last << [ref, pos, cM]
end

# Remove spurious lgs.
groups.reject! { |group| group.length < 10 }

# lgs are sorted by number of markers. We want to sort by cM.
groups.sort_by! { |group| -group.last.last }

# We have double the number of expected linkage groups.
# Select the first of each pair.
groups.select!.with_index { |_, i| i.even? }

# Print
groups.each_with_index do |group, i|
    group.each do |marker|
        ref, pos, cM = marker
        puts [ref, pos, pos + 1, "lg#{i}:#{cM}"].join("\t")
    end
end
