#!/usr/bin/env ruby

fasta_index = ARGV.shift
agouti_paths = ARGV.shift

intervals = {}
IO.foreach(fasta_index) do |line|
  line = line.split
  contig_id = line[0]
  last_base = line[1].to_i - 1
  intervals[contig_id] = [0, last_base]
end

paths = []
IO.foreach(agouti_paths) do |line|
  next if line =~ /^>/
  paths << line.strip.split(',')
end

paths.each_with_index do |contigs, i|
  contigs.each_with_index do |contig, j|
    contig_id = contig.gsub /\-/, ''
    contig_interval = intervals[contig_id]
    contig_interval.reverse! if contig[0] == '-'

    contig_interval.each_with_index do |pos,k|
      puts [contig_id,  pos, pos + 1, "ag#{i+1}:#{j}.#{k}"].join("\t")
    end
  end
end
