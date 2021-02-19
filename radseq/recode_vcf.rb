#!/usr/bin/env ruby
# Copyright 2020 Anurag Priyam, Queen Mary U London

# The input VCF file.
vcf_file = ARGV.shift
max_dp = Integer(ARGV.shift)
min_gq = Integer(ARGV.shift)
min_num_observations = Integer(ARGV.shift)
min_maf = Float(ARGV.shift)

# Check the files exist and are not empty.
if !File.exist?(vcf_file) || File.zero?(vcf_file)
  puts 'File does not exist, or it is empty.'
  exit!
end

# Parse VCF.
IO.foreach(vcf_file) do |line|
    # Skip all but the last line of the header.
    next if line[0..1] == '##'

    # Subsequent lines will either contain sample names or genotypes.
    # Information relevant to us is present from 10th column onwards.
    line = line.split
    if line[0][0] == '#'
        puts line[9..-1].unshift('#').join("\t")
        next
    end

    # Extract ref, pos and calls.
    ref = line[0]
    pos = line[1]
    calls = line[9..-1]

    # How many individuals do we have?
    num_individuals = calls.length

    # What are we calling this site?
    site_id = "#{ref.gsub('|', '_')}_#{pos}"

    # Let's get mean depth, mean gq, and genotype counts for this
    # site. In addition to '.', also ignore heterygous calls in
    # the calculation. While at it, also transform the calls to
    # 0, 1 or nil (where nil is for ./. and heterozygous call).
    mean_dp = 0
    mean_gq = 0
    geno_count = Hash.new(0)
    calls.map! do |call|
        gt, dp, gq = call.split(':').values_at(0, 1, 3)
        gt = gt.split('/').uniq
        dp = dp.to_i
        gq = gq.to_i
        next if gt.include? '.'
        next if gt.length > 1
        gt = gt.first
        mean_dp += dp
        mean_gq += gq
        geno_count[gt] += 1
        gt
    end
    # Remove sites that have no valid genotype call or are not polymorphic.
    # geno_count.length will respectively be 0 and 1 in these cases.
    next unless geno_count.length == 2

    # How many observations do we have for this site?
    num_observations = geno_count.map(&:last).reduce(:+)

    # What is the mean DP and GQ of this site?
    mean_dp = mean_dp / num_observations
    mean_gq = mean_gq / num_observations

    # Minor allele frequency:
    geno_count = geno_count.sort_by(&:last)
    minor_allele_frequency = (geno_count[0][1] / num_observations.to_f).round(2)    

    $stderr.puts [num_observations, mean_dp, mean_gq, minor_allele_frequency].join("\t")
    next if num_observations < min_num_observations ||
            mean_dp > max_dp || mean_gq < min_gq ||
            minor_allele_frequency < min_maf
     
    # Re-code the site.
    # a1 = '0'
    # a2 = '1'
    a1 = rand(2).to_s
    a2 = ['0', '1'].find { |f| f != a1 }
    # a1 = geno_count[1][0] # major allele
    # a2 = geno_count[0][0] # minor allele
    calls.map! do |call|
        case call
        when a1 then 'A'
        when a2 then 'B'
        else 'U'
        end
    end

    puts calls.unshift(site_id).join(' ')
end
