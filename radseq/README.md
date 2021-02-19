## Usage

To use the scripts you need Ruby 2.0 or higher.

The first step is to obtain number of genotype calls, mean depth, mean genotype quality, and minor allele frequency for each site.

```
./recode_vcf.rb vcf_file 2> site_stats > /dev/null
```

Next, visualise the site data using `analyse_sites.Rmd` R markdown notebook and decide the following cutoffs:

1. Minimum number of observations per site
2. Maximum mean site depth
3. Minimum mean genotype quality per site
4. Minimum minor allele frequency

Next, filter and recode the VCF:

```
./recode_vcf.rb vcf_file depth_cutoff genotype_quality_cutoff min_observations_cutoff maf_cutoff > genotypes_matrix 2> /dev/null
```

Convert the genotypes matrix to MSTmap:

```
./to_mstmap.rb genotypes_matrix p_value_cutoff missing_threshold > mstmap_input
```

Run MSTmap

```
MSTmap mstmap_input > lgs 2> mstmap.log
```

Convert to ALLMAPS

```
./to_allmaps.rb lgs > map.bed
```
