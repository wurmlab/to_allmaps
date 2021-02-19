#!/usr/bin/env perl
# Version 1.0
# - https://github.com/tangerzhang/my_script/blob/master/bionano2Allmaps.pl,
#   Accessed: 2019-04-11 12:24 pm BST
#
# Version 1.1 (Anurag Priyam)
# - Convert coordinates to int before printing
# - Take XMAP and key file as arguments
# - Output BED instead of CSV
# - Prefix 'bn' to marker ids

use FindBin;
use Getopt::Std;
use File::Basename;

getopts "x:k:o:";

if ((!defined $opt_x) || (!defined $opt_k) || (!defined $opt_o)) {
  die "************************************************************************
  Usage: perl $0 -x xmap -i ref.fasta -e BspQI
  -h : help and usage
  -x : xmap from hybrid assembly
  -k : _key file from in silico digestion of the reference
  -o : output file name
  ************************************************************************\n";
}else{
  print "*******************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger, tanger.zhang\@gmail.com\nPriyam, anurag.priyam\@qmul.ac.uk\n";
  print "RUNNING...\n";
  print "*******************************************************************************\n";

}
$xmap = $opt_x;
$key  = $opt_k;
$bed = $opt_o;

open(IN, $key) or die"";
while(<IN>){
  chomp;
  next if(/#/);
  next if(/CompntId/);
  @data = split(/\s+/,$_);
  $id   = $data[0];
  $scaf = $data[1];
  $iddb{$id} = $scaf;
}
close IN;

open(IN, $xmap) or die"";
while(<IN>){
  chomp;
  next if(/#/);
  @data      = split(/\s+/,$_);
  $chr_id    = $data[1];
  $chr_posia = $data[3];
  $chr_posib = $data[4];
  $lg        = $data[2];
  $lg_posia  = $data[5];
  $lg_posib  = $data[6];
  $infordb{$lg}->{$lg_posia} = $chr_id."_".$chr_posia;
  $infordb{$lg}->{$lg_posib} = $chr_id."_".$chr_posib;
}
close IN;
open(OUT, "> ${bed}") or die"";
foreach $lg(sort {$a<=>$b} keys %infordb){
  foreach $lg_posi (sort {$a<=>$b} keys %{$infordb{$lg}}){
    ($chr_id,$chr_posi) = split(/_/,$infordb{$lg}->{$lg_posi});
    $chrn = $iddb{$chr_id};

    # Turn to integer
    $chr_pos1 = int($chr_posi);
    # pos+1 for 3rd column
    $chr_pos2 = $chr_pos1 + 1;
    print OUT "$chrn\t$chr_pos1\t$chr_pos2\tbn$lg:$lg_posi\n";
  }
}
close OUT;
