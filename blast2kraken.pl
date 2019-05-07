#!/usr/bin/env perl

# Give a Kraken-style report from a Blast output
#
# Xiangzhi Liang <xiangzhi.liang@simceredx.com>
# version 0.1  2019-05

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


my ($blast_result, $fasta, $taxdb, $min_ident, $min_length);
my $show_zeros = 0;
my $PROG = "blast2kraken.pl";

GetOptions(
  "help" => \&display_help,
  "x=s" => \$blast_result,
  "q=s" => \$fasta,
  "t=s" => \$taxdb,
  "min-ident=i" => \$min_ident,
  "min-length=i"=> \$min_length
) or usage();

usage() unless defined $blast_result;

sub usage {
  my $exit_code = @_ ? shift : 64;
  print STDERR "
Usage: blast2kraken.pl -x <blast result>  -q <fasta>  -t <taxdb> OPTIONS  > <kraken-style.out>

blast2kraken.pl creates Kraken-style reports from blast out files.

Options:
    -x  Blast            (REQUIRED) Blast result 
    -q  Fasta            (REQUIRED) Fasta input
    -t  TaxDB            (REQUIRED) Taxdb from taxonomy
    -min-ident  Score           Require a minimum identity score for reads alignment
    -min-length Score           Require a minimum lentgh for reads query

";
  exit $exit_code;
}

sub display_help {
  usage(0);
}

# taxonomy
my (%child_lists, %name_map, %rank_map, %parent_map);
print STDERR "Loading taxonomy ...\n";
load_taxonomy($taxdb,\%name_map, \%child_lists, \%rank_map, \%parent_map);

# fasta reads count
my $seq_count = readpipe("cat $fasta|grep '>'|wc -l");

# blast result
my (%hash,%taxo_counts,$prevReadID,$prevTaxID);
my $classified_count = 0;
open BLAST, "$blast_result"||die;
while (<BLAST>) {
    chomp;
    my @arr = split /\t/ ;
    my ($readID,$qlen,$seqID,$pident,$taxID)= @arr[0,1,2,5,15] ;
    next if defined $min_ident && $pident < $min_ident;
    next if defined $min_length && $qlen < $min_length;

    if( ! exists $hash{$readID}){
        $hash{$readID} = 1;
        $classified_count += 1;
        $taxID = 1 if ( !isTaxIDInTree( $taxID )) ;
        if ( ( defined $prevReadID ) && ( $readID eq $prevReadID ) ) {
            --$taxo_counts{$prevTaxID};
            $prevTaxID = lca($prevTaxID, $taxID);
            ++$taxo_counts{$prevTaxID};
        } else {
            ++$taxo_counts{$taxID};
            $prevTaxID = $taxID;
        }
        $prevReadID = $readID ;
    }
}

my %clade_counts = %taxo_counts;
dfs_summation(1);

for (keys %name_map) {
  $clade_counts{$_} ||= 0;
}
die "No sequence matches with given settings" unless $seq_count > 0;

# output
my $unclassified = $seq_count - $classified_count;
printf "%6.2f\t%d\t%d\t%s\t%d\t%s%s\n",
  $unclassified * 100 / $seq_count,
  $unclassified, $unclassified, "U",
  0, "", "unclassified";
dfs_report(1, 0);

######################
sub isTaxIDInTree {
  my $a = $_[0] ;

  while ( $a > 1 )
  {
    if ( !defined $parent_map{ $a } )
    {
      #print STDERR "Couldn't find parent of taxID $a - directly assigned to root.\n";
      return 0 ;
    }
    last if ( $a eq $parent_map{$a} ) ;
    $a = $parent_map{ $a } ;
  }
  return 1 ;
}

sub dfs_report {
  my $node = shift;
  my $depth = shift;
  if (! $clade_counts{$node} && ! $show_zeros) {
    return;
  }
  printf "%6.2f\t%d\t%d\t%s\t%d\t%s%s\n",
    ($clade_counts{$node} || 0) * 100 / $seq_count,
    ($clade_counts{$node} || 0),
    ($taxo_counts{$node} || 0),
    rank_code($rank_map{$node}),
    $node,
    "  " x $depth,
    $name_map{$node};
  my $children = $child_lists{$node};
  if ($children) {
    my @sorted_children = sort { $clade_counts{$b} <=> $clade_counts{$a} } @$children;
    for my $child (@sorted_children) {
      dfs_report($child, $depth + 1);
    }
  }
}

sub lca {
  my ($a, $b) = @_;
  return $b if $a eq 0;
  return $a if $b eq 0;
  return $a if $a eq $b;
  my %a_path;
  while ($a ge 1) {
    $a_path{$a} = 1;
    if (!defined $parent_map{$a}) {
      print STDERR "Couldn't find parent of taxID $a - directly assigned to root.\n";
      last;
    }
    last if $a eq $parent_map{$a};
    $a = $parent_map{$a};
  }
  while ($b > 1) {
    return $b if (defined $a_path{$b});
    if (!defined $parent_map{$b}) {
      print STDERR "Couldn't find parent of taxID $b - directly assigned to root.\n";
      last;
    }
    last if $b eq $parent_map{$b};
    $b = $parent_map{$b};
  }
  return 1;
}

sub rank_code {
  my $rank = shift;
  for ($rank) {
    $_ eq "species" and return "S";
    $_ eq "genus" and return "G";
    $_ eq "family" and return "F";
    $_ eq "order" and return "O";
    $_ eq "class" and return "C";
    $_ eq "phylum" and return "P";
    $_ eq "kingdom" and return "K";
    $_ eq "superkingdom" and return "D";
  }
  return "-";
}

sub dfs_summation {
  my $node = shift;
  my $children = $child_lists{$node};
  if ($children) {
    for my $child (@$children) {
      dfs_summation($child);
      $clade_counts{$node} += ($clade_counts{$child} || 0);
    }
  }
}


sub load_taxonomy {

  my($taxdb, $name_map, $child_lists, $rank_map, $parent_map) = @_;

  open TAXDB, "$taxdb"||die;
  while (<TAXDB>) {
    chomp;
    my @fields = split /\t/;
    my ($node_id, $parent_id, $name, $rank) = @fields[0,1,2,3];
    if ($node_id == 1) {
      $parent_id = 0;
    }
    $$name_map{$node_id} = $name;
    $$child_lists{$parent_id} ||= [];
    push @{ $$child_lists{$parent_id} }, $node_id;
    $$rank_map{$node_id} = $rank;
    $$parent_map{$node_id} = $parent_id;
  }
  close TAXDB;
}