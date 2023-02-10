#! /usr/bin/env perl

use strict;
use warnings;

# ============================ Check ARGS and set Variables
my $USAGE = "USAGE: $0 <positions> <alnment.fa> > subseq.txt\n";
$USAGE=$USAGE."    positions - comma separated positions \n";


if(@ARGV<2){
    die $USAGE;
    exit;
}

# ============================ Initial Variables
#my ($fn2,$col,$fn1)=@ARGV;
my ($positions,$fn)=@ARGV;

my @pos=split(/,/,$positions);
# ======================== Variables

# ======================== Functions

my $header="";
my $seq="";
#my $start=30;
#my $end=60;

sub printSubSeq{
    foreach my $i (@pos){
	if($i =~ /\-/){
	    my ($start,$end)=split(/\-/,$i);
	    print substr($seq,$start-1,$end-$start+1);
	}else{
	    print substr($seq,$i-1,1);
	}
    }
#    print substr($seq,$start,$end-$start+1);
    print "\t${positions}\t$header\n";
}

# ======================== Process
my $fh;
open($fh, '<:encoding(UTF-8)', $fn)
    or die "Could not open file '$fn' $1";

while(<$fh>){
    chomp;
    if(/^>(.+)/){
	if(length($header)>0){
	    printSubSeq;
	}
	$header=$1;
	$seq="";
    }else{
	$seq=$seq.$_;
	$seq=~s/ //g;  # remove any spaces
    }
}

if(length($header)>0){
    printSubSeq;
}
