#! /usr/bin/env perl
# Auth: Jennifer Chang
# Date: 2018/05/14

use strict;
use warnings;

# ===== Check ARGS
my $USAGE = "USAGE: $0 <input.gb> > <tabular|fasta>\n";
$USAGE = $USAGE . "    input.gb - genbank file\n";
$USAGE =
  $USAGE . "    Edit the printEntry function to select tab or fasta output\n";

if ( @ARGV != 1 ) {
    die $USAGE;
    exit;
}

# ===== Variables
my $fn = $ARGV[0];

my $unknown = "-";

my $col_country = "Collection_Country";
my $col_date    = "Collection_Date";
my $consrtm     = "Consrtm";
my $gb          = "GenBank";
my $host        = "Host";
my $iso_source  = "Isolation_Source";
my $sample      = "Sample_Identifier";
my $state       = "State";
my $strain      = "Strain";
my $subtype     = "Subtype";
my $segment     = "Segment";
my $finprot     = "";

# = Fasta output
my $fasta = "";
my $prot  = "";
my $seq   = -1;
my $yr    = "ColYear";
my $mo    = "ColMon";
my $dd    = "ColDay";

my %months = (
    "Jan", "01", "Feb", "02", "Mar", "03", "Apr", "04",
    "May", "05", "Jun", "06", "Jul", "07", "Aug", "08",
    "Sep", "09", "Oct", "10", "Nov", "11", "Dec", "12"
);

# ===== Print line of GB data
sub printEntry() {

    # = Recover A0 number from strain name if not in isoalte
    if ( $sample eq $unknown ) {
        $sample = ( split( /\//, $strain ) )[3];
        $sample =~ s/ //g;
    }

    $state       = ( split( /\//, $strain ) )[2]      || $state;
    $col_country = ( split( /:/,  $col_country ) )[0] || $col_country;

    # = Numerical dates
    my @temp = reverse( split( /-/, $col_date ) );
    if ( $temp[0] eq $col_date ) {
        $yr = $temp[0];
    }
    $yr       = $temp[0]     || $yr;
    $mo       = $temp[1]     || $mo;
    $dd       = $temp[2]     || $dd;
    $mo       = $months{$mo} || $mo;
    $col_date = join( "/", $yr, $mo, $dd );
    $col_date =~ s/\/00/\/01/g;

    $strain =~ s/ /_/g;
    $strain = ( split( /\(/, $strain ) )[0];

# = Tabular output
#    print join("\t",$sample,$strain,$host,$subtype,$col_date,$yr,$mo,$dd,$col_country,$iso_source,$gb),"\n";
#    print join("\t",$sample,$col_date,$yr,$mo,$dd,$state,$gb,$subtype,$strain,$iso_source),"\n";

    # = Fasta output
    $fasta =~ s/[0-9]//g;
    $fasta =~ s/ //g;
    $segment = "_$segment";

    print ">", join( "|", $gb, $strain, $subtype, $col_date ), "\n";
    print $finprot, "\n";

    # = Reset variables
    $col_country = $unknown;
    $col_date    = $unknown;
    $yr          = "0000";
    $mo          = "00";
    $dd          = "00";
    $consrtm     = $unknown;
    $gb          = $unknown;
    $host        = $unknown;
    $iso_source  = $unknown;
    $sample      = $unknown;
    $state       = $unknown;
    $strain      = $unknown;
    $subtype     = $unknown;
    $segment     = $unknown;
    $fasta       = "";
    $seq         = -1;
}

# ===== Main
my $fh;
open( $fh, "<:encoding(UTF-8)", $fn )
  or die "Could not open file '$fn'";

# = Reset Variable names for fasta output

while (<$fh>) {
    if (/^\/\//) {
        printEntry;
        $finprot = "";
    }
    if ( $seq == 1 ) {
        $fasta = $fasta . $_;
    }
    elsif ( $seq == 2 ) {
        $prot = $prot . $_;
        $prot =~ s/ //g;
        $prot =~ s/"//g;

        if (/(.+)"/) {
            $seq = -1;
            if ( length($prot) > length($finprot) ) {
                $finprot = $prot;
            }
        }
    }
    elsif (/LOCUS\s+(\S+)/) {
        $gb = $1;
    }
    elsif (/isolate="(.+)"/) {
        $sample = $1;
    }
    elsif (/host="(.+)"/) {
        $host = $1;
    }
    elsif (/country="(.+)"/) {
        $col_country = $1;
    }
    elsif (/isolation_source="(.+)"/) {
        $iso_source = $1;
    }
    elsif (/strain="(.+)"/) {
        $strain = $1;
    }
    elsif (/serotype="(.+)"/) {
        $subtype = $1;
    }
    elsif (/segment="(.+)"/) {
        $segment = $1;
    }
    elsif (/collection_date="(.+)"/) {
        $col_date = $1;
    }
    elsif (/^ORIGIN/) {
        $seq = 1;
    }
    elsif (/CONSRTM\s+(\S.+)/) {
        $consrtm = $1;
    }
    elsif (/translation="(.+)"/) {
        $prot = $1;
        if ( length($prot) > length($finprot) ) {
            $finprot = $prot;
        }
    }
    elsif (/translation="(.+)/) {
        $prot = "$1\n";
        $seq  = 2;
    }
}
