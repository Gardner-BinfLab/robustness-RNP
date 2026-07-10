#!/usr/bin/perl

#script to analyse sam2var.pl output -- train phred thresholds for variant calls outside of overlapping regions 

use warnings;
use strict;
use Getopt::Long;
 
my $fdr = 0.05;
my ($help, $verbose, $input, $from, $to, $phred); 

$from  = 1;
$to    = 200;
$phred = 12; 

&GetOptions(
    "f|fdr=s"            => \$fdr,
    "i|input=s"          => \$input,
    "f|from=s"           => \$from,
    "t|to=s"             => \$to,
    "p|phred=s"          => \$phred,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}


if (not defined $input or not -s $input) {
    print "-i <filename> is a compulsory input requirement!\n";
    &help();
    exit(1);
}


my $utFile = $input . '.fixed';
open(IN, "tac $input | " ) or die "FATAL: failed to open pipe [tac $input]\n[$!]";
open(UT, "> $utFile") or die "FATAL: failed to open file for writing [$utFile]\n[$!]";
my %seen;
#Fix false 'F' labels 
while(my $in = <IN>){
    chomp($in);
    my @in = split(/\t/, $in);
    my $key = $in[0] . '/' . $in[1];
    if($in[6] eq 'T'){
	$seen{$key} = 1;
    }
    elsif (defined($seen{$key})){
	$in[6] = 'T';
    }

    if( $from <= $in[1] && $in[1] <= $to){
	my $printLine = join("\t", @in);
	print UT $printLine . "\n";
    }
    
}
close(IN);
close(UT);



system("tac $utFile > blah && mv blah $utFile"  );

#Read counts in:
open(IN, "cat $utFile | sort -k5,5nr -k7,7d | cut -f 5,7 | uniq -c  | " ) or die "FATAL: failed to open pipe [sort $utFile]\n[$!]";
my %scoreCounts;
my ($trueTot, $falseTot) = (0, 0);
while(my $in = <IN>){
    chomp($in);
    my @in = split(/\s+/, $in);
    $scoreCounts{$in[2]}{$in[3]} = $in[1];

#    print "0[$in[0]] 1[$in[1]] 2[$in[2]] 3[$in[3]]\n";
#    exit(0);
    
    if($in[3] eq 'T'){
	$trueTot+= $in[1];
    }
    elsif($in[3] eq 'F'){
	$falseTot+= $in[1];
    }
}
close(IN);

my ($tp, $fp, $fn, $tn) = (0, 0, 0, 0);
print "phred\tsens\tspec\tfdr\tfpr\tf1\tmcc\ttp\tfp\tfn\ttn\n";
foreach my $thresh (sort {$b <=> $a} keys %scoreCounts ){
    $scoreCounts{$thresh}{'T'} = 0 if (not defined $scoreCounts{$thresh}{'T'});
    $scoreCounts{$thresh}{'F'} = 0 if (not defined $scoreCounts{$thresh}{'F'});

    $tp += $scoreCounts{$thresh}{'T'};
    $fp += $scoreCounts{$thresh}{'F'};
    
    $fn = $trueTot  - $tp;
    $tn = $falseTot - $fp;
    
    my ( $sens, $spec, $fdr, $fpr, $f1, $mcc ) = accuracyStats($tp, $fp, $fn, $tn);
    
    printf "$thresh\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%d\t%d\t%d\t%d\n",
	$sens, $spec, $fdr, $fpr, $f1, $mcc,
	$tp, $fp, $fn, $tn;
}



exit(0); 

sub accuracyStats {
    
    my ($tp, $fp, $fn, $tn) = @_; 

    my ( $sens, $spec, $fdr, $fpr, $f1, $mcc ) = (0, 0, 0, 0, 0, 0);
    if(isPositiveInteger($tp) && isPositiveInteger($fp) && isPositiveInteger($fn) && isPositiveInteger($tn)){
    
	$sens = $tp / ($tp + $fn) if ($tp + $fn);
	$spec = $tn / ($tn + $fp) if ($tn + $fp);
	$fdr  = $fp / ($fp + $tp) if ($fp + $tp);
	$fpr  = $fp / ($fp + $tn) if ($fp + $tn);
	$f1   = 2 * $tp / (2 * $tp + $fp + $fn) if ($tp + $fp + $fn);
	my $numerator = sqrt( ($tp + $fp)*($tp + $fn)*($tn + $fp)*($tn + $fn) );
	$mcc  = ($tp * $tn - $fp * $fn) / $numerator if ($numerator);
    }
    return ( $sens, $spec, $fdr, $fpr, $f1, $mcc );
}

######################################################################
sub isPositiveInteger {
    my $num = shift;
    if ($num=~/^\d+$/) { 
        return 1; 
    }
    else {
        return 0;
    }
}



######################################################################
sub help {
    print STDERR <<EOF;

    $0: Analyse sam2var.pl output -- train phred thresholds for variant calls outside of overlapping regions.
	Return a matrix of thresholds for each 

Usage:   $0 -i <seqfile> 
Options:
    -f|--fdr    <float>       We want thresholds that limit the FDR to a specified value (Default: 0.05).
    -i|--input  <filename>    A filename for *.variation output from sam2var.pl. 
    -j|--from   <int>         Restrict analysis to co-ordinates between "from" and "to".
    -t|--to     <int>         Restrict analysis to co-ordinates between "from" and "to".
    -h|--help                 Show this help.
    -v|--verbose              Print lots of stuff.

  Dependencies:
  Perl 

  TODO:

EOF
}



