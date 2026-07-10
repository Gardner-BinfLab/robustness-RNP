#!/usr/bin/perl

#script to extract probable variants from a SAM file derived from mutagenesis experiments

use warnings;
use strict;
use Getopt::Long;

my ($refFile, $barcodeFile, $samFile, $from, $to, $conservative, $phred, $help, $verbose);
my $minFreq = 2;
$phred = 12; #low default phred score threshold 
&GetOptions(
    "b|barcodes=s"       => \$barcodeFile,
    "r|ref=s"            => \$refFile, 
    "s|sam=s"            => \$samFile,
    "m|minFreq=s"        => \$minFreq,
    "f|from=s"           => \$from,
    "t|to=s"             => \$to,
    "c|conservative"     => \$conservative,
    "p|phred=s"          => \$phred,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}

######################################################################
#Read barcodes:
my (%barcodes, %barcodesIds, @barcodes, @barcodeIds); #
my $barcodeLen = 0;
open(BIN, "< $barcodeFile") or die "FATAL: failed to open barcode file [$barcodeFile]\n$!";
while(my $bin = <BIN>){
    chomp($bin);
    next if $bin =~ /^#/; #skip comments
    if($bin =~ /^(\S+)\t(\S+)$/){
	push(@barcodes,   $2);
	push(@barcodeIds, $1) if not defined($barcodesIds{$1});
	$barcodes{$2}   =$1;
	$barcodesIds{$1}=$2   if not defined($barcodesIds{$1}); #only record the first barcode sequence for each ID (others are low freq variants)
	$barcodeLen = max(length($2), $barcodeLen); 
    }
    else {
	print "WARNING: barcode file line [$bin], does not conform!\n";
    }
}
close(BIN);
######################################################################
#Read reference file:
my $id;
my $refSeq = '';
open(RIN, "< $refFile") or die "FATAL: failed to open reference sequence file [$refFile]\n$!";
while(my $rin = <RIN>){
    chomp($rin);
    #print "rin:[$rin]\n";
    if($rin =~ /^>\s*(\S+)/){
	$id = $1;
    }
    else {
	$refSeq .= $rin;
    }
}
close(RIN);

if(not defined $from and not defined $to){
    $from = 1;
    $to = length($refSeq);
}

#print "#Printing variation between [$from,$to].\n";

######################################################################
#Open SAM file, decode, initialise...:
#assign sequences to barcodes in first round
my %seqId2barcode; 
open(SIN, "< $samFile") or die "FATAL: failed to open reference sequence file [$samFile]\n$!";
while(my $sin = <SIN>){
    
    chomp($sin); 
    my @sin = split(/\t/, $sin);
    if(scalar(@sin) >= 11){

	next if defined($seqId2barcode{$sin[0]});

	$seqId2barcode{$sin[0]} = 'unknown';
	foreach my $bar (@barcodes){
	    if ( $sin[9] =~ /^$bar/ ){
		$seqId2barcode{$sin[0]} = $barcodes{$bar};
		last;
	    }
	}
    }    
}
close(SIN);


# Col Field   Type Regexp/Range Brief description
# 1   QNAME   String [!-?A-~]{1,254}             Query template NAME
# 2   FLAG    Int [0, 2^{16} − 1]                bitwise FLAG
# 3   RNAME   String \*|[:rname:∧*=][:rname:]*   Reference sequence NAME11
# 4   POS     Int [0, 231 − 1]                   1-based leftmost mapping POSition
# 5   MAPQ    Int [0, 28 − 1]                    MAPping Quality
# 6   CIGAR   String \*|([0-9]+[MIDNSHPX=])+     CIGAR string
# 7   RNEXT   String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
# 8   PNEXT   Int [0, 231 − 1]                   Position of the mate/next read
# 9   TLEN    Int [−231 + 1, 231 − 1]            observed Template LENgth
# 10  SEQ     String \*|[A-Za-z=.]+              segment SEQuence
# 11  QUAL    String [!-~]+                      ASCII of Phred-scaled base QUALity+33

#https://stackoverflow.com/questions/20251428/perl-converting-special-characters-to-corresponding-numerical-values-in-qualit
my %decodeQuals;
my $phredString = q{!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~};
foreach my $asc ( split( //, $phredString) ) {
    $decodeQuals{$asc} = ord($asc) - 33;
}
#Q = −10 log10 Pr{base is wrong}
#Pr{base is wrong} = 10 ^ (Q / (-10))
#
#e.g. Pr{base is wrong} = 1/1000 --> Q = 30

#initialise:
my @cooccurrence;
for (my $i = $from; $i<=$to; $i++){
    for (my $j = $from; $j<=$to; $j++){
	$cooccurrence[$i][$j] = 0; 
    }
}
######################################################################
#Parse CIGARS & count variants
my %variants; 
my $cnt=0;
my %seenVariants; #hash of hashes: restrict to variants observed on both pairs!
my %variantCounts;
open(SIN, "< $samFile") or die "FATAL: failed to open reference sequence file [$samFile]\n$!";
my $mutFile = $samFile . ".variation";
open my $mutFH, "> $mutFile" or die "FATAL: failed to open mutation file for writing [$mutFile]\n$!";
while(my $sin = <SIN>){
    chomp($sin); 
    my @sin = split(/\t/, $sin);
    if(scalar(@sin) >= 11){
	#        MAPQ
	next if ($sin[4] < 20 or $sin[4] == 255);       #skip poor mapping quality (usually due to a high number of errors)
	next if ($sin[8] == 0);                         #skip unmapped seqs
	next if ($seqId2barcode{$sin[0]} eq 'unknown'); #skip reads with unknown barcode ADD A COUNTER OF THESE!!!!!
	#                               POS        TLEN
	my $subRefSeq = substr($refSeq, $sin[3]-1, abs($sin[8])+1);
	#print "substr($refSeq, $sin[3]-1, abs($sin[8])+1);\n";
	#                                                       SEQ      CIGAR
	my ($gappedRef, $gappedSeq) = cigarSeqs2aln($subRefSeq, $sin[9], $sin[5]); 
	
	#print $mutFH "# $sin[0]\n";
	
	my $varCount = aln2vars($sin[0], $gappedRef, $gappedSeq, $sin[10], \%variants, \%seenVariants, $sin[3]-1, $barcodeLen, $seqId2barcode{$sin[0]}, \%decodeQuals, $from, $to, \@cooccurrence, $mutFH, $conservative, $phred);	
	#printf "$sin[5]\n$id/%d-%d\t$gappedRef\n$sin[0]\t$gappedSeq\nvarCount=[$varCount]\n", $sin[3]-1, abs($sin[8]);
	$variantCounts{$varCount} = 0 if not defined($variantCounts{$varCount});
	$variantCounts{$varCount}++;
    }
    $cnt++;
    
    #last if ($cnt>200);
}
close(SIN);
close($mutFH);
######################################################################
#PRINT MAIN RESULTS TO STDOUT:
print "pos\trefnuc\tquerynuc\t"; #headers 
foreach my $barcodeId (@barcodeIds){
    print "$barcodeId\t";
}
print "total\tvpos\n";

#VARIANTS	
foreach my $var (sort keys %variants) {
    my $sum=0;
    my ($vpos, $rnuc, $qnuc) = split(/:/, $var);
    next if (($vpos < $from) or ($to < $vpos) );
    my $printPos = $vpos - $from + 1;
    print "$printPos\t$rnuc\t$qnuc\t";    
    foreach my $barcodeId (@barcodeIds){
	$variants{$var}{$barcodeId}=0 if not defined($variants{$var}{$barcodeId});
        printf "%d\t", $variants{$var}{$barcodeId};
	$sum+=$variants{$var}{$barcodeId};
    }
    print "$sum\t$vpos\n";
}

###################################
#Print summary stats:
my $totalVariantCounts=0;
foreach my $cnt (keys %variantCounts) {
    $totalVariantCounts += $variantCounts{$cnt};
}
my $statFile = $samFile . ".stats";
open(UT, "> $statFile");
print UT "numVarsPerSeq\tcount\tpercentage\tcumulativePercentage\n";
my $cumPercentVariantCounts=0;
foreach my $cnt (sort { $a <=> $b } keys %variantCounts) {
	my $percentVariantCounts     = 100*$variantCounts{$cnt}/$totalVariantCounts;
	$cumPercentVariantCounts += $percentVariantCounts;
	printf UT "$cnt\t$variantCounts{$cnt}\t%0.2f\t%0.2f\n", $percentVariantCounts, $cumPercentVariantCounts;
}
close(UT);
###################################

my $cutFile = $samFile . ".co-occurance";
open(UT, "> $cutFile");
#Print co-occurance matrix

print UT "row\t"; #colnames
for (my $i = $from; $i<$to; $i++){
    my $printPos = $i - $from + 1;    
    print UT "$printPos\t";
}
printf UT "%d\n", $to - $from + 1;  #"$to\n";

for (my $i = $from; $i<=$to; $i++){
    my $printPos = $i - $from + 1;    
    print UT "$printPos"; #rowname
    for (my $j = $from; $j<=$to; $j++){
	print UT "\t" . $cooccurrence[$i][$j]; 
    }
    print UT "\n";
}
close(UT);

exit(0); 

######################################################################

sub aln2vars {
    
    my ($seqId, $reference, $query, $qual, $ptrVarHash, $ptrSeenVarHash, $pos, $barcodeLen, $barcodeGrp, $decodeQuals, $from, $to, $cooccurrencePtr, $mutFilePtr, $conservative, $phred) = @_;

    #printf "#aln2vars: Printing variation between [$from,$to]. Ref length=[%d]\n", length($reference);

    my @reference = split(//, $reference); 
    my @query     = split(//, $query);
    my @qual      = split(//, $qual);
    #walk along alignment, find the mismatches,
    #--check that the sequence quality is still high****
    #--count the number of observations of each inclusion variant
    #--keys are on the position+type of mismatch and  barcode
    #--e.g. '122A:C' & 'presort'
    my $varCount=0;
    my $numRefGaps=0;
    my @sites;
    my @cooccurrence = @{$cooccurrencePtr};
    for (my $i = 0; $i < scalar(@reference); $i++){
	
	$numRefGaps++ if (defined($reference[$i]) && not is_nucleotide($reference[$i]));
	if( defined($reference[$i]) && defined($query[$i])    && defined($qual[$i]) #do the data points all exist 
	    && $reference[$i] ne $query[$i]                                         #is this a bona fide change?
	    && is_nucleotide($reference[$i])   &&   is_nucleotide($query[$i])       #are the characters legal nucleotides?
	    && ($decodeQuals->{$qual[$i]}) >= $phred                                #is the sequence quality high?
	    ){                                       

	    
	    my $nucPos = $pos + $i - $numRefGaps;
	    my $printPos = $nucPos - $from + 1;
	    
	    printf $mutFilePtr "$seqId\t$printPos\t$reference[$i]\t$query[$i]\t%d\t$nucPos", $decodeQuals->{$qual[$i]} if ($printPos >= 0);
	    if (defined($conservative) && not defined($ptrSeenVarHash->{$seqId}{$nucPos}) ){
		$ptrSeenVarHash->{$seqId}{$nucPos} = 1;
		printf $mutFilePtr "\tF\n" if ($printPos >= 0);
		next;
	    }
	    printf $mutFilePtr "\tT\n" if ($printPos >= 0);
	    my $varString = sprintf "%d:$reference[$i]:$query[$i]", $nucPos;
	    if(defined( $ptrVarHash->{$varString}{$barcodeGrp} )){
		($ptrVarHash->{$varString}{$barcodeGrp})++;
	    }
	    else{
		$ptrVarHash->{$varString}{$barcodeGrp} = 1; 
	    }

	    if ( $from <= $nucPos && $nucPos <= $to ){#is the variant in the right location?
		push(@sites, $nucPos);
		$varCount++;
	    }
	}
    }
    #print $mutFilePtr "#$varCount\n";
    for (my $i = 0; $i<scalar(@sites); $i++){
	for (my $j = 0; $j<scalar(@sites); $j++){
	    $cooccurrence[$sites[$i]][$sites[$j]]++;
	}
    }
    return $varCount; 
}


######################################################################
#https://davetang.org/muse/2011/01/28/perl-and-sam/
sub cigarSeqs2aln {

    my $reference = shift;
    my $query     = shift;
    my $cigar     = shift;
    
    #from the CIGAR string, fill in the insertions and deletions
    #print "cigar[$cigar] ref[$reference] q[$query]\n";
    my $position = 0;
    while ($cigar !~ /^$/){
	if ($cigar =~ /^([0-9]+)([MIDS])/){
            my $cigar_part = $1 . $2;
	    #print "cigElementLen[$1]\tcigElementType[$2]\n";
	    if(length($reference) > int($1)){
		if ($cigar_part =~ /(\d+)M/){
		    $position += $1;
		} elsif ($cigar_part =~ /(\d+)I/){
		    my $insertion = '-' x $1;
		    #print "substr($reference,$position,0,$insertion); [$cigar_part]\n";
		    substr($reference,$position,0,$insertion);
		    $position += $1;
		} elsif ($cigar_part =~ /(\d+)D/){
		    my $insertion = '-' x $1;
		    substr($query,$position,0,$insertion);
		    $position += $1;
		} elsif ($cigar_part =~ /(\d+)S/){
		    print "Not ready for this!\n";
		    #my $insertion = 'x' x $1;
		    #substr($new_ref,$position,0,$insertion);
		    #$position += $1;
		}
	    }
	    else {
		printf "Warning: reference is too short for the cigar[$cigar] int($1) refLen[%d] ref[$reference] q[$query]\n", length($reference);
	    }
            $cigar =~ s/$cigar_part//;
	} else {
            print "Warning: unexpected cigar string: [$cigar]\n";
	}
    }
    
    return ($reference, $query); 
}

######################################################################
#returns true if input character is a nucleotide (IUPAC codes):
sub is_nucleotide {
    my $a = shift;
    if (defined($a) ){
	$a =~ tr/a-z/A-Z/;
	if (length($a) && ($a =~ /[ACGUTRYWSMKBDHVN]/) ){
	    return 1;
	}
	else {
	    return 0;
	}
    }
    else {
	return 0;
    }
}

######################################################################
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}

######################################################################
sub help {
    print STDERR <<EOF;

$0: extract the probable variants from a SAM file derived from mutagenesis experiments.

Usage:   $0 -i <seqfile> 
  Options:
    -b|--barcodes [str]       Tab delimited file of IDs & barcode seqs. 
    -r|--ref      [str]       Reference sequence file that was used to generate the SAM file.
    -s|--sam      [str]       SAM file. 
    -f|--from     [int]       Only report variation after  this coordinate in the reference sequence.
    -t|--to       [int]       Only report variation before this coordinate in the reference sequence.
    -c|--conservative         Only report variation found on both the forward and reverse paired-end reads. 
    -h|--help                 Show this help.
    -v|--verbose              Print lots of stuff.

  Dependencies:
  Perl 

  TODO:

EOF
}



