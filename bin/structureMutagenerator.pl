#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use Statistics::RankCorrelation;
use File::Copy;
use Data::Dumper;
    
my ($infile, $isProtein, $help, $verbose);
my ($numMutants, $mutRate, $indelRate, $outdir) = (1000, 0.01, 0.01, '');
&GetOptions (
    'h'           => \$help,
    'help'        => \$help,
    'p'           => \$isProtein,
    'i=s'         => \$infile,
    'm=s'         => \$numMutants,
    'indelrate=s' => \$indelRate,
    'mrate=s'     => \$mutRate,
    'dir=s'       => \$outdir,
    'v|verbose'   => \$verbose
        );
if ($help)
{
print STDERR <<EOF;

structureMutagenerator.pl: 
A programme to generate mutants for either a ncRNA or mRNA sequence, 
    fold the sequences (with either RNAfold -p for RNA, I-TASSER for protein),
    compute a correlation coefficient with structure probabilities between the 
    wild-type and mutant. 

    
input: Fasta RNA sequence, either a ncRNA or mRNA 
output: 
     
Options:
    -h or --help    This help screen
    -p              Is this a mRNA or not? 
    -i              A fasta file
    -m              Number of mutations to generate (default: 1000)
    -indelrate      The INDEL rate (default: 0)
    -mrate          Mutation rate (default: 1/100)
    -v              Verbose, print lots of stuff (debugging)
  Usage:
    structureMutagenerator.pl -p  -i [fasta file]
    structureMutagenerator.pl     -i [fasta file]

   TODO:


EOF
exit;
}

my ($name,$seq);
die "input file [$infile] doesn't exist! Use '-i <filename>'!" if (! -s $infile);
open(IN, "< $infile");
while(my $in = <IN>){
    chomp($in);
    #next if ($in =~ /^\#/);
    if($in =~ /^>(\S+)/){
	$name=$1; 
	#print "[$2]\n";
    }
    elsif(defined($in) && length($in) && ($in =~ /[ACGUTRYWSMKBDHVN]/i)){
	$seq.=$in; 
    }
}
close(IN);
$seq =~ tr/a-z/A-Z/;
$seq =~ tr/T/U/;


#print ">$name\n$seq\n";
my $probabilitiesWT;
if(not defined($isProtein)){
    #ncRNA
    $probabilitiesWT = RNAfoldP($seq); 
    if(-s "dot.ps"){
	copy("dot.ps","data/$outdir/ss-predictions/dot.ps-mrate". 0.00 . "-indelrate" . 0.00 . "-num" . 0) or die "Copy failed: $!";
    }#

    print "ID\tstructCorrSpearman\tstructCorrPearson\tnumPointMutatns\tnumIns\tnumDels\tnumAlignedResidues\tgappedMutantSeq\tgappedWTSeq\n";
    printf "$name\t1.00\t1.00\t0\t0\t0\t%d\t$seq\t$seq\n", length($seq);
    for (my $i=0; $i<$numMutants; $i++ ){

	my ($gappedMutSeq,$gappedSeq,$numMuts,$numIns,$numDels) = mutateSequence($seq, $mutRate, $indelRate);	
	my $mutSeq = $gappedMutSeq; 
	$mutSeq =~ s/\-//g;
	my $probabilitiesMut = RNAfoldP($mutSeq);
	my ($alignedProbabilitiesWT, $alignedProbabilitiesMut, $numResidues) = alignArraysRNA($gappedMutSeq,$gappedSeq,$probabilitiesMut, $probabilitiesWT );
	copy("dot.ps","data/$outdir/ss-predictions/dot.ps-mrate". $mutRate . "-indelrate" . $indelRate . "-num" . $i) or die "Copy failed: $!";
	copy("alignedRNAprobs.tsv", "data/$outdir/ss-predictions/alignedRNAprobs.tsv-mrate". $mutRate . "-indelrate" . $indelRate . "-num" . $i) or die "Copy failed: $!";

	#print "[@{$alignedProbabilitiesWT}]\n[@{$alignedProbabilitiesMut}]\n";
	my $corr = Statistics::RankCorrelation->new( $alignedProbabilitiesWT, $alignedProbabilitiesMut, sorted => 1 );
	my $pcorr = pearsonCorr($alignedProbabilitiesWT, $alignedProbabilitiesMut); 
	
	printf "mutseq$i\t%0.2f\t%0.2f\t$numMuts\t$numIns\t$numDels\t$numResidues\t$gappedMutSeq\t$gappedSeq\n", $corr->spearman, $pcorr;

	unlink("dot.ps","rna.ps", "RNAfold-p.out","RNAfold-p.err");
	
    }
}
else {
    #mRNA/PROTEIN! 
    
    my $protSeq = translateMRNA($seq);
    
    print ">wt (nuc)\n$seq\n"      if defined($verbose);
    print ">wt (prot)\n$protSeq\n" if defined($verbose);
    my $ssArrayWT = [];
    ($probabilitiesWT,$ssArrayWT) = ProtfoldP($protSeq,$ssArrayWT); 
    printf "probsWT:@{$probabilitiesWT}\n" if defined($verbose);
    
    copy("seq.dat.ss","data/$outdir/ss-predictions/seq.dat.ss-mrate". 0.00 . "-indelrate" . 0.00 . "-num" . 0) or die "Copy failed: $!";
    unlink("blast.out","freqccw","freqccwG","in-seq.fasta","mtx","output1.ss","output2.ss","output3.ss","output4.ss","output5.ss","output6.ss","output7.ss","profw","psitmp.aux","psitmp.chk","psitmp.fasta","psitmp.mn","psitmp.mtx","psitmp.pn","psitmp.sn","pssm.txt");
    
    print "ID\tstructCorrSpearman\tstructCorrPearson\tnumPointMutatns\tnumIns\tnumDels\tnumAlignedResidues\tgappedMutantSeq\tgappedWTSeq\ttranslatedMutantSeq\n";
    printf "$name\t1.00\t1.00\t0\t0\t0\t%d\t$seq\t$seq\t$protSeq\n", length($seq);
    for (my $i=0; $i<$numMutants; $i++ ){
	
	my ($gappedMutSeq,$gappedSeq,$numMuts,$numIns,$numDels) = mutateSequence($seq, $mutRate, $indelRate);
	my $mutSeq = $gappedMutSeq; 
	$mutSeq =~ s/\-//g;
	print "#" x 70 . "\n"  if defined($verbose);
	print ">mut$i (nuc)\n$mutSeq\n" if defined($verbose);
	my $mutProtSeq = translateMRNA($mutSeq);
	print ">mut$i (prot)\n$mutProtSeq\n\nGapped:\n"        if defined($verbose);
	print "WT  $gappedMutSeq\nMT  $gappedSeq\n\n" if defined($verbose);
	if(($numMuts+$numIns+$numDels) == 0){
	    #If no mutations, then save lots of compute!
	    #could hash protein seqs and reuse mutns too...
	    printf "mutseq$i\t%0.2f\t%0.2f\t$numMuts\t$numIns\t$numDels\t%d\t$gappedMutSeq\t$gappedSeq\t$mutProtSeq\n", 1.000, 1.000, length($seq);
	    next;
	}
	my ($probabilitiesMut,$ssArrayMut) = ProtfoldP($mutProtSeq,$ssArrayWT);
	printf "ungapped protWT[ @{$probabilitiesWT } ]\t{%d}\n",   scalar(@{$probabilitiesWT})  if (defined($verbose)); 
	printf "ungapped protMT[ @{$probabilitiesMut} ]\t{%d}\n\n", scalar(@{$probabilitiesMut}) if (defined($verbose));
	
	my ($alignedProbabilitiesWT, $alignedProbabilitiesMut, $numResidues) = alignArraysProtein($gappedMutSeq,$gappedSeq,$mutProtSeq, $protSeq, $probabilitiesMut, $probabilitiesWT );
	printf "alignedWT:[@{$alignedProbabilitiesWT}] {%d}\nalignedMT:[@{$alignedProbabilitiesMut}] {%d}\n", scalar(@{$alignedProbabilitiesWT}), scalar(@{$alignedProbabilitiesMut}) if (defined($verbose));
	
	copy("seq.dat.ss","data/$outdir/ss-predictions/seq.dat.ss-mrate". $mutRate . "-indelrate" . $indelRate . "-num" . $i                     ) or die "Copy failed: $!";
	copy("alignedPROTprobs.tsv", "data/$outdir/ss-predictions/alignedPROTprobs.tsv-mrate". $mutRate . "-indelrate" . $indelRate . "-num" . $i) or die "Copy failed: $!";
	
	my $corr = Statistics::RankCorrelation->new( $alignedProbabilitiesWT, $alignedProbabilitiesMut, sorted => 1 );   
	#print "Statistics::RankCorrelation:" . Data::Dumper->Dump([$corr]) . "\n" if (defined($verbose));
	
	my $lcorr = ($corr->spearman);
	if ( scalar(@{$alignedProbabilitiesMut}) == 1 or length($mutProtSeq) == 1 ){ #weird edge case -- RankCorrelation says R = 1.0 when N=1.0.
	    $lcorr = 0.00;
	}
	my $pcorr = pearsonCorr($alignedProbabilitiesWT, $alignedProbabilitiesMut); 
	printf "mutseq$i\t%0.2f\t%0.2f\t$numMuts\t$numIns\t$numDels\t$numResidues\t$gappedMutSeq\t$gappedSeq\t$mutProtSeq\n", $lcorr, $pcorr;
	unlink("blast.out","freqccw","freqccwG","in-seq.fasta","mtx","output1.ss","output2.ss","output3.ss","output4.ss","output5.ss","output6.ss","output7.ss","profw","psitmp.aux","psitmp.chk","psitmp.fasta","psitmp.mn","psitmp.mtx","psitmp.pn","psitmp.sn","pssm.txt");
	
    }

    
}


exit();
######################################################################
#
sub alignArraysProtein {
    my ($gappedMutSeq,$gappedSeq,$mutProtSeq, $wtProtSeq, $probabilitiesMut, $probabilitiesWT ) = @_;
    #NB: $probabilitiesMut, $probabilitiesWT are concatenated coil, alpha and beta probabilities
    my (@alignedProbsMut, @alignedProbsWT);
    my @gappedMutSeq = split(//,$gappedMutSeq);
    my @gappedSeq    = split(//,$gappedSeq);
    my @mutProtSeq   = split(//,$mutProtSeq);
    my @wtProtSeq    = split(//,$wtProtSeq);
    my $cnt=0;
    
    my ($mutNucCnt, $mutProtCnt, $wtNucCnt, $wtProtCnt) = (0,0,0,0);
    my %seenIndices; #only report non-redundant values 
    printf "WARNING: WT sequence length (= %d) and number probability in array (%d) don\'t match!\n", length( $wtProtSeq), scalar(@{$probabilitiesWT})  if( length($wtProtSeq)  != scalar(@{$probabilitiesWT} )  );  
    printf "WARNING: MT sequence length (= %d) and number probability in array (%d) don\'t match!\n", length($mutProtSeq), scalar(@{$probabilitiesMut}) if( length($mutProtSeq) != scalar(@{$probabilitiesMut}) );  

    for (my $i = 0; $i < length($gappedMutSeq); $i++){
	
	if(isNucleotide($gappedMutSeq[$i]) && isNucleotide($gappedSeq[$i]) && 
	   isNumeric($probabilitiesMut->[$mutProtCnt]) && isNumeric($probabilitiesWT->[$wtProtCnt]) 
	    && not defined($seenIndices{"$mutProtCnt:$wtProtCnt"})){
	    #Aligned nucleotides
	    push(@alignedProbsMut, ($probabilitiesMut->[$mutProtCnt]) );
	    push(@alignedProbsWT,  ($probabilitiesWT->[  $wtProtCnt]) );
	    $seenIndices{"$mutProtCnt:$wtProtCnt"}=1;
	    $cnt++;
	}
	#This was a really bad idea -- created strong negative correlations which made R^2 > 0 
	#elsif(isNucleotide($gappedSeq[$i]) && not defined($mutProtSeq[$mutProtCnt]) ){ 
	#    #catch truncated mutant sequences:
	#    push(@alignedProbsMut, (0.000)                           );
	#    push(@alignedProbsWT,  ($probabilitiesWT->[$wtProtCnt] ) );
	#    $cnt++;
	#}
	
	$mutNucCnt++  if (isNucleotide($gappedMutSeq[$i]));
	$wtNucCnt++   if (isNucleotide($gappedSeq[$i]   ));
	$mutProtCnt++ if ((($mutNucCnt)>0) && ((($mutNucCnt) % 3) == 0) );
	$wtProtCnt++  if (( ($wtNucCnt)>0) && (( ($wtNucCnt) % 3) == 0) );
    }
    
    open(UT, "> alignedPROTprobs.tsv");
    for (my $i = 0; $i < $cnt; $i++){
	printf UT "%0.9f\t%0.9f\n", $alignedProbsWT[$i], $alignedProbsMut[$i];
    }
    close(UT);
    
    return(\@alignedProbsWT, \@alignedProbsMut, $cnt); 
}

######################################################################
#
sub alignArraysRNA {
    my ($gappedMutSeq,$gappedSeq,$probabilitiesMut, $probabilitiesWT )=@_;
    my (@alignedProbsMut, @alignedProbsWT);
    my @gappedMutSeq = split(//,$gappedMutSeq);
    my @gappedSeq    = split(//,$gappedSeq);
    my $cnt=0;
    for (my $i = 0; $i < length($gappedMutSeq); $i++){
	if(isNucleotide($gappedMutSeq[$i]) && isNucleotide($gappedSeq[$i]) && isNumeric($probabilitiesMut->[$i]) && isNumeric($probabilitiesWT->[$i]) ){
	    push(@alignedProbsMut, $probabilitiesMut->[$i]);
	    push(@alignedProbsWT,  $probabilitiesWT->[$i] );
	    $cnt++;
	}
    }
    
    open(UT, "> alignedRNAprobs.tsv");
    for (my $i = 0; $i < $cnt; $i++){
	printf UT "%0.9f\t%0.9f\n", $alignedProbsWT[$i], $alignedProbsMut[$i];
    }
    close(UT);
    return(\@alignedProbsWT, \@alignedProbsMut, $cnt); 
}

######################################################################
#
sub RNAfoldP{
    my $seq = shift; 

    system("echo $seq | RNAfold -p > RNAfold-p.out 2> RNAfold-p.err") and die "FATAL: failed to run system command: [echo $seq | RNAfold -p > RNAfold-p.out 2> RNAfold-p.err]\n[$!]";
    
    my @probsArray = (0) x length($seq);
    open(DOT, "< dot.ps");
    while(my $dp = <DOT>){
	if($dp =~ /^(\d+)\s+(\d+)\s+(\S+)\s+ubox/){
	    $probsArray[$1-1]+=($3*$3); #sqrt(p(i,j)) stored in PS file. 
	    $probsArray[$2-1]+=($3*$3); 
	}
	
    }
    close(DOT);
    
    return \@probsArray;
}

######################################################################
#
sub ProtfoldP{
    my ($seq, $structureArray) = @_; 
    my @structureArray = @{$structureArray};
    
    my $file = "in-seq.fasta";
    open(UT, "> in-seq.fasta");
    printf UT ">prot\n$seq\n";
    close(UT);
    
    system("~/inst/I-TASSER5.1/PSSpred/mPSSpred.pl $file ~/inst/I-TASSER5.1 ~/inst/I-TASSER5.1 > /dev/null") 
	and print "FATAL: failed to run: [~/inst/I-TASSER5.1/PSSpred/mPSSpred.pl $file  ~/inst/I-TASSER5.1 ~/inst/I-TASSER5.1]\n[$!]";

    my @probsArray = ( (0) x (length($seq)) );
    #hash to convert secondary structure assignment to column of file
    my %ss2col = (
	'C' => 0,
	'H' => 1,
	'E' => 2
	);
    open(SS, "< seq.dat.ss");
    while(my $ss = <SS>){
	if($ss =~ /^\s+(\d+)\s+(\w)\s+(\w)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    my @probs=($4,$5,$6); #coil helix beta
	    my $ssAss = $3;
	    if(defined($structureArray[$1-1])){
		$ssAss = $structureArray[$1-1];
	    }
	    $probsArray[$1-1] += $probs[ $ss2col{$ssAss} ];
	    push(@structureArray, $3);
	}
    }
    close(SS);
    
    return (\@probsArray, \@structureArray);
}


######################################################################

sub translateMRNA {

    my $nucSeq = shift; 
    $nucSeq =~ s/a-z/A-Z/g;
    $nucSeq =~ s/U/T/g;
    my $len = length($nucSeq);
    my $protSeq='';
    
    # if ( ($len % 3) ){
    # 	print "WARNING: length of [$nucSeq] is not a factor 3! Frameshifts may cause errors!\n";
    # }
    
    my @nucSeq = split(//, $nucSeq); 
    for (my $i=0; $i<$len; $i+=3){

	if(defined($nucSeq[$i+2])){
	    my $aa = codon2aa( $nucSeq[$i] . $nucSeq[$i+1] . $nucSeq[$i+2] );
	    last if ($aa  eq '*'); #Truncate sequence if a stop is hit! 
	    $protSeq .= $aa;
	}
	
    }
    
    #print "  protseq:[$protSeq]\n" if defined($verbose);
    return $protSeq;
    
}



######################################################################
#Given a sequence, generate a mutant with a given error rate (option for indels too)
sub mutateSequence{
    my ($seq, $mutRate, $indelRate) = @_; 
    my @seq = split(//, $seq);
    my @mutSeq = ();
    my ($numMuts,$insTot,$delTot)=(0,0,0);
    my $gappedSeq = $seq; 
    my $num;
    for (my $i=0; $i<scalar(@seq); $i++ ){
	#Point mutation
	if(isNucleotide($seq[$i]) ){
	    $num = rand(1);
	    if($num < $mutRate){
		$seq[$i] = mutateChar($seq[$i]);
		push(@mutSeq, $seq[$i]); 
		$numMuts++; 
	    }
	    else {
		push(@mutSeq, $seq[$i]); 
	    }
	}

	# #IF INDEL -- INSERT OR DELETE? -- return a "gapped" version of both original and mutant sequence
	$num = rand(1);
	if($num < $indelRate){
	    $num = rand(1);
	    my $indelLength=int(rand(5))+1;
	    if($num < 0.5){#INSERTION
		my $ins = '';
		
		for(my $j=0; $j<$indelLength; $j++){ 
		    $ins .= grabNuc();
		}
		push(@mutSeq, $ins);
		#Jesus. This is ugly:
		#Gapping sequence so that can return an alignment:
		$gappedSeq = substr($gappedSeq,0,($i+$insTot+1)) . ('-' x $indelLength) . substr($gappedSeq,($i+$insTot+1), length($gappedSeq)-1);
		
		$insTot+=$indelLength;
	    }
	    else {#DELETION
		next if (($i+$indelLength) > scalar(@seq)-1); #Can't delete sequence that isn't there! 
		for(my $j=0; $j<$indelLength; $j++){ 
		    push(@mutSeq, '-');
		}
		$i+=$indelLength;
		$delTot+=$indelLength;
	    }
	}
    }
    my $mutSeq=join('', @mutSeq);
#Forcing mutations, bad idea:
#    if(($numMuts+$delTot+$insTot)==0){#infinite loop if mutation & indel rate are 0...
#	($mutSeq,$gappedSeq, $numMuts, $insTot, $delTot)=mutateSequence($seq, $mutRate, $indelRate);
#    }
    #print "WARNING: mutant sequence length not equal to gapped seq! [length($mutSeq) != length($gappedSeq)]\nmutd:$mutSeq\ngapd:$gappedSeq\n" if (length($mutSeq) != length($gappedSeq));
    return ($mutSeq, $gappedSeq, $numMuts, $insTot, $delTot);
}

######################################################################
#
sub mutateChar {
    my $char = shift;
    my @nucs = qw(A C G U);
    my $newChar = $char;
    while($newChar eq $char){
	$newChar = $nucs[int(rand($#nucs))];
    }
    return $newChar; 
}

######################################################################
#
sub grabNuc {
    my @nucs = qw(A C G U);
    return $nucs[int(rand(4))]; 
}

######################################################################
#returns true if input character is a nucleotide (IUPAC codes):
sub isNucleotide {
    my $a = shift;
    
    if (defined($a)){
	$a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ACGUTRYWSMKBDHVN]/) ){
	return 1;
    }
    else {
	return 0;
    }
    
}

######################################################################
#returns true if input character is an amino-acid (IUPAC codes):
sub isAA {
    my $a = shift;
    
    if (defined($a)){
	$a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ACDEFGHIKLMNPQRSTVWXY]/) ){
	return 1;
    }
    else {
	return 0;
    }
    
}


######################################################################
sub isNumeric {
    my $num = shift;
    if (defined($num) && $num=~/^-?\d+\.?\d*$/) { 
	return 1; 
    }
    else {
	return 0;
    }
}

######################################################################
#Internet is out -- implementing my own Pearson correlation function!
sub pearsonCorr {
    my ($x, $y) = @_;
    my @x=@{$x};
    my @y=@{$y};

    if (scalar(@x) != scalar(@y) or scalar(@x) < 2){
	printf STDERR "WARNING: pearsonCorr arrays [@x] and [@y] are differnt length or have less than 2 elements!\n";
	return 0.0;
    }
    #my $xm = mean($x);
    #my $ym = mean($y);

    #my ($sumx2,$sumy2, $sumxy)  = (0.0,0.0,0.0);
    my ($sumxsqr)            = (0.0);
    my ($sumysqr)            = (0.0);
    my ($sumxy,$sumy, $sumx) = (0.0,0.0,0.0);

    my $n = scalar(@x);
    for (my $i = 0; $i<$n; $i++){

	if (isNumeric($x[$i]) && isNumeric($y[$i])){
	    $sumxsqr += $x[$i] * $x[$i];
	    $sumysqr += $y[$i] * $y[$i];
	    $sumxy   += $x[$i] * $y[$i];
	    
	    $sumx    += $x[$i];
	    $sumy    += $y[$i];
	}
    }
    my $ssx  = $sumxsqr - ($sumx * $sumx)/$n;
    my $ssy  = $sumysqr - ($sumy * $sumy)/$n;
    my $ssxy = $sumxy   - ($sumx * $sumy)/$n;
    my $rho = 0.0; 
    $rho = $ssxy / sqrt($ssx * $ssy) if (($ssx * $ssy) > 0.0); 
    
    return $rho;
    
}
######################################################################
sub mean {
    my ($x) = shift;
    my @x=@{$x};
    
    if ( scalar(@x) == 0){
	printf STDERR "WARNING: mean() array [@x] is empty!\n";
	return 0.0;
    }
    my ($n,$sum)=(0,0);
    while( my $xi = shift @x ) {
	$sum += $xi;
	$n++;
    }    
    return $sum/$n;
}

######################################################################
sub codon2aa {
    
    my (%geneticCode) = (
    'AAA' => 'K',    # Lysine
    'AAC' => 'N',    # Asparagine
    'AAG' => 'K',    # Lysine
    'AAT' => 'N',    # Asparagine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AGA' => 'R',    # Arginine
    'AGC' => 'S',    # Serine
    'AGG' => 'R',    # Arginine
    'AGT' => 'S',    # Serine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ATT' => 'I',    # Isoleucine
    'CAA' => 'Q',    # Glutamine
    'CAC' => 'H',    # Histidine
    'CAG' => 'Q',    # Glutamine
    'CAT' => 'H',    # Histidine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'GAA' => 'E',    # Glutamic Acid
    'GAC' => 'D',    # Aspartic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'TAA' => '*',    # Stop
    'TAC' => 'Y',    # Tyrosine
    'TAG' => '*',    # Stop
    'TAT' => 'Y',    # Tyrosine
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TGA' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGG' => 'W',    # Tryptophan
    'TGT' => 'C',    # Cysteine
    'TTA' => 'L',    # Leucine
    'TTC' => 'F',    # Phenylalanine
    'TTG' => 'L',    # Leucine
    'TTT' => 'F',    # Phenylalanine
    );
    
    my $codon = shift; 
    
    if (length($codon)!=3){
	print "WARNING: [$codon] should be length 3!\n";
	return $codon; 
    }
    
    if(defined($geneticCode{$codon})){
	return $geneticCode{$codon};
    }
    else {
	return 'X'; 
    }
}

