#!/usr/bin/perl 


use warnings;
use strict;
use Getopt::Long;


my ($infile, $isRNA, $help, $verbose);
&GetOptions (
        'h'         => \$help,
        'help'      => \$help,
        'r=s'       => \$isRNA,
        'i=s'       => \$infile,
        'v|verbose' => \$verbose
        );
if ($help)
{
print STDERR <<EOF;

 computeSynonNonsynon: 
A programme to count the number of synonymous and non-synonymous variation between two nucleotide sequences. 
input: Stockholm alignment of 2 RNA sequences (mRNA or ncRNA)
output: classification of variation between the sequences
     --code ncRNA: classify as basepair/loop neutral
     --code mRNA:  classify as synonymous/nonsynonymous variation 
     --biochemical ncRNA: R<->R and Y<->Y variation is neutral, otherwise not! 
     --biochemical mRNA:  B/A/N/P neutral variation is neutral, otherwise not!
     --biochemical mRNA:  BLOSUM62 (>=0) variation is neutral, otherwise not!
     --structural ncRNA:  fold with RNAfold -p, then convert matrix to probability of pairing for each nuc.
                          if value changes (more than 0.2[?]) then variation is not neutral  
     --structural mRNA:   fold protein with iTasser, if any probability of (alpha, beta or loop changes by more than 0.2) then variation is not neutral
     
Program to extend RNA stockholm format multiple alignments
Options:
        -h or --help    This help screen
        -r              Is this a ncRNA or not? 1 for yes, 0 for no.
	-i              A Stockholm alignment file, 2 sequences, both nucleotide
computeSequenceDistances.pl -r [0|1] -i [Stockholm alignment file]


   TODO:
print "|||+|||+|||||||||||||||||++++|-+||--|+||||++|-||||+||||||" strings for each synonymous type
implement structure 
     
     CONTROLS???????
--randomly introduce the same number of variation in: (e.g. SNPs, slide gap-blocks or regenerate alignments?) 


computeSynonNonsynon.pl -v -r 0 -i  mrna-seqs/PF00281.stk
computeSynonNonsynon.pl -v -r 1 -i ncrna-seqs/5S_rRNA.stk     
     
EOF
exit;
}

my (@names,@seqs);
my $sstruct='';
die "input file [$infile] doesn't exist! Use '-i <filename>'!" if (! -s $infile);
open(IN, "< $infile");
while(my $in = <IN>){
    chomp($in);
    #next if ($in =~ /^\#/);
    if($in =~ /^(\S+)\s+(\S+)$/){
	push(@names, $1); 
	push(@seqs,  $2); 
	#print "[$2]\n";
    }
    elsif($in =~ /^#=GC\s+SS_cons\s+(\S+)$/){
	$sstruct .= $1; 
	#print "[$1]\n";
    }
    
    
}
close(IN);
$seqs[0] =~ tr/a-z/A-Z/;
$seqs[1] =~ tr/a-z/A-Z/;

print "   seq1:[$seqs[0]] \n" if defined($verbose);
print "   seq2:[$seqs[1]] \n" if defined($verbose);
print "sstruct:[$sstruct] \n" if (defined($verbose) and length($sstruct));

#nucleotide distance
my $hDist = hammingDistance($seqs[0], $seqs[1]); 
printf "Total variation:\t%d\tlen=%d\n", $hDist, length($seqs[0]);

if(not defined($isRNA)){
    print "-r needs to be defined, you muppet!\n";
    exit(1);
}
elsif($isRNA){
    #RNA
    #if a ncRNA, encode in R/Y
    my $rySeq0 = nucleotide2ry($seqs[0]);
    my $rySeq1 = nucleotide2ry($seqs[1]);
    my $ryDist=hammingDistance($rySeq0, $rySeq1);
    #my $blsmDist=nuc2BLOSUMdist($rySeq0, $rySeq1);

    #take the difference between $hDist and $ryDist to get the synon:nonsynon counts:
    #printf "    ryDist:\t%d\n", $ryDist;
    printf "        rySynon:\t%d (%0.3f)\tryNonsynon:\t%d (%0.3f)\n", $hDist-$ryDist, ($hDist-$ryDist)/$hDist, $ryDist, $ryDist/$hDist;

    #structure distances:
    my $bpSyn = countStructurallySynonVariation($seqs[0], $seqs[1], $sstruct);
    printf "        bpSynon:\t%d (%0.3f)\tbpNonsynon:\t%d (%0.3f)\n", $bpSyn, $bpSyn/$hDist, $hDist-$bpSyn, ($hDist-$bpSyn)/$hDist;
    
    #printf "blosumDist:\t%0.2f\n", $blsmDist;    
#    foldingRNASynon($seqs[0], $seqs[1])
}
else {
    #PROTEIN
    #if a mRNA, translate, count synonymous w.r.t. AA seq variation:
    my $protSeq0 = translateMRNA($seqs[0]);
    my $protSeq1 = translateMRNA($seqs[1]);
    my $aaSyn = countAminoSynonVariation($seqs[0], $seqs[1], $protSeq0, $protSeq1);
    printf "        aaSynon:\t%d (%0.3f)\t     aaNonsynon:\t%d (%0.3f)\n", $aaSyn, $aaSyn/$hDist, $hDist-$aaSyn, ($hDist-$aaSyn)/$hDist;

    #encode in BANP alphabet, count synonymous w.r.t. BANP seq variation:
    my $banpSeq0 = protein2banp($protSeq0);
    my $banpSeq1 = protein2banp($protSeq1);
    my $banpSyn = countAminoSynonVariation($seqs[0], $seqs[1], $banpSeq0, $banpSeq1);
    printf "      banpSynon:\t%d (%0.3f)\t   banpNonsynon:   \t%d (%0.3f)\n", $banpSyn, $banpSyn/$hDist, $hDist-$banpSyn, ($hDist-$banpSyn)/$hDist;

    my $blossumSeq2 = protein2blossum62($protSeq0, $protSeq1);
    my $blossumSyn = countAminoSynonVariation($seqs[0], $seqs[1], $protSeq0, $blossumSeq2);
    printf "   blossumSynon:\t%d (%0.3f)\tblossumNonsynon:\t%d (%0.3f)\n", $blossumSyn, $blossumSyn/$hDist, $hDist-$blossumSyn, ($hDist-$blossumSyn)/$hDist;

    #foldingProteinSynon($seqs[0], $seqs[1])
}




exit(0); 

######################################################################
#compute the hamming distance between 2 sequences 
sub hammingDistance {
    
    my @seqs = @_; 
    
    my $len = length($seqs[0]);

    if ( ($len != length($seqs[0])) or (scalar(@seqs) != 2) ){
	print "ERROR: input seqs need to be the same length, and only 2 of them! MUPPET!\n"; 
	return -1;
    }
    
    my @seq1 = split(//, $seqs[0]);
    my @seq2 = split(//, $seqs[1]);
    my $dist = 0; 
    for (my $i=0; $i<$len; $i++){
	
	if($seq1[$i] ne $seq2[$i]){
	    $dist++;
	}
	
    }
    
    return $dist;
    
}

######################################################################
#fold RNA sequences, compute probabilities of paired/unpaired, count aligned characters that don't change  
sub foldingRNASynon {
    my ($seq1, $seq2) = @_; 
    
    my @seq1 = split(//, $seq1);
    my @seq2 = split(//, $seq2);
    
    my ($countAll, $countSynon) = (0,0);
    my ($seq1Ungap,$seq2Ungap)=('','');
    for (my $i=0; $i<length($seq1); $i++){

	if(isNucleotide($seq1[$i])){
	    $seq1Ungap .= $seq1[$i];
	}

	if(isNucleotide($seq2[$i])){
	    $seq2Ungap .= $seq2[$i];
	}
	
    }
    
    print "Ungapped seq1:[$seq1Ungap]\n" if defined($verbose);
    print "Ungapped seq2:[$seq2Ungap]\n" if defined($verbose);

    foldRNA($seq1Ungap);
    
    return 1; 
    
}

######################################################################
#fold RNA sequence
sub foldRNA {
    my $seq = shift;
    system("echo $seq | RNAfold -p");
    
    
    
}

######################################################################
#given 2 aligned protein sequences, return a sequence derived from seq2, with only changes for 0, or negative blosum62 scores   
sub protein2blossum62 {
    my ($seq1, $seq2) = @_; 

    my @blosum62 = qw(
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
	);
 #A R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  
    my @aas = qw(A R N D C Q E G H I L K M F P S T W Y V);

    my %blosum62; 
    for (my $i = 0; $i < scalar(@aas); $i++) {
        for (my $j = 0; $j < scalar(@aas); $j++) {
            $blosum62{$aas[$i]}{$aas[$j]} = $blosum62[$i * scalar(@aas) + $j];
        }
    }
    
    my @seq1 = split(//, $seq1); 
    my @seq2 = split(//, $seq2);
    my @blossumSeq = @seq2;
    print "blossum62:[" if defined($verbose);
    for (my $i=0; $i<length($seq1); $i++){

	if( ($seq1[$i] ne $seq2[$i]) and isAminoacid($seq1[$i]) and isAminoacid($seq2[$i]) and defined($blosum62{$seq1[$i]}{$seq2[$i]}) ){
	    if($blosum62{$seq1[$i]}{$seq2[$i]} > 0){
		$blossumSeq[$i]=$seq1[$i];
		print "+" if defined($verbose);
	    }
	    else{
		print "-" if defined($verbose);

	    }
	}
	else{
	    print "|" if defined($verbose);
	}
    }
    print "]\n" if defined($verbose);


    return join('', @blossumSeq);
}

######################################################################
#count the number of variable sites that are neutral w.r.t. amino acid sequence
sub countAminoSynonVariation{
    my ($seq1, $seq2, $aa1, $aa2) = @_; 
    
    my $len1 = length($seq1);

    if ($len1 != length($seq2)){
	print "ERROR: lengths of mRNA sequences are unequal, check Stockholm file!\n";
	print "  seq1:[$seq1]\n";
	print "  seq2:[$seq2]\n";
	print "   aa1:[$aa1]\n";
	print "   aa2:[$aa2]\n";
	return -1; 
    }
    
    
    my @seq1 = split(//, $seq1); 
    my @seq2 = split(//, $seq2); 
    my @aa1  = split(//, $aa1 ); 
    my @aa2  = split(//, $aa2 );
    my ($aaCnt, $countSynon) = (0,0);;
    for (my $i=0; $i<$len1; $i+=3){
	my $numVars = hammingDistance( $seq1[$i] . $seq1[$i+1] . $seq1[$i+2], $seq2[$i] . $seq2[$i+1] . $seq2[$i+2] ) if(defined($seq1[$i+2]) and defined($seq2[$i+2]));
	if($aa1[$aaCnt] eq $aa2[$aaCnt]){
	    $countSynon += $numVars;
	}
	$aaCnt++;
    }

    return $countSynon;

}

######################################################################
#count the number of variable sites that are structurally neutral:
sub countStructurallySynonVariation{
    my ($seq1, $seq2, $struct) = @_; 
    
    if (length($seq1) != length($seq2) or length($seq1) != length($struct)){
	print "ERROR: lengths of sequences and structure strings are unequal, check Stockholm file!\n";
	print "  seq1:[$seq1]\n";
	print "  seq2:[$seq2]\n";
	print "struct:[$struct]\n";
	return -1; 
    }

    my @seq1 = split(//, $seq1);
    my @seq2 = split(//, $seq2);
    my @struct = split(//, $struct);
    my @bpTable = @{ make_pair_table($struct) };
    #my @bpTable = @{$bpTable};
    my ($countAll, $countSynon) = (0,0);
    for (my $i=0; $i<length($seq1); $i++){

	if($seq1[$i] ne $seq2[$i]){
	    $countAll++;

	    if (not $bpTable[$i]){#unpaired nucleotide 
		$countSynon++;
	    }
	    elsif( isComplementary( $seq1[$i], $seq1[ $bpTable[$i]-1 ]) && isComplementary( $seq2[$i], $seq2[ $bpTable[$i]-1 ]) ){#covarying
		$countSynon++;
	    }
	    elsif( (isGap($seq1[$i]) && isGap($seq1[ $bpTable[$i]-1 ]) ) or (isGap($seq2[$i]) && isGap($seq2[ $bpTable[$i]-1 ]) ) ){#gap-gap -- count as synonymous
		$countSynon++;
	    }
	}
	
    }

    return $countSynon;
}

######################################################################
#make_pair_table: takes a structure string and returns an array/pair_table.
#                 pair_table[0] contains the length of the string.
#                 pair_table[1..pair_table[0]] contains 0's if no basepair exists 
#                 for this position, otherwise it contains the index for the 
#                 corresponding pair.
#                 Eg. make_pair_table(".((...))") returns:
#                 (8,0,8,7,0,0,0,3,2)
######################################################################
sub make_pair_table {

    my $str = shift;
        
    my (%bpsymbs_posns, @pair_table, $prime5, $prime3, %bpsymbs3p5p, %bpsymbs3p5p_counts);
    my ($count, $unbalanced) = (0,0);

    #Match up the 5' and 3' basepair simbols:
    my %bpsymbs5p3p = (
	'(' => ')',
	'<' => '>',
	'[' => ']',
	'{' => '}',
	A => 'a',
	B => 'b',
	C => 'c',
	D => 'd',
	E => 'e',
	F => 'f',
	G => 'g',
	H => 'h',
	I => 'i'
	);

    my %bpsymbs5p3p_counts = (
	'(' => 0,
	'<' => 0,
	'[' => 0,
	'{' => 0,
	A => 0,
	B => 0,
	C => 0,
	D => 0,
	E => 0,
	F => 0,
	G => 0,
	H => 0,
	I => 0
	);
    
    #We also need the reverse of the above hashes (ie. Match the 3' symbols with the 5'):
    foreach my $symb5p (keys %bpsymbs5p3p){
	my $symb3p = $bpsymbs5p3p{$symb5p};
	$bpsymbs3p5p{$symb3p} = $symb5p;
	$bpsymbs3p5p_counts{$symb3p} = 0;
    }
    
    my %unpairedsymbs = (
	'.' => 1,
	',' => 2,
	'-' => 3,
	':' => 4,
	'_' => 5,
	'~' => 6
	);
    
    my @ss = split(//,$str);
    $pair_table[0]  = length($str);
    
    my $j=0;
    foreach my $char (@ss){
	$j++;
	$pair_table[$j] = 0;
	
	if ( defined( $unpairedsymbs{$char} ) ) {
	    next; #Boring unpaired region.
	}
	elsif ( defined( $bpsymbs5p3p{$char}) ){#Record position of open bps:
	    push( @{ $bpsymbs_posns{$char} }, $j);
	    ++$count;
	    $bpsymbs5p3p_counts{$char}++;
	}
	elsif ( defined( $bpsymbs3p5p{$char}) ){#close bp, save positions of matches in pair_table:
	    my $mchar = $bpsymbs3p5p{$char};
	    $prime5 = pop( @{ $bpsymbs_posns{$mchar} } );
	    $prime3 = $j;
	    $pair_table[$prime3] = $prime5;
	    $pair_table[$prime5] = $prime3;
	    --$count;
	    $bpsymbs3p5p_counts{$char}++;
	}
       	else {
	    printf STDERR "Strange character \"$char\" in secondary structure:\n$str\n";
	}
    }
    
    #Check basepair symbols are all matched:
    foreach my $symb5p (keys %bpsymbs5p3p_counts){
	my $symb3p = $bpsymbs5p3p{$symb5p};
	my $diff = $bpsymbs5p3p_counts{$symb5p} - $bpsymbs3p5p_counts{$symb3p};
	if ($diff!=0){
	    printf STDERR "BAD EVIL AND NASTY: Unbalanced brackets in secondary structure:\n$bpsymbs5p3p_counts{$symb5p}x\'$symb5p\' and $bpsymbs3p5p_counts{$symb3p}x\'$symb3p\'\n";
	    $unbalanced = 1;
	}
    }
    
    
    if ($count != 0 || $unbalanced){
	printf STDERR "Unbalanced basepair symbols in secondary structure:\n$str\n";
	return ();
    }    
    else {
	return \@pair_table;
    }
    
}

    
######################################################################
#returns true if the two input characters can form a canonical/watson-crick basepair
sub isComplementary {
    my $a = shift;
    my $b = shift;

#    if (defined($a)){
#	$a =~ tr/a-z/A-Z/;
#    }

#    if (defined($b)){
#	$b =~ tr/a-z/A-Z/;
#    }
    
    my $ab = $a . $b;
    
    my %canonicalBasepair = (
	AU => 1,
	AT => 1,
	UA => 1,
	TA => 1,
	CG => 1,
	GC => 1,
	UG => 1,
	GU => 1,
	TG => 1,
	GT => 1,
	RY => 1,
	YR => 1,
	MK => 1,
	KM => 1,
	SS => 1,
	WW => 1
	);
    
    if ( defined($canonicalBasepair{$ab} ) ) {
	return 1;
    }
    else {
        return 0;
    }
	 
}

######################################################################
#convert to R/Y
sub nucleotide2ry {
    
    my $seq = shift; 
    
    my $len = length($seq);
    my $rySeq='';
    my @seq = split(//, $seq); 
    for (my $i=0; $i<$len; $i++){
	$rySeq .= nuc2RY($seq[$i] )
	
    }    
    print "  ryseq:[$rySeq] len:[$len]\n" if defined($verbose);
    
    return $rySeq;

}

######################################################################
#nucleotide to R/Y
sub nuc2RY {
    
    my $nuc = shift;
    if (length($nuc)!=1){
	print "WARNING: [$nuc] should be a single character!\n";
	return $nuc; 
    }

    if ($nuc =~ /[AG]/){
	return 'R'; 
    }
    elsif ($nuc =~ /[TUC]/){
	return 'Y'; 
    }
    else {
	return $nuc;
    }
    
}


######################################################################
#nucleotide BLOSUM score:
# sub nuc2BLOSUMdist {
# #+5/-4 match, mis-match scores 

# #D(x, y) = S(x, x) − S(x, y)   
# #Baussand and Carbone (2008) Inconsistent Distances in Substitution Matrices can be Avoided by Properly Handling Hydrophobic Residues

    
#     my @seqs = @_; 
    
#     my $len = length($seqs[0]);
    
#     if ( ($len != length($seqs[0])) or (scalar(@seqs) != 2) ){
# 	print "ERROR: input seqs need to be the same length, and only 2 of them! MUPPET!\n"; 
# 	return -1;
#     }
    
#     my @seq1 = split(//, $seqs[0]);
#     my @seq2 = split(//, $seqs[1]);
#     my ($Sxx, $Sxy) = (0,0); 
#     for (my $i=0; $i<$len; $i++){
	
# 	if(isNucleotide($seq1[$i]) and isNucleotide($seq2[$i])){
# 	    if($seq1[$i] ne $seq2[$i]){
# 		$Sxy-=4;
# 	    }
# 	    else {
# 		$Sxy+=5;
# 	    }
# 	    $Sxx+=5;
# 	}
	
#     }
    
#     return ($Sxx-$Sxy)/$len;
    
# }


######################################################################

sub translateMRNA {

    my $nucSeq = shift; 
    
    my $len = length($nucSeq);
    my $protSeq='';
    
    if ( ($len % 3) ){
	print "WARNING: length of [$nucSeq] is not a factor 3! Frameshifts may cause errors!\n";
    }
    
    my @nucSeq = split(//, $nucSeq); 
    for (my $i=0; $i<$len; $i+=3){
	$protSeq .= codon2aa( $nucSeq[$i] . $nucSeq[$i+1] . $nucSeq[$i+2] ) if(defined($nucSeq[$i+2]));
    }
    
    print "  protseq:[$protSeq]\n" if defined($verbose);
    
    return $protSeq;
    
}

sub codon2aa {
    
    my (%geneticCode) = (
    'AAA' => 'K',    # Lysine
    'AAA' => 'K',    # Lysine
    'AAC' => 'N',    # Asparagine
    'AAC' => 'N',    # Asparagine
    'AAG' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AAT' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'ACA' => 'T',    # Threonine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AGA' => 'R',    # Arginine
    'AGA' => 'R',    # Arginine
    'AGC' => 'S',    # Serine
    'AGC' => 'S',    # Serine
    'AGG' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'AGT' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'ATA' => 'I',    # Isoleucine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ATG' => 'M',    # Methionine
    'ATT' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'CAA' => 'Q',    # Glutamine
    'CAA' => 'Q',    # Glutamine
    'CAC' => 'H',    # Histidine
    'CAC' => 'H',    # Histidine
    'CAG' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CAT' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CCA' => 'P',    # Proline
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CGA' => 'R',    # Arginine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'CTA' => 'L',    # Leucine
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'GAA' => 'E',    # Glutamic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAC' => 'D',    # Aspartic Acid
    'GAC' => 'D',    # Aspartic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GCA' => 'A',    # Alanine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GGA' => 'G',    # Glycine
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    'GTA' => 'V',    # Valine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'TAA' => '_',    # Stop
    'TAA' => '_',    # Stop
    'TAC' => 'Y',    # Tyrosine
    'TAC' => 'Y',    # Tyrosine
    'TAG' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TAT' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TCA' => 'S',    # Serine
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TGA' => '_',    # Stop
    'TGA' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGC' => 'C',    # Cysteine
    'TGG' => 'W',    # Tryptophan
    'TGG' => 'W',    # Tryptophan
    'TGT' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TTA' => 'L',    # Leucine
    'TTA' => 'L',    # Leucine
    'TTC' => 'F',    # Phenylalanine
    'TTC' => 'F',    # Phenylalanine
    'TTG' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TTT' => 'F',    # Phenylalanine
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

sub protein2banp {

    my $protSeq = shift; 
    my $len = length($protSeq);
    my $banpSeq='';

    my @protSeq = split(//, $protSeq); 

    for (my $i=0; $i<$len; $i++){
	$banpSeq .= aa2banp( $protSeq[$i] );
    }
    print "  banpseq:[$banpSeq] len:[$len]\n" if defined($verbose);
    
    return $banpSeq;
}

sub aa2banp {
    
    my (%aa2banp) = (
	'H' => 'B',
	'R' => 'B',
	'K' => 'B',
	'D' => 'A',
	'E' => 'A',
	'G' => 'P',
	'S' => 'P',
	'Y' => 'P',
	'C' => 'P',
	'T' => 'P',
	'N' => 'P',
	'Q' => 'P',	
	'F' => 'N',
	'L' => 'N',
	'W' => 'N',
	'P' => 'N',
	'I' => 'N',
	'M' => 'N',
	'V' => 'N',
	'A' => 'N'
       );
    
    my $aa = shift; 

    if(defined($aa2banp{$aa})){
	return $aa2banp{$aa};
    }
    else {
	return $aa; 
    }
    
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
#returns true if input character is an amino acid (IUPAC codes):
sub isAminoacid {
    my $a = shift;
    
    if (defined($a)){
	$a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ARNDCQEGHILKMFPSTWYV]/) ){
	return 1;
    }
    else {
	return 0;
    }
    
}



######################################################################
#returns true if the input character is a gap character:
sub isGap {
    my $a = shift;

    my %gap = (
	'.' => 1,
	',' => 2,
	'-' => 3,
	':' => 4,
	'_' => 5,
	'~' => 6
	);

    if( defined($gap{$a}) ){
	return 1;
    }
    else {
	return 0;
    }
}

