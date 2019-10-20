#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;
use Statistics::RankCorrelation;
use File::Copy;
use Data::Dumper;
    
my ($ncRNAFile, $mRNAFile, $proteinFile, $pdb2model, $pdb2cmscore, $pdb2hmmscore) = qw(ncRNA-pdb-seqs.fa  mRNA-pdb-seqs.fa  protein-pdb-seqs.fa hmm-cm-models/model2pdb-maps.txt hmm-cm-models/ncRNA-pdb-seqs.rfam14.cmscan.tblout hmm-cm-models/protein-pdb-seqs.hmmscan.tblout);
my $numMutations     = 1;   #number of mutations per sequence 
my $numMutants       = 100;  #number of mutant sequences to generate per sequence 
my ($help, $verbose) = (0,0);

&GetOptions (
    'h'              => \$help,
    'help'           => \$help,
    'n=s'            => \$ncRNAFile,
    'm=s'            => \$mRNAFile, 
    'p=s'            => \$proteinFile,
    'numMutants=s'   => \$numMutants,
    'numMutations=s' => \$numMutations,
    'v|verbose'      => \$verbose
        );
if ($help)
{
print STDERR <<EOF;

 computeDeltaDelta: 
A programme to assess the impact of randomnly generate mutations on RNA/protein structures. 

input: 
     Sequence files: [DEFAULT filenames]
              ncRNA-pdb-seqs.fa  
              mRNA-pdb-seqs.fa  
              protein-pdb-seqs.fa 
     
     Mapping between PDB structures and HMM/CM models: 
                          hmm-cm-models/model2pdb-maps.txt 
     
     Precomputed HMM/CM scores for unmutated sequences: 
                          hmm-cm-models/ncRNA-pdb-seqs.rfam14.cmscan.tblout 
                          hmm-cm-models/protein-pdb-seqs.hmmscan.tblout

     PDB files:
              in a directory called 'pdbfiles/'
output: 
     TAB delimited output of scores for each simulated mutant:
                          ID     - an identifier for each mutant  	
                          ddG	 - estimated delta delta G value 
                          dBS	 - delta bitscore 
                          type	 - rna or protein
                          mutSeq - the mutated sequence
     
Program to extend RNA stockholm format multiple alignments
Options:
     -h or --help    This help screen
     -v              print verbose output to STDERR
     -n              use a different ncRNA sequence file
     -m              use a different mRNA sequence file
     -p              use a different protein sequence file
     -numMutants     how many mutant sequences to generate per sequence [DEFAULT: 100] 
     -numMutations   how many mutations to introduce into each sequence [DEFAULT: 1] 
   TODO:

     computeDeltaDelta.pl -v
     computeDeltaDelta.pl -v  -n 6nd4_2.fa  -m 6nd4_c-mRNA.fa   -p  6nd4_2-protein.fa     
     
     
EOF
exit;
}



######################################################################
#Read in model/PDB pairs & precomputed datasets:
my (%pdb2model, %wildTypeBitscore); 
if (-s $pdb2model && -s $pdb2cmscore && -s $pdb2hmmscore){
    open(MIN, "< $pdb2model") or die "FATAL: failed to open [$pdb2model]!\n [$!]";
    while (my $min = <MIN>){
        if($min =~ /^(\S+)\t(\S+)/){
            $pdb2model{$2} = $1;
            #print "$1\t$2\n";
        }
    }
    close(MIN);

    open(NIN, "< $pdb2cmscore") or die "FATAL: failed to open [$pdb2cmscore]!\n [$!]";
    while (my $nin = <NIN>){
        next if $nin =~ /^#/;
        my @nin = split(/\s+/, $nin); 
        if( defined($nin[3]) && defined( $pdb2model{ $nin[3] } ) &&  $pdb2model{$nin[3]} eq $nin[1]){
            #print "$nin[3]\t$nin[16]\n";
            $wildTypeBitscore{$nin[3]} = $nin[16]; 
        }
    }
    close(NIN);

    open(PIN, "< $pdb2hmmscore") or die "FATAL: failed to open [$pdb2hmmscore]!\n [$!]";
    while (my $pin = <PIN>){
        next if $pin =~ /^#/;
        my @pin = split(/\s+/, $pin); 
        if( defined($pin[2]) && defined( $pdb2model{$pin[2]} ) &&  $pdb2model{$pin[2]} eq $pin[0]){
            $wildTypeBitscore{$pin[2]} = $pin[5];
        }
    }
    close(PIN);
}
else {
    print STDERR "WARNING: missing [$pdb2model], [$pdb2cmscore] or [$pdb2hmmscore] input files!\n"; 
}

######################################################################
#Mutate RNAs 
print "ID\tddG\tdBS\ttype\tmutSeq\n"; 
if (-s $ncRNAFile){

    my ($ncRNAptr, $offsetN) = readFasta($ncRNAFile);
    my %ncRNAseqs = %{ $ncRNAptr };
    
    foreach my $n (sort keys %ncRNAseqs){
        
        #print "n[$n]\tseq[$ncRNAseqs{$n}]\tmodel[$pdb2model{$n}]\n";

        my $WTmfe = foldRNA($ncRNAseqs{$n});
        for (my $i=0; $i<$numMutants; $i++ ){
            my $mutSeq     = mutateSequence($ncRNAseqs{$n}, $numMutations);	
            my $MTmfe = foldRNA($mutSeq);
	    my $ddG = 'NA';
	    $ddG    = sprintf("%0.2f", $MTmfe-$WTmfe) if ($MTmfe !~ /NA/ && $WTmfe !~ /NA/);
            my $dbs = computeDeltaCmBS( $mutSeq, $pdb2model{$n}, $wildTypeBitscore{$n} );
            
            printf "$n\.mut$i\t$ddG\t$dbs\trna\t$mutSeq\n";
        }
        
    }
    
}
else {
    print STDERR "WARNING: missing [$ncRNAFile] input files!\n"; 
}

######################################################################
#Mutate proteins 
if (-s $mRNAFile && -s $proteinFile){

    my ($mRNAptr, $offsetM) = readFasta($mRNAFile);
    my ($protptr, $offsetP) = readFasta($proteinFile);
    my %mRNAseqs = %{ $mRNAptr };
    my %protseqs = %{ $protptr };
    my %offsetsPDB  = %{ $offsetP  };
    
    foreach my $n (sort keys %mRNAseqs){
        my ($s,$e, $pdbId, $chain);
        my $name = $n;
        if($n =~ /(\S+)\/(\d+)-(\d+)/){
            ($name,$s,$e)=($1,$2,$3);
        }
        $offsetsPDB{$name} = 0 if not defined($offsetsPDB{$name}); 
        ($pdbId, $chain) = split(/_/, $name);
        
        if (defined($protseqs{$name})){
            #printf ">$name";
            my $prot        = $protseqs{$name};
            my $translation = translateMRNA($mRNAseqs{$n});
            if (defined ($s)){ #if mRNA is a partial match to the 3D structure:
                #printf "\t$s\t$e";
                $prot = substr($prot, $s-1, ($e - $s +1)); 
            }
            else {
                $s=1; 
            }
            
            #check mRNA and protein sequences are consistent:
            if($prot ne $translation && length($prot) == length($translation) ){
                print STDERR "WARNING: [$name] protein seq and translation of [$n] have point mutations!\n";
            }
            elsif ($prot ne $translation && length($prot) != length($translation) ){
                print STDERR "ERROR: [$name] protein seq and translation of [$n] have different lengths!\n";
            }
            
            #printf "\n$mRNAseqs{$n}\n%s\n%s\n", $translation, $prot;
            #generate mutant mRNAs & corresponding proteins here: 
            for (my $i=0; $i<$numMutants; $i++ ){
                my $mutSeq     = mutateSequence($mRNAseqs{$n}, $numMutations);	
                my $mutProtSeq = translateMRNA($mutSeq);
                #return variant sites in maestro format: "<wild type amino acid><residue id[icode]>[.chains]{<substitutions>}" 
                #e.g. I102.A{L}
                my ($numSites, $maestroString) = seqDiffs($mutProtSeq, $prot, $chain, $s-1, $offsetsPDB{$name});
                if ($numSites){
                    #print ">mut$i (prot)\n$mutProtSeq\n$maestroString\n";
                    my $printString = computeDeltaDeltaG($pdbId, $chain, $i, $maestroString ); 
                    if (not defined($printString) or length($printString) < 10){ 
                        print STDERR "WARNING: computeDeltaDeltaG failed to return a result!\n";
                        next;
                    }
                    ######################################################################
                    $printString .= computeDeltaHmmBS( $mutProtSeq, $pdb2model{$name}, $wildTypeBitscore{$name} );
                    
                    ######################################################################
                    print $printString . "\tprotein\t$mutProtSeq\n";
                }
                else {
                    print "$name\.mut$i\t0.000\t0.0\tprotein\t$mutProtSeq\n";
                }
            }
            
            
            
        }
        else {
            print STDERR "ERROR: missing protein sequence for [$name]  [$n]!\n";
        }
        
    }
}
else {
    print STDERR "WARNING: missing [$mRNAFile] &/or [$proteinFile] input files!\n"; 
}

#print "\n\n" . ("#" x 70) . "\n\n";

#Fin.
exit(0);

######################################################################

sub readFasta {
    
    my $file = shift;
    my %seq;
    my $name;
    my %offsets; 
    open(IN, "< $file") or die "FATAL: failed to open [$file] in readFasta!\n [$!]";
    while (my $in = <IN>){
        if($in =~ /^>\s*(\S+)/){
            $name = $1;
            print STDERR "WARNING: duplicate sequence name [$name] in [$file]!" if (defined($seq{$name})); 
            $seq{$name}='';

            if($in =~ /offset=(\S+)/){
                $offsets{$name}=$1; 
            }
        }
        elsif($in =~ /(\S+)/){
            $in =~ s/\s+//g;
            chomp($in);
            $seq{$name} .= $in;
        }
    }
    close(IN);
    return (\%seq, \%offsets); 
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
	print STDERR "WARNING: [$codon] should be length 3!\n";
	return $codon; 
    }
    
    if(defined($geneticCode{$codon})){
	return $geneticCode{$codon};
    }
    else {
	return 'X'; 
    }
}


###REWRITE THIS!!!!!!
######################################################################
#Given a sequence, generate a mutant with a given number of variable sites 
sub mutateSequence{
    my ($seq, $numMuts) = @_;
    
    my @seq = split(//, $seq);

    for (my $i=0; $i<$numMuts; $i++ ){
        
        my $pos = rand( length($seq) );
	#Point mutation
	if(isNucleotide($seq[$pos]) ){
            $seq[$pos] = mutateChar($seq[$i]);
	}
    }
    my $mutSeq=join('', @seq);

    return ($mutSeq);
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
#seqDiffs: compare 2 protein strings, return differences in maestro format
sub seqDiffs {
    #return variant sites in maestro format: "<wild type amino acid><residue id[icode]>[.chains]{<substitutions>}" 
    #e.g. I102.A{L}
    my ($mutProtSeq, $prot, $chain, $offset, $offsetPdb) = @_;
    my ($numSites, $maestroString)=(0,'');
    my @maestroStrings;
    my @mutProtSeq = split(//, $mutProtSeq); 
    my @prot       = split(//, $prot); 
    for (my $i=1; $i<scalar(@prot)+1; $i++){
        if(defined($mutProtSeq[$i-1]) && $mutProtSeq[$i-1] ne $prot[$i-1]){
            my $pos = $i + $offset + $offsetPdb;
            next if $mutProtSeq[$i-1] eq 'X';
            push(@maestroStrings, "$prot[$i-1]$pos\.$chain\{$mutProtSeq[$i-1]\}");
            $numSites++;
        }
        elsif( (not defined($mutProtSeq[$i-1])) && ($prot[$i-1] !~ /A/i) ){
            my $pos = $i + $offset + $offsetPdb;
            push(@maestroStrings, "$prot[$i-1]$pos\.$chain\{A\}");
            $numSites++;
        }
    }
    $maestroString = join(',', @maestroStrings);
    return($numSites, $maestroString);
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
#
sub computeDeltaDeltaG {
    my ($pdbId, $chain, $i, $maestroString) = @_; 
    
    my $runMe = "\$MAESTRO_HOME/maestro \$MAESTRO_HOME/config.xml ./pdbfiles/$pdbId\.pdb --evalmut=\47$maestroString\47";
    my $printMe = "$pdbId\_$chain\.mut$i\tNA";; 
    open(MAESTRO, "$runMe |")  or die "FATAL: failed to open pipe: [$runMe]\n[$!]";
    while (my $ms = <MAESTRO>){
        next if ($ms =~ /structure|wildtype/);
        my @ms = split(/\s+/, $ms);
        if( isNumeric($ms[5]) ){
            $printMe = "$pdbId\_$chain\.mut$i\t$ms[5]";
        }
    }
    close(MAESTRO);
    return $printMe; 
}

######################################################################
#
sub computeDeltaHmmBS {

    my ($mutProtSeq, $model, $wildTypeBitscore ) = @_;

    open(UT, "> seq");
    print UT ">mutSeq\n$mutProtSeq\n";
    close(UT);
    
    my $runMe = "hmmscan -T 5 --tblout tab hmm-cm-models/$model" . "_hmm.txt seq > /dev/null";
    system( "$runMe" ) and print STDERR "FAILED: failed to run [$runMe]!\n[$!]";
    my $printMe = "\tNA";
    open(TB, "grep $model tab |") or print STDERR  "FATAL: failed to open pipe: [grep $model tab]\n[$!]";
    while (my $tb = <TB>){
        my @tb = split(/\s+/, $tb);
        if( defined($tb[5]) && isNumeric($tb[5]) ){
            $printMe = "\t" . sprintf("%0.2f", $tb[5] - $wildTypeBitscore);
            last;
        }
    }
    close(TB);
    
    print STDERR "Warning: computeDeltaHmmBS failed! mutSeq[$mutProtSeq], HMM[$model], WTbitscore[$wildTypeBitscore]\n" if $printMe =~ /NA/;
    return $printMe;
}
######################################################################
#
sub computeDeltaCmBS { 
    
    my ($mutSeq, $model, $wildTypeBitscore ) = @_;
    
    open(UT, "> seq");
    print UT ">mutSeq\n$mutSeq\n";
    close(UT);
    
    my $runMe = "cmscan -T 5 --tblout tab hmm-cm-models/$model" . ".cm seq > /dev/null";
    #print "$runMe\n";
    system( "$runMe" ) and print STDERR "WARNING: failed to run [$runMe]!\nmutSeq[$mutSeq], CM[$model], WTbitscore[$wildTypeBitscore]\n[$!]";
    my $printMe = "NA";
    #print "grep $model tab\n";
    open(TB, "grep $model tab |")  or die "FATAL: failed to open pipe: [grep $model tab]\n[$!]";
    while (my $tb = <TB>){
        my @tb = split(/\s+/, $tb);
        if( isNumeric($tb[14]) ){
            $printMe = sprintf("%0.2f", $tb[14] - $wildTypeBitscore);
            last;
        }
    }
    close(TB);
    return $printMe;
    print STDERR "$printMe\n";
    
}

######################################################################
#
sub foldRNA{
    my $seq = shift; 

    system("echo $seq | RNAfold > RNAfold-p.out 2> RNAfold-p.err") and die "FATAL: failed to run system command: [echo $seq | RNAfold -p > RNAfold-p.out 2> RNAfold-p.err]\n[$!]";

    my $mfe='NA';
    open(RF, "echo $seq | RNAfold 2> RNAfold-p.err |")  or die "FATAL: failed to open pipe: [RNAfold]\n[$!]";
    while (my $rf = <RF>){
        if( $rf =~ /^\S+\s+\((\S+)\)/){
            $mfe = $1; 
        }
    }
    close(RF);
    return $mfe;
}
