#!/bin/sh

#Fetch annotations of selected Pfam domain from GFFs:
if [ "$2" ]; then
    
    echo "cat   U00096-ncRNAs-nonredundant.gff   U00096-pfam31.gff  > blah && sortGffs.pl blah | grep -C 15 $2 | grep \"ID=$1;\" | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42\"\'   | sh  > mrna-seqs/e1.fasta"
    echo "cat AE014613-ncRNAs-nonredundant.gff AE014613-pfam31.gff  > blah && sortGffs.pl blah | grep -C 15 $2 | grep \"ID=$1;\" | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d AE014613.fasta -r \42AE014613/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AE014613;\42\"\' | sh >> mrna-seqs/n1.fasta"
    
    cat   U00096-ncRNAs-nonredundant.gff   U00096-pfam31.gff  > blah && sortGffs.pl blah | grep -C 15 $2 | grep "ID=$1;"  | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42"'   | sh  > mrna-seqs/e1.fasta
    cat AE014613-ncRNAs-nonredundant.gff AE014613-pfam31.gff  > blah && sortGffs.pl blah | grep -C 15 $2 | grep "ID=$1;" | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d AE014613.fasta -r \42AE014613/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AE014613;\42"' | sh >> mrna-seqs/n1.fasta
    
else 
    
    echo "grep \"ID=$1;\" U00096-pfam31.gff   | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42\"\'   | sh  > mrna-seqs/e1.fasta"
    echo "grep \"ID=$1;\" AE014613-pfam31.gff | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d AE014613.fasta -r \42AE014613/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AE014613;\42\"\' | sh >> mrna-seqs/n1.fasta"
    grep "ID=$1;" U00096-pfam31.gff   | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42"'   | sh  > mrna-seqs/e1.fasta
    grep "ID=$1;" AE014613-pfam31.gff | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d AE014613.fasta -r \42AE014613/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AE014613;\42"' | sh >> mrna-seqs/n1.fasta
    
fi

#Fetch nucleotide sequences:
echo "cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta"
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta > mrna-seqs/$1.fasta

#Translate them (check translations are correct):
~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-stranslate -r 10 mrna-seqs/e1.fasta | perl -lane 'if(/^(>\S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  > mrna-seqs/$1.pep
~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-stranslate -r 10 mrna-seqs/n1.fasta | perl -lane 'if(/^(>\S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}' >> mrna-seqs/$1.pep

cat mrna-seqs/$1.pep
rm mrna-seqs/e1.fasta mrna-seqs/n1.fasta

#Build protein sequence alignment (could use MAFFT/MUSCLE/CLUSTAL instead)...
hmmfetch ~/data/pfam/pfam31/Pfam-A.hmm  $1               > $1.hmm
hmmalign                        $1.hmm mrna-seqs/$1.pep  > mrna-seqs/$1.pep.stk
cat mrna-seqs/$1.pep.stk

#Convert to Clustal format:
~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat --mingap clustal mrna-seqs/$1.pep.stk | tr "." "-" > mrna-seqs/$1.pep.aln 
cat mrna-seqs/$1.pep.aln

#Reverse translate alignment (i.e. align mRNA sequences consistently with the protein sequence alignments):
echo "pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa"
pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa

#FAILS, FUCKING USELESS: ~/inst/RevTrans-1.4/revtrans.py

#DONE!!!
sreformat -u --informat a2m pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk
head -n 50 mrna-seqs/$1.nuc.stk

######################################################################




