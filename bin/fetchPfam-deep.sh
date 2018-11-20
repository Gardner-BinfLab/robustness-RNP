#!/bin/sh

echo "grep \"ID=$1;\" U00096-pfam31.gff   | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42\"\'   | sh  > mrna-seqs/e1.fasta"
echo "grep \"ID=$1;\" AL157959-pfam31.gff | sort -k6nr | head -n 1 | perl -lane \'($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq \"-\");  print \"sfetch -d AL157959.fasta -r \42AL157959/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AL157959;\42\"\' | sh >> mrna-seqs/n1.fasta"
grep "ID=$1;" U00096-pfam31.gff   | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42"'   | sh  > mrna-seqs/e1.fasta
grep "ID=$1;" AL157959-pfam31.gff | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d AL157959.fasta -r \42AL157959/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AL157959;\42"' | sh >> mrna-seqs/n1.fasta

echo "cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta"
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta > mrna-seqs/$1.fasta

~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-stranslate -r 50 mrna-seqs/e1.fasta | perl -lane 'if(/^(>\S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  > mrna-seqs/$1.pep
~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-stranslate -r 50 mrna-seqs/n1.fasta | perl -lane 'if(/^(>\S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}' >> mrna-seqs/$1.pep
#echo "/home/ppg15/inst/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-translate -l 50 -c 11 --watson mrna-seqs/e1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1..\d+ length=\d+ frame=1/){print \">$1\"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  > mrna-seqs/$1.pep"
#echo "/home/ppg15/inst/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-translate -l 50 -c 11 --watson mrna-seqs/n1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1..\d+ length=\d+ frame=1/){print \">$1\"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}' >> mrna-seqs/$1.pep"
#/home/ppg15/inst/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-translate -l 50 -c 11 --watson mrna-seqs/e1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1..\d+ length=\d+ frame=1/){print ">$1"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  > mrna-seqs/$1.pep
#/home/ppg15/inst/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-translate -l 50 -c 11 --watson mrna-seqs/n1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1..\d+ length=\d+ frame=1/){print ">$1"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}' >> mrna-seqs/$1.pep

cat mrna-seqs/$1.pep

rm mrna-seqs/e1.fasta mrna-seqs/n1.fasta

hmmfetch ~/data/pfam/pfam31/Pfam-A.hmm  $1               > $1.hmm
hmmalign                        $1.hmm mrna-seqs/$1.pep  > mrna-seqs/$1.pep.stk
cat mrna-seqs/$1.pep.stk

~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat --mingap clustal mrna-seqs/$1.pep.stk | tr "." "-" > mrna-seqs/$1.pep.aln 
cat mrna-seqs/$1.pep.aln

echo "pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa"
pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa

#FAILS, FUCKING USELESS:
#echo "~/inst/RevTrans-1.4/revtrans.py -match name -vvv mrna-seqs/$1.fasta mrna-seqs/$1.pep.aln -O aln > mrna-seqs/$1.nuc.aln"
#~/inst/RevTrans-1.4/revtrans.py -match name -vvv mrna-seqs/$1.fasta mrna-seqs/$1.pep.aln -O aln > mrna-seqs/$1.nuc.aln

#~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk
sreformat -u --informat a2m pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk

head -n 50 mrna-seqs/$1.nuc.stk
