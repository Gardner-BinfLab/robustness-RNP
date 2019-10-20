#!/bin/sh

#Dependencies:
#**esl-sfetch
#**esl-translate
#**esl-reformat
#**pal2nal.pl


ID=$1
GFFFILE1=$2
SEQFILE1=$3
GFFFILE2=$4
SEQFILE2=$5

#return error if inputs don't exist

mkdir -p mrna-seqs
esl-sfetch --index $SEQFILE1
esl-sfetch --index $SEQFILE2

#fetch highest scoring match
grep        "ID=$ID;"  $GFFFILE1 | sort -k6nr | head -n 1 | perl -lane  '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq  "-" );  print  "esl-sfetch  -n \42$F[0]/$F[3]-$F[4]\42 -c $F[3]-$F[4] '"$SEQFILE1"' \42$F[0]\42"' | sh  > mrna-seqs/e1.fasta || { echo "sfetch failed for failed $SEQFILE1" ; exit 1; }
grep        "ID=$ID;"  $GFFFILE2 | sort -k6nr | head -n 1 | perl -lane  '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq  "-" );  print  "esl-sfetch  -n \42$F[0]/$F[3]-$F[4]\42 -c $F[3]-$F[4] '"$SEQFILE2"' \42$F[0]\42"' | sh  > mrna-seqs/n1.fasta || { echo "sfetch failed for failed $SEQFILE2" ; exit 1; }


echo "## mRNA sequences: [cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta]"
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta
cat mrna-seqs/e1.fasta mrna-seqs/n1.fasta > mrna-seqs/"$ID".fasta

#/home/paulgardner/inst/squid-1.9g/translate -q -l 50 mrna-seqs/e1.fasta | perl -lane 'if(/^(> \S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  > mrna-seqs/$1.pep
#/home/paulgardner/inst/squid-1.9g/translate -q -l 50 mrna-seqs/n1.fasta | perl -lane 'if(/^(> \S+)_s0_1to\d+/){print $1; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}' >> mrna-seqs/$1.pep

#~/inst/infernal-1.1.2/easel/miniapps/
esl-translate --watson -l  13 mrna-seqs/e1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1\.\..*frame=1\s+/){print ">$1"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  >  mrna-seqs/$1.pep
#~/inst/infernal-1.1.2/easel/miniapps/
esl-translate --watson -l  13 mrna-seqs/n1.fasta | perl -lane 'if(/^>orf\d+ source=(\S+) coords=1\.\..*frame=1\s+/){print ">$1"; $printme=1}elsif(/^>/){$printme=0;}elsif($printme){print;}'  >> mrna-seqs/$1.pep

echo "\n\n## Peptide sequences: [cat mrna-seqs/$1.pep]"
cat mrna-seqs/$1.pep

rm mrna-seqs/e1.fasta mrna-seqs/n1.fasta

hmmfetch ~/data/pfam/pfam31/Pfam-A.hmm  $1                       > $1.hmm
hmmalign                                $1.hmm mrna-seqs/$1.pep  > mrna-seqs/$1.pep.stk
cat mrna-seqs/$1.pep.stk

#~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat --mingap clustal mrna-seqs/$1.pep.stk | tr "." "-" > mrna-seqs/$1.pep.aln 
#~/inst/infernal-1.1.2/easel/miniapps/
esl-reformat --mingap clustal mrna-seqs/$1.pep.stk | tr "." "-" > mrna-seqs/$1.pep.aln 

echo "\n\n## Peptide alignment: [cat mrna-seqs/$1.pep.aln]"
cat mrna-seqs/$1.pep.aln

echo "## pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa"
pal2nal.pl mrna-seqs/$1.pep.aln mrna-seqs/$1.fasta  -output  fasta > mrna-seqs/$1.nuc.afa

#FAILS, FUCKING USELESS:
#echo "~/inst/RevTrans-1.4/revtrans.py -match name -vvv mrna-seqs/$1.fasta mrna-seqs/$1.pep.aln -O aln > mrna-seqs/$1.nuc.aln"
#~/inst/RevTrans-1.4/revtrans.py -match name -vvv mrna-seqs/$1.fasta mrna-seqs/$1.pep.aln -O aln > mrna-seqs/$1.nuc.aln

#~/inst/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk
#~/inst/infernal-1.1.2/easel/miniapps/
esl-reformat -u  --informat afa pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk
#sreformat -u --informat a2m pfam mrna-seqs/$1.nuc.afa > mrna-seqs/$1.nuc.stk

echo "\n\n## mRNA alignment (reverse translated): [cat mrna-seqs/$1.nuc.stk]"
cat mrna-seqs/$1.nuc.stk

