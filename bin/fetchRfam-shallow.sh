#!/bin/sh


grep "Name=$1;" U00096-ncRNAs-nonredundant.gff   | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4]\42   -f $F[3] -t $F[4] \42U00096;\42"'   | sh > ncrna-seqs/$1.fasta
grep "Name=$1;" AE014613-ncRNAs-nonredundant.gff | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d AE014613.fasta -r \42AE014613/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AE014613;\42"' | sh >> ncrna-seqs/$1.fasta


~/inst/infernal-1.1.1/src/cmfetch  ~/data/rfam/rfam12.2/Rfam.cm $1 > $1.cm
~/inst/infernal-1.1.1/src/cmalign $1.cm ncrna-seqs/$1.fasta   > ncrna-seqs/blah

~/inst/infernal-1.1.1/easel/miniapps/esl-reformat pfam ncrna-seqs/blah > ncrna-seqs/$1.stk

rm ncrna-seqs/blah
