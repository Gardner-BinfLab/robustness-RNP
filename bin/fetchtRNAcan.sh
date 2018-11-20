#!/bin/sh

grep "type=$1;anti" U00096-ncRNAs.gff   | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d U00096.fasta   -r \42U00096/$F[3]-$F[4] \42   -f $F[3] -t $F[4] \42U00096;\42"'   | sh > ncrna-seqs/$1.fasta
grep "type=$1;anti" AL157959-ncRNAs.gff | sort -k6nr | head -n 1 | perl -lane '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq "-");  print "sfetch -d AL157959.fasta -r \42AL157959/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42AL157959;\42"' | sh >> ncrna-seqs/$1.fasta

~/inst/infernal-1.1.1/src/cmfetch  ~/data/rfam/rfam12.2/Rfam.cm tRNA > tRNA.cm
~/inst/infernal-1.1.1/src/cmalign tRNA.cm ncrna-seqs/$1.fasta   >> ncrna-seqs/$1.stk

#Realign Ser tRNA!!!
