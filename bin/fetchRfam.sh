#!/bin/bash

ID=$1
GFFFILE1=$2
SEQFILE1=$3
GFFFILE2=$4
SEQFILE2=$5

if [[ $ID == *tRNA* ]];
then
    #Using arrays in bash!: 
    IFS='-'
    read -ra ADDR <<< "$ID"   
    ID="${ADDR[0]}"
    FID="${ADDR[0]}-tRNA"
    CMID="${ADDR[1]}"
    GREPSTR="type=$ID;anti"
    
else
    ID=$ID
    FID=$ID
    CMID=$ID
    GREPSTR="Name=$ID;|rfam-id=$ID;"
fi
IFS=' '
echo -ne "\n\nID[$ID] GFFFILE1[$GFFFILE1] SEQFILE1[$SEQFILE1] GFFFILE2[$GFFFILE2] SEQFILE2[$SEQFILE2] CMID[$CMID] GREPSTR[$GREPSTR]\n"

mkdir -p ncrna-seqs
egrep "$GREPSTR" $GFFFILE1 | sort -k6nr | head -n 1 | perl -lane  '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq  "-" );  print  "sfetch -d '"$SEQFILE1"' -r \42$F[0]/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42$F[0]\42"' | sh  > ncrna-seqs/$FID.fasta
egrep "$GREPSTR" $GFFFILE2 | sort -k6nr | head -n 1 | perl -lane  '($F[3],$F[4])=($F[4],$F[3]) if ($F[6] eq  "-" );  print  "sfetch -d '"$SEQFILE2"' -r \42$F[0]/$F[3]-$F[4]\42 -f $F[3] -t $F[4] \42$F[0]\42"' | sh >> ncrna-seqs/$FID.fasta

~/inst/infernal-1.1.1/src/cmfetch  ~/data/rfam/rfam12.2/Rfam.cm $CMID  > $CMID.cm
~/inst/infernal-1.1.1/src/cmalign $CMID.cm ncrna-seqs/$FID.fasta        > ncrna-seqs/blah
#esl-alimanip --detrunc 5 --trim
~/inst/infernal-1.1.1/easel/miniapps/esl-reformat pfam ncrna-seqs/blah > ncrna-seqs/$FID.stk

echo -ne "ncrna-seqs/$ID.stk\n"
cat ncrna-seqs/$FID.stk
rm ncrna-seqs/blah


