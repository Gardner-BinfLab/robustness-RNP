#!/usr/bin/bash

head -n 1 gofBroc-appended.Rdata | cut -d ',' -f 1-12 > gofBroc-appended2.Rdata  &&
    grep -v "pos" gofBroc-appended.Rdata | cut -d ',' -f 1-12  | tr "T" "U" >> gofBroc-appended2.Rdata;
head -n 1 gofGFP-appended.Rdata | cut -f 1-12 > gofGFP-appended2.Rdata &&
    grep -v "pos" gofGFP-appended.Rdata  | cut -d ',' -f 1-12 >> gofGFP-appended2.Rdata

#READ DATA FROM:
#lofBroc-appended.Rdata
# lofGFP-appended.Rdata

#WRITE ANNOTATED DATA:
#lofBroc-appended2.Rdata
# lofGFP-appended2.Rdata

#Counts of LoF variants for each model:
#lofBroc-counts.tsv
# lofGFP-counts.tsv

../../bin/summariseRFFire.pl
R CMD BATCH ../../bin/summariseRFFire.R



