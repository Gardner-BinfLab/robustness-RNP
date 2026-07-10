#!/usr/bin/bash

# Execute in a loop with:
# rm *appended.Rdata
# for i in {1..100}; do echo "ROUND $i:"; ./runTrainClassifier2.sh; done

rm -f .RData trainClassifier2.Rout;
R CMD BATCH ../../bin/trainClassifier2.R &&
    cat gofBroc.Rdata    >>   gofBroc-appended.Rdata &&
    cat  gofGFP.Rdata    >>    gofGFP-appended.Rdata &&
    cat lofBroc.Rdata    >>   lofBroc-appended.Rdata &&
    cat  lofGFP.Rdata    >>    lofGFP-appended.Rdata &&
    cat  brocProbs.Rdata >> brocProbs-appended.Rdata &&
    cat   gfpProbs.Rdata >>  gfpProbs-appended.Rdata ;
#
tail trainClassifier2.Rout > blah;
#
echo "Brocolli GoF" &&
    grep -v pos gofBroc-appended.Rdata | cut -f 2 -d "," | sort | uniq -c | sort -nr | head ;
echo "Brocolli LoF" &&
    grep -v pos lofBroc-appended.Rdata | cut -f 2 -d "," | sort | uniq -c | sort -nr | head ;
echo "GFP GoF" &&
    grep -v pos gofGFP-appended.Rdata | cut -f 2 -d "," | sort | uniq -c | sort -nr | head;
echo "GFP LoF" &&
    grep -v pos lofGFP-appended.Rdata | cut -f 2 -d "," | sort | uniq -c | sort -nr | head;


