#!/usr/bin/perl

#Generate new datafiles for R from "lofBroc-appended.Rdata" and "lofGFP-appended.Rdata" without all the headers.
#Additionally provide counts of LoF variants for each model


use warnings;
use strict;

open(bIN,     "< lofBroc-appended.Rdata" ) or die "FATAL: failed to open [lofBroc-appended.Rdata]\n[$!]";
open(gIN,     "<  lofGFP-appended.Rdata" ) or die "FATAL: failed to open  [lofGFP-appended.Rdata]\n[$!]";

open(bOUT,    "> lofBroc-appended2.Rdata" ) or die "FATAL: failed to open [lofBroc-appended2.Rdata]\n[$!]";
open(gOUT,    ">  lofGFP-appended2.Rdata" ) or die "FATAL: failed to open  [lofGFP-appended2.Rdata]\n[$!]";

open(bOUTcnt, "> lofBroc-counts.tsv" ) or die "FATAL: failed to open [lofBroc-counts.tsv]\n[$!]";
open(gOUTcnt, ">  lofGFP-counts.tsv" ) or die "FATAL: failed to open  [lofGFP-counts.tsv]\n[$!]";

processLoF(*bIN, *bOUT, *bOUTcnt);
processLoF(*gIN, *gOUT, *gOUTcnt);

exit 0;

sub processLoF {
    
    my ($in, $out, $outCnt) = @_;
    my ($cnt, $modelCnt)=(0,0);
    while(my $in = <$in>){
	chomp($in);
	if($in =~ /pos/){
	    my @in = split(/,/, $in);
	    my $printLine = join(",", @in[0..10]);
	    print $out "$printLine,\42LoF\42,\42model\42\n" if  $modelCnt == 0;
	    print $outCnt "$cnt\n" if $modelCnt != 0;
	    $modelCnt++;
	    $cnt=0;
	}
	else{
	    print $out "$in,\42$modelCnt\42\n";
	    $cnt++; 
	}
	
	
    }
    print $outCnt "$cnt\n";
    close($in); 
    close($out); 
    close($outCnt); 
    return 0;

}

