
## Experiment 2

## \Delta \Delta G and \Delta bitscore analysis of PDB RNPs  

18 July 2019

* Fetch all PDB sequences from here:
 * ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz

* Selecting candidate RNPs:
```
 * cat pdb_seqres.txt | perl -lane 'if(/^>.*(protein|DNA|ribosomal|tRNA|transfer|mRNA|rRNA|T\-RNA)/i){$p=0;}elsif(/^>.*length:(\d+).*RNA/){if (50<$1 && $1<150){$p=1;}else{$p=0}} print if $p'
```

* selection criteria 
 * structured RNA
 * structured protein
 * "unique" i.e. not too many tRNAs and synthetases
 * ncRNA length >50 and <150, protein length <500 and >50 (although some exceptions have been made)
 * RNA sequence maps to Rfam, and protein sequence maps to EggNOG
 
The selected protein sequences are saved to "protein-pdb-seqs.fa" and corresponding ncRNA sequences are saved to "ncRNA-pdb-seqs.fa". 

The (best) matching messenger RNA sequences corresponding to the
protein-pdb-seqs were determined using tblastn searches of the NR
database. If necessary, mRNAs were manually corrected to match the
corresponding PDB sequences. 

### QC selected proteins, check mRNA sequences, ...
```
esl-seqstat -a protein-pdb-seqs.fa > protein-pdb-seqs.seqstat
esl-seqstat -a mRNA-pdb-seqs.fa    > mRNA-pdb-seqs.seqstat
paste protein-pdb-seqs.seqstat  mRNA-pdb-seqs.seqstat | perl -lane 'if(/^=\s+(\S+)\s+(\d+).*=\s+(\S+)\s+(\d+)/){printf "%d\t$1\t$2\t$4\t$3\n", 3*$2-$4;}'
```
### fetch pdb files: 
```
cat protein-pdb-seqs.fa | perl -lane 'if(/^>(\S+)_/){$p=$1;  $p =~ tr/a-z/A-Z/;  print "curl -G https://files.rcsb.org/download/$p.pdb.gz > pdbfiles/$1.pdb.gz"}' | sh
```

### Find mappings between the RNA sequences and Rfam (covariance models):
```
cd hmm-cm-models
```
* Rfam matches:
 * "--cut_ga" removed - not using GA thresholds as some matches are partial: 
```
~/inst/infernal-1.1.2/src/cmscan -o ncRNA-pdb-seqs.rfam14.cmscan  --tblout  ncRNA-pdb-seqs.rfam14.cmscan.tblout --fmt 2 -T 10 --clanin ~/data/rfam/rfam14.1/Rfam.clanin --oskip ~/data/rfam/rfam14.1/Rfam.cm  ../ncRNA-pdb-seqs.fa
```

* Edited list of unique matches saved to "model2pdb-maps.txt"

* fetch cm models: 
```
cat model2pdb-maps.txt | cut -f 1 | sort -d | uniq | awk '{print "~/inst/infernal-1.1.2/src/cmfetch ~/data/rfam/rfam14.1/Rfam.cm "$1" > "$1".cm"}' | sh
```

* 5u30_B has no corresponding Rfam model, building a local one: 
```
~/inst/infernal-1.1.2/src/cmbuild sgRNA.cm 5u30_B.stk
~/inst/infernal-1.1.2/src/cmcalibrate --cpu 8  sgRNA.cm
echo -e "sgRNA\t5u30_B" >> model2pdb-maps.txt
```

### Find mappings between the protein sequences and EggNOG (profile HMMs): 
```
curl -G http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_hmms.tar > bacteria.tar
curl -G http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/2759_hmms.tar > eukaryota.tar
```

* batch search: http://eggnog-mapper.embl.de/job_status?jobname=MM_3_8x_otq
 * Version 1 & 2 of the eggNOG mapper: eggnog2pdb-proteins.tsv  eggnog-mapper2pdb-proteins.tsv

* Manually downloaded HMMs matching each protein, and corrected the names in the HMM files. E.g. http://eggnogapi5.embl.de/nog_data/file/hmm/COG0050

* concatenated HMMs are available in this file:
```
hmmpress HMM
```
* truncated sequences may score differently compared to full length seqs: 
```
grep ^">" ../mRNA-pdb-seqs.fa  |  perl -lane 'if(/^>(\S+)\/(\d+)-(\d+)/){print "esl-sfetch -c $2..$3 -n $1  ../protein-pdb-seqs.fa $1"}elsif(/^>(\S+)/){print "esl-sfetch ../protein-pdb-seqs.fa $1"}' | sh  >  protein-pdb-seqs-trunced.fa
hmmscan --tblout protein-pdb-seqs.hmmscan.tblout  HMM protein-pdb-seqs-trunced.fa > protein-pdb-seqs.hmmscan.out
```

### Removed these problematic structures:
* 5u30 (CRISPR-associated endonuclease C2c1) sequence is too long
* 4lck problems running maestro
```
../../bin/computeDeltaDelta.pl > computeDeltaDelta-1mut.out
```
* remove failed runs: 
```
egrep '^ID|protein|rna' computeDeltaDelta-1mut.out | egrep -v 'Failed|Warning' | perl -lane 'print if($F[2] !~ /NA/)' > computeDeltaDelta-1mut.dat
```

* update "\$numMutations" & rerun:
```
../../bin/computeDeltaDelta.pl > computeDeltaDelta-4mut.out
egrep '^ID|protein|rna' computeDeltaDelta-4mut.out | egrep -v 'Failed|Warning' | perl -lane 'print if(scalar(@F) == 5)' > computeDeltaDelta-4mut.dat
```

### Small scale test (big runs take ~24 hours on a high-spec laptop)
```
grep -A 1   6nd4_2 ncRNA-pdb-seqs.fa   > 6nd4_2.fa
grep -A 11  6nd4_c mRNA-pdb-seqs.fa    > 6nd4_c-mRNA.fa   
grep -A 1   6nd4_c protein-pdb-seqs.fa > 6nd4_2-protein.fa     
../../bin/computeDeltaDelta.pl -v -numMutants 10 -numMutations 1 -n 6nd4_2.fa  -m 6nd4_c-mRNA.fa   -p  6nd4_2-protein.fa     

#2fmt_A   2fmt_C
grep -A 1   2fmt_C ncRNA-pdb-seqs.fa   > 2fmt_C.fa
grep -A 19  2fmt_A mRNA-pdb-seqs.fa    > 2fmt_A-mRNA.fa   
grep -A 1   2fmt_A protein-pdb-seqs.fa > 2fmt_A-protein.fa     
../../bin/computeDeltaDelta.pl -v -numMutants 10 -numMutations 1 -n 2fmt_C.fa  -m 2fmt_A-mRNA.fa   -p  2fmt_A-protein.fa     
```



### The LARGE scale test
```
../../bin/computeDeltaDelta.pl -v -numMutations 1 > computeDeltaDelta-1mut.out &&
../../bin/computeDeltaDelta.pl -v -numMutations 4 > computeDeltaDelta-4mut.out
```

* make a summary table of the results:
```
cat protein-pdb-seqs.fa | perl -lane 'if(/^>(\S+)\s+\S+\s+length:(\d+)\s+(.*)/ or /^>(\S+)\s+offset=\S+\s+\S+\s+length:(\d+)\s+(.*)/  ){print "$1\t$2\t$3"}elsif(/^>/){print "WARNING: $_"}' > blah1 
cat   ncRNA-pdb-seqs.fa | perl -lane 'if(/^>(\S+)\s+\S+\s+length:(\d+)\s+(.*)/ or /^>(\S+)\s+offset=\S+\s+\S+\s+length:(\d+)\s+(.*)/  ){print "$1\t$2\t$3"}elsif(/^>/){print "WARNING: $_"}' > blah2
echo -e "Protein.PDBID\tProtein.Length\tProtein.Description\tRNA.PDBID\tRNA.Length\tRNA.Description" > RNP-summary-table.tsv
paste blah1 blah2 >> RNP-summary-table.tsv
```
