# robustness-RNP
A collection of scripts, datasets, results and documents from an investigation of the robustness of ncRNAs and proteins to mutation 

---

## Experiment 1

## Figure 1: variation of recent (shallow) vs old (deep) divergence of RNA and protein families  

- this is to pull out linked RNA & protein families from E. coli (ENA ID: U00096.3), S. enterica (ENA ID:AE014613.1) and N. meningitidis (ENA ID:AL157959.1):

* Dependencies include:
  * hmmer-3.1b2
  * pal2nal.pl
  * infernal-1.1.2
  * tRNAscan-SE-1.3.1 

-Some error messages from pal2nal.pl due to phenylalanine/leucine bug in esl-translate.
-Alignments are unaffected though. 

```
cd data/genomes/

# Generate alignments for the deep and shallow ncRNAs and protein domains (bacteria first, fungi second):

# Proteins, deep
grep deep    rna-protein-pairs.tsv | cut -f 2 | sort -d | uniq | awk '{print "../../bin/fetchPfam.sh "$1" U00096-pfam31.gff  U00096.fasta  AL157959-pfam31.gff    AL157959.fasta"}' | sh
# Proteins, shallow
grep shallow rna-protein-pairs.tsv | cut -f 2 | sort -d | uniq | awk '{print "../../bin/fetchPfam.sh "$1" U00096-pfam31.gff  U00096.fasta  AE014613-pfam31.gff    AE014613.fasta"}' | sh

# RNAs, deep
grep deep    rna-protein-pairs.tsv | cut -f 1 | sort -d | uniq | awk '{print "../../bin/fetchRfam.sh "$1" U00096-ncRNAs.gff  U00096.fasta  AL157959-ncRNAs.gff    AL157959.fasta"}' | sh
# RNAs, shallow
grep shallow rna-protein-pairs.tsv | cut -f 1 | sort -d | uniq | awk '{print "../../bin/fetchRfam.sh "$1" U00096-ncRNAs.gff  U00096.fasta  AE014613-ncRNAs.gff    AE014613.fasta"}' | sh

#cd -
cd ../genomes-fungi

# Proteins, deep
grep deep    rna-protein-pairs.tsv | cut -f 2 | sort -d | uniq | awk '{print "../../bin/fetchPfam.sh "$1" cerevisiae-pfam32.gff  cerevisiae.fa        pombe-pfam32.gff         pombe.fa "}' | sh
# Proteins, shallow
grep shallow rna-protein-pairs.tsv | cut -f 2 | sort -d | uniq | awk '{print "../../bin/fetchPfam.sh "$1" cerevisiae-pfam32.gff  cerevisiae.fa kudriavzevii-pfam32.gff  kudriavzevii.fa "}' | sh

# RNAs, deep
grep deep    rna-protein-pairs.tsv | cut -f 1 | sort -d | uniq | awk '{print "../../bin/fetchRfam.sh "$1" cerevisiae-ncRNAs.gff  cerevisiae.fa        pombe-ncRNAs.gff         pombe.fa "}' | sh
# RNAs, shallow
grep shallow rna-protein-pairs.tsv | cut -f 1 | sort -d | uniq | awk '{print "../../bin/fetchRfam.sh "$1" cerevisiae-ncRNAs.gff  cerevisiae.fa kudriavzevii-ncRNAs.gff  kudriavzevii.fa "}' | sh

# The Rfam Group I Intron misses this sequence, so here's a work around:
~/inst/infernal-1.1.2/src/cmalign  introns/IC2.cm   ncrna-seqs/Intron_gpI-IC2.fasta > ncrna-seqs/blah && ~/inst/infernal-1.1.2/easel/miniapps/esl-reformat pfam ncrna-seqs/blah > ncrna-seqs/Intron_gpI-IC2.stk


######################################################################

# Run "computeSynonNonsynon.pl" over the alignments:
ls -1 ncrna-seqs/*.stk     | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 1 -i "$1}' | sh > synonNonsynon-ncRNA-results.txt
ls -1  mrna-seqs/*.nuc.stk | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 0 -i "$1}' | sh > synonNonsynon-mRNA-results.txt

# Some hoop-jumping to get into tabular format:
echo -e "ID\tdistID\ttotVar\tlen\taaSynTot\taaSynProp\taaNonSynTot\taaNonSynProp\tbanpSynTot\tbanpSynProp\tbanpNonSynTot\tbanpNonSynProp\tblossSynTot\tblossSynProp\tblossNonSynTot\tblossNonSynProp" > synonNonsynon-mRNA-results.tsv
cat synonNonsynon-mRNA-results.txt | perl -lane 'if(  /^mrna-seqs\/(\S+)\.nuc/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/aaSynon:\s+(\d+) \((\S+)\)\s+aaNonsynon:\s+(\d+) \((\S+)\)/){($aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP)=($1,$2,$3,$4);}elsif(/banpSynon:\s+(\d+) \((\S+)\)\s+banpNonsynon:\s+(\d+) \((\S+)\)/){($banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP)=($1,$2,$3,$4);}elsif(/blossumSynon:\s+(\d+) \((\S+)\)\s+blossumNonsynon:\s+(\d+) \((\S+)\)/){($blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP)=($1,$2,$3,$4); $gg = `grep $id rna-protein-pairs.tsv `; @gg = split(/\t/, $gg); $dist=$gg[2]; @a=($id,$dist,$tV,$len,$aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP,$banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP,$blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP); $a=join("\t",@a); print $a}' >> synonNonsynon-mRNA-results.tsv

echo -e "ID\tdistID\ttotVar\tlen\trySynTot\trySynProp\tryNonSynTot\tryNonSynProp\tbpSynTot\tbpSynProp\tbpNonSynTot\tbpNonSynProp" > synonNonsynon-ncRNA-results.tsv
cat synonNonsynon-ncRNA-results.txt | perl -lane 'if(/^ncrna-seqs\/(\S+)\.stk/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/rySynon:\s+(\d+) \((\S+)\)\s+ryNonsynon:\s+(\d+) \((\S+)\)/){($rySynT,$rySynP,$ryNonSynT,$ryNonSynP)=($1,$2,$3,$4);}elsif(/bpSynon:\s+(\d+) \((\S+)\)\s+bpNonsynon:\s+(\d+) \((\S+)\)/){($bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP)=($1,$2,$3,$4); $gg = `grep $id rna-protein-pairs.tsv `; @gg = split(/\t/, $gg); $dist=$gg[2]; @a=($id,$dist,$tV,$len,$rySynT,$rySynP,$ryNonSynT,$ryNonSynP,$bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP); $a=join("\t",@a); print $a }' >> synonNonsynon-ncRNA-results.tsv

# Repeat for bacteria and fungi

# Clean up files:
# rm *cm *hmm blah ncrna-seqs/*fasta mrna-seqs/*pep mrna-seqs/*fasta mrna-seqs/*afa mrna-seqs/*aln
```

* FINAL RESULTS:
  * synonNonsynon-ncRNA-results.tsv
  * synonNonsynon-mRNA-results.tsv

* Plot graphs:
```
cd ../../
./bin/plotConservation.R
convert  -flatten -density 300 -trim  manuscript/figures/figure1.pdf     -quality 100   manuscript/figures/figure1.png
convert  -flatten -density 300 -trim  manuscript/figures/figure1c.pdf    -quality 100   manuscript/figures/figure1c.png
convert  -flatten -density 200 -trim  manuscript/figures/suppfigure1.pdf -quality 100   manuscript/figures/suppfigure1.png
######################################################################
```

---

## Experiment 2

## Figure 2A&B: 

See data/delta-delta/README.md

* Dependencies include:
  * RNAfold (v2.4.9)
  * MAESTRO (v1.1)
  * HMMER (v3.1b2)
  * INFERNAL (v1.1.2)

```
cd data/delta-delta/

../../bin/computeDeltaDelta.pl -v -numMutations 1 > computeDeltaDelta-1mut.out &&
../../bin/computeDeltaDelta.pl -v -numMutations 4 > computeDeltaDelta-4mut.out
```

* Plot graphs:

```
cd ../..
./bin/plotDelta.R
```

## Figure 2C&D: Robustness of structure predictions to random in silico mutagenesis for a protein (SgrT) and non-coding RNA (SgrS).

* Dependencies include:
  * RNAfold (v2.4.9)
  * I-TASSER (v5.1)

## SgrS/SgrT structure 
### RNA & protein point mutations

```
#WARNING: this simulation uses a lot of CPU...
for name in 0.001 0.01 0.05 0.1 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50
do
echo "./bin/structureMutagenerator.pl -dir sgr-structural -i data/sgr-structural/sgrS-sRNA.fasta  -m 100 -mrate $name  -indelrate 0.0    > data/sgr-structural/tmp/sgrS-sRNA-mrate$name.txt"
echo "./bin/structureMutagenerator.pl -dir sgr-structural -i data/sgr-structural/sgrT-mRNA.fasta  -m 100 -mrate $name  -indelrate 0.0 -p > data/sgr-structural/tmp/sgrT-mRNA-mrate$name.txt"
done | sh

### RNA & protein INDEL mutations
for name in 0.0001 0.001 0.01 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50
do
echo "./bin/structureMutagenerator.pl -dir sgr-structural -i data/sgr-structural/sgrS-sRNA.fasta  -m 100 -indelrate $name  -mrate 0.0    > data/sgr-structural/tmp/sgrS-sRNA-indelrate$name.txt"
echo "./bin/structureMutagenerator.pl -dir sgr-structural -i data/sgr-structural/sgrT-mRNA.fasta  -m 100 -indelrate $name  -mrate 0.0 -p > data/sgr-structural/tmp/sgrT-mRNA-indelrate$name.txt"
done | sh

cd data/sgr-structural/tmp
ls -1 | awk '{print "egrep -v \42^ID|^ECOL\42 "$1" > ../"$1}' | sh

```
* Plot graphs:
```
cd data/sgr-structural/
../../bin/plotStructureMutagenerator.R
```

---

## Experiment 3

## Figure 3: Relative fluorescence intensities of mutated RNA Broccoli and mutated protein mCherry.

* Plot graphs:

```
../bin/plotFluoro.R
```

