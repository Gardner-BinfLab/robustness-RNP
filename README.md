# robustness-RNP
A collection of scripts, datasets, results and documents from an investigation of the robustness of ncRNAs and proteins to mutation 

---

######################################################################

## Experiment 1

## Figure 1: mutagenesis of a fluorescent RNA and protein 

- This is to map 

* Dependencies include:
  * bowtie2 (version 2.4.4)
  
```
cd data/fluoro-II

# create bowtie2 index databases
bowtie2-build broccoli.fasta  broccoli
bowtie2-build     egfp.fasta  egfp

# bowtie2 mapping of sequencing reads to broccoli/eGFP fasta files 
bowtie2 -x broccoli -1 fastq/broccoli_R1.fastq -2 fastq/broccoli_R2.fastq  -S broccoli.sam
bowtie2 -x     egfp -1     fastq/egfp_R1.fastq -2     fastq/egfp_R2.fastq  -S     egfp.sam
#NOTE: fastq and sam files are managed with git-lfs (large file management)

# Call variants (low PHRED threshold in the first round "-p 25", and give the region where variants are expected "-f 118 -t 318")
# Conservative mode (-c) used to identify variants found in both forward and reverse reads (labelled "T" for true, otherwise "F")
../../bin/sam2var.pl -b broccoli-barcodes.txt -r ./broccoli.fasta  -s ./broccoli.sam -p 25 -f  118 -t  318 -c  | sort -n   > sam2var.broccoli-118-318.txt
../../bin/sam2var.pl -b     egfp-barcodes.txt -r     ./egfp.fasta  -s     ./egfp.sam -p 25 -f 2410 -t 2610 -c  | sort -n   > sam2var.egfp-2410-2610.txt

# Calculate accuracy statistics to identify 
../../bin/trainPhredThresh.pl -i  broccoli.sam.variation  --from 1 --to 200
../../bin/trainPhredThresh.pl -i      egfp.sam.variation  --from 1 --to  75

#BROCCOLI
# phred	sens	spec	fdr	fpr	f1	mcc	tp	fp	fn	tn
# 40	0.00	1.00	0.25	0.00	0.00	0.00	6	2	229730	108726
# 39	0.26	0.77	0.30	0.23	0.38	0.03	60809	25496	168927	83232
# 38	0.48	0.54	0.32	0.46	0.56	0.01	109162	50490	120574	58238
# 37	0.68	0.44	0.28	0.56	0.70	0.12	156254	60792	73482	47936
# 36	0.73	0.39	0.28	0.61	0.73	0.13	168469	65817	61267	42911
# 35	0.76	0.37	0.28	0.63	0.74	0.14	175129	68378	54607	40350
# 34	0.80	0.33	0.28	0.67	0.76	0.15	184296	72536	45440	36192
# 33	0.85	0.30	0.28	0.70	0.78	0.18	196406	76034	33330	32694
# 32	0.89	0.28	0.28	0.72	0.80	0.21	203569	78249	26167	30479
# 31	0.90	0.26	0.28	0.74	0.80	0.21	206561	80479	23175	28249
# 30	0.93	0.21	0.29	0.79	0.81	0.21	213236	85558	16500	23170 <--  
# 29	0.94	0.19	0.29	0.81	0.81	0.20	215030	87763	14706	20965 <--  
# 28	0.94	0.17	0.29	0.83	0.81	0.18	216504	90184	13232	18544 <--  
# 27	0.97	0.06	0.32	0.94	0.80	0.06	222332	102279	7404	6449
# 26	0.98	0.04	0.32	0.96	0.81	0.06	225998	104709	3738	4019
# 25	1.00	0.00	0.32	1.00	0.81	0.00	229736	108728	0	0

#eGFP
# phred	sens	spec	fdr	fpr	f1	mcc	tp	fp	fn	tn
# 39	0.39	0.55	0.04	0.45	0.56	-0.02	36602	1426	57220	1770
# 38	0.65	0.41	0.03	0.59	0.78	0.02	60902	1870	32920	1326
# 37	0.96	0.21	0.03	0.79	0.97	0.15	90308	2535	3514	661
# 36	0.97	0.19	0.03	0.81	0.97	0.16	91090	2579	2732	617
# 35	0.98	0.17	0.03	0.83	0.97	0.16	91614	2646	2208	550
# 34	0.98	0.16	0.03	0.84	0.98	0.16	92018	2693	1804	503
# 33	0.99	0.10	0.03	0.90	0.98	0.16	93022	2871	800	325
# 32	0.99	0.09	0.03	0.91	0.98	0.16	93235	2907	587	289
# 31	0.99	0.09	0.03	0.91	0.98	0.16	93332	2918	490	278
# 30	1.00	0.08	0.03	0.92	0.98	0.17	93470	2949	352	247 <--  
# 29	1.00	0.07	0.03	0.93	0.98	0.16	93488	2962	334	234 <--  
# 28	1.00	0.07	0.03	0.93	0.98	0.16	93520	2979	302	217 <--  
# 27	1.00	0.04	0.03	0.96	0.98	0.12	93617	3059	205	137
# 26	1.00	0.01	0.03	0.99	0.98	0.05	93724	3158	98	38
# 25	1.00	0.00	0.03	1.00	0.98	0.00	93822	3196	0	0

# Both report PHRED score of 29 as a good balance between true and false positives based upon MCC: 
../../bin/sam2var.pl -b broccoli-barcodes.txt -r ./broccoli.fasta  -s ./broccoli.sam -p 29 -f  118 -t  318   | sort -n   > sam2var.broccoli-118-318.txt
../../bin/sam2var.pl -b     egfp-barcodes.txt -r     ./egfp.fasta  -s     ./egfp.sam -p 29 -f 2410 -t 2610   | sort -n   > sam2var.egfp-2410-2610.txt

# Plot summaries of variant calls:
R CMD BATCH ../../bin/plotMutagenesis.R ; tail -n 20 plotMutagenesis.Rout

rm -f *appended.Rdata
#Simulate data, train a RF, predict variant classes based on bin freqs. Repeat 100x
for i in {1..100}; do echo "ROUND $i:"; ../../bin/runTrainClassifier2.sh; done

../../bin/summariseRFFire.sh

```


######################################################################

## Experiment 2

## Figure 2: variation of recent (shallow) vs old (deep) divergence of RNA and protein families  

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

## Experiment 3

## Figure 3A&B: 

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

## Figure 3C&D: Robustness of structure predictions to random in silico mutagenesis for a protein (SgrT) and non-coding RNA (SgrS).

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

## Supplemental Experiment

## Supplemental Figure 3: Relative fluorescence intensities of mutated RNA Broccoli and mutated protein mCherry.

* Plot graphs:

```
../bin/plotFluoro.R
```

