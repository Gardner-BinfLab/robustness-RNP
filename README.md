# robustness-RNP
A collection of scripts, datasets, results and documents from an investigation of the genetic robustness of ncRNAs and proteins

## Figure 1: variation of recent (shallow) vs old (deep) divergence of RNA and protein families  

- this is sausage making, a lot of munging to pull out linked RNA & protein families from E. coli (ENA ID: U00096.3), S. enterica (ENA ID:AE014613.1) and N. meningitidis (ENA ID:AL157959.1):

* Dependencies include:
  * hmmer-3.1b2
  * pal2nal.pl
  * infernal-1.1.1

-Some error messages from pal2nal.pl due to phenylalanine/leucine bug in esl-translate.
-Alignments are unaffected though. 

```
cd data/genomes/

# Generate alignments for the different cohorts:

# Proteins, deep
../../bin/fetchPfam-deep.sh Ribosomal_L5 
../../bin/fetchPfam-deep.sh Ribosomal_L5_C
../../bin/fetchPfam-deep.sh tRNA-synt_1
../../bin/fetchPfam-deep.sh tRNA-synt_1_2 
../../bin/fetchPfam-deep.sh tRNA-synt_1d
../../bin/fetchPfam-deep.sh Arg_tRNA_synt_N
../../bin/fetchPfam-deep.sh B3_4
../../bin/fetchPfam-deep.sh Seryl_tRNA_N 
../../bin/fetchPfam-deep.sh tRNA-synt_2b 
../../bin/fetchPfam-deep.sh Ribonuclease_P
../../bin/fetchPfam-deep.sh Sigma70_ner
../../bin/fetchPfam-deep.sh GTP_EFTU
../../bin/fetchPfam-deep.sh S1
../../bin/fetchPfam-deep.sh SRP54_N
../../bin/fetchPfam-deep.sh SRP54
../../bin/fetchPfam-deep.sh SRP_SPB
../../bin/fetchPfam-deep.sh Ribosomal_S21
../../bin/fetchPfam-deep.sh Ribosomal_S15 
../../bin/fetchPfam-deep.sh Ribosomal_S17
../../bin/fetchPfam-deep.sh Ribosomal_L6
../../bin/fetchPfam-deep.sh ribosomal_L24
../../bin/fetchPfam-deep.sh Ribosomal_L31
../../bin/fetchPfam-deep.sh Ribosomal_S15 
../../bin/fetchPfam-deep.sh ThiC-associated
../../bin/fetchPfam-deep.sh ThiC_Rad_SAM 

# Proteins, shallow
../../bin/fetchPfam-shallow.sh Leader_Trp    Trp_leader
../../bin/fetchPfam-shallow.sh Leader_Thr    Thr_leader
../../bin/fetchPfam-shallow.sh SgrT          SgrS
../../bin/fetchPfam-shallow.sh HSP20         ROSE_2
../../bin/fetchPfam-shallow.sh MGTL          Mg_sensor
../../bin/fetchPfam-shallow.sh Leu_leader    Leu_leader
../../bin/fetchPfam-shallow.sh His_leader    His_leader
../../bin/fetchPfam-shallow.sh DHBP_synthase FMN
../../bin/fetchPfam-shallow.sh CSD           cspA
../../bin/fetchPfam-shallow.sh Plug          Cobalamin
../../bin/fetchPfam-shallow.sh TonB_dep_Rec  Cobalamin
../../bin/fetchPfam-shallow.sh Hfq
../../bin/fetchPfam-shallow.sh RNase_E_G 
../../bin/fetchPfam-shallow.sh ATP_bind_2
../../bin/fetchPfam-shallow.sh CsrA

# RNAs, deep
../../bin/fetchRfam-deep.sh RNaseP_bact_a
../../bin/fetchRfam-deep.sh 6S
../../bin/fetchRfam-deep.sh tmRNA
../../bin/fetchRfam-deep.sh Bacteria_small_SRP
../../bin/fetchRfam-deep.sh SSU_rRNA_bacteria
../../bin/fetchRfam-deep.sh LSU_rRNA_bacteria
../../bin/fetchRfam-deep.sh S15
../../bin/fetchRfam-deep.sh TPP
#tRNAs, deep
../../bin/fetchtRNAcan.sh Leu
../../bin/fetchtRNAcan.sh Arg
../../bin/fetchtRNAcan.sh Phe
../../bin/fetchtRNAcan.sh Ser

# RNAs, shallow
../../bin/fetchRfam-shallow.sh Cobalamin
../../bin/fetchRfam-shallow.sh cspA
../../bin/fetchRfam-shallow.sh CsrB
../../bin/fetchRfam-shallow.sh FMN
../../bin/fetchRfam-shallow.sh GcvB
../../bin/fetchRfam-shallow.sh GlmY_tke1
../../bin/fetchRfam-shallow.sh GlmZ_SraJ
../../bin/fetchRfam-shallow.sh His_leader
../../bin/fetchRfam-shallow.sh Leu_leader
../../bin/fetchRfam-shallow.sh Mg_sensor
../../bin/fetchRfam-shallow.sh MicA
../../bin/fetchRfam-shallow.sh OxyS
../../bin/fetchRfam-shallow.sh ROSE_2
../../bin/fetchRfam-shallow.sh rseX
../../bin/fetchRfam-shallow.sh SgrS
../../bin/fetchRfam-shallow.sh Thr_leader
../../bin/fetchRfam-shallow.sh Trp_leader

# Run "computeSynonNonsynon.pl" over the alignments:
ls -1 ncrna-seqs/*.stk     | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 1 -i "$1}' | sh > synonNonsynon-ncRNA-results.txt
ls -1  mrna-seqs/*.nuc.stk | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 0 -i "$1}' | sh > synonNonsynon-mRNA-results.txt

# Evil hoop-jumping to get into tabular format:
echo -e "ID\tdistID\ttotVar\tlen\taaSynTot\taaSynProp\taaNonSynTot\taaNonSynProp\tbanpSynTot\tbanpSynProp\tbanpNonSynTot\tbanpNonSynProp\tblossSynTot\tblossSynProp\tblossNonSynTot\tblossNonSynProp" > synonNonsynon-mRNA-results.tsv
cat synonNonsynon-mRNA-results.txt | perl -lane 'if(  /^mrna-seqs\/(\S+)\.nuc/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/aaSynon:\s+(\d+) \((\S+)\)\s+aaNonsynon:\s+(\d+) \((\S+)\)/){($aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP)=($1,$2,$3,$4);}elsif(/banpSynon:\s+(\d+) \((\S+)\)\s+banpNonsynon:\s+(\d+) \((\S+)\)/){($banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP)=($1,$2,$3,$4);}elsif(/blossumSynon:\s+(\d+) \((\S+)\)\s+blossumNonsynon:\s+(\d+) \((\S+)\)/){($blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP)=($1,$2,$3,$4); %salm=(ATP_bind_2=>1,CSD=>1,CsrA=>1,DHBP_synthase=>1,Hfq=>1,His_leader=>1,HSP20=>1,Leader_Thr=>1,Leader_Trp=>1,Leu_leader=>1,MGTL=>1,Plug=>1,RNase_E_G=>1,SgrT=>1,TonB_dep_Rec=>1); $dist='ecol2nmen'; $dist='ecol2styp' if defined(%salm{$id});   @a=($id,$dist,$tV,$len,$aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP,$banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP,$blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP); $a=join("\t",@a); print $a}' >> synonNonsynon-mRNA-results.tsv

echo -e "ID\tdistID\ttotVar\tlen\trySynTot\trySynProp\tryNonSynTot\tryNonSynProp\tbpSynTot\tbpSynProp\tbpNonSynTot\tbpNonSynProp" > synonNonsynon-ncRNA-results.tsv
cat synonNonsynon-ncRNA-results.txt | perl -lane 'if(/^ncrna-seqs\/(\S+)\.stk/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/rySynon:\s+(\d+) \((\S+)\)\s+ryNonsynon:\s+(\d+) \((\S+)\)/){($rySynT,$rySynP,$ryNonSynT,$ryNonSynP)=($1,$2,$3,$4);}elsif(/bpSynon:\s+(\d+) \((\S+)\)\s+bpNonsynon:\s+(\d+) \((\S+)\)/){($bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP)=($1,$2,$3,$4); %salm=(Cobalamin=>1,cspA=>1,CsrB=>1,FMN=>1,GcvB=>1,GlmY_tke1=>1,GlmZ_SraJ=>1,His_leader=>1,Leu_leader=>1,Mg_sensor=>1,MicA=>1,OxyS=>1,ROSE_2=>1,rseX=>1,SgrS=>1,Thr_leader=>1,Trp_leader=>1); $dist='ecol2nmen'; $dist='ecol2styp' if defined(%salm{$id}); @a=($id,$dist,$tV,$len,$rySynT,$rySynP,$ryNonSynT,$ryNonSynP,$bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP); $a=join("\t",@a); print $a }' >> synonNonsynon-ncRNA-results.tsv

# Clean up files:
# rm *cm *hmm blah ncrna-seqs/*fasta mrna-seqs/*pep mrna-seqs/*fasta mrna-seqs/*afa mrna-seqs/*aln
```

* FINAL RESULTS:

..* synonNonsynon-ncRNA-results.tsv

..* synonNonsynon-mRNA-results.tsv

* Plot graphs:
```
../../bin/plotConservation.R
######################################################################
```

## Figure 2: Robustness of structure predictions to random in silico mutagenesis for a protein (SgrT) and non-coding RNA (SgrS).

```
## SgrS/SgrT structure 
### RNA & protein point mutations

#WARNING: these structure predictions use a lot of CPU...
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
../../bin/plotStructureMutagenerator.R
######################################################################

```

## Figure 3: Relative fluorescence intensities of mutated RNA Broccoli and mutated protein mCherry.

* Plot graphs:

```
../bin/plotFluoro.R
```

