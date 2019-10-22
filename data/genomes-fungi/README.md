## Experiment 1

## Figure 1B: variation of recent (shallow) vs old (deep) divergence of RNA and protein families  


* Compare 3 fungal genomes, find shared RNA and protein families with linked functions
 * *Saccharomyces cerevisiae*,  *Saccharomyces kudriavzevii* and *Schizosaccharymocyes pombe*


* ***Saccharomyces cerevisiae***

* https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/

```grep ^chr README  | awk '{print "wget https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/"$1}'  | sh
cat chr*fsa > s-cerevisiae-S288C_reference.fsa && rm chr*fsa
```

* ***Schizosaccharymocyes pombe*** 
 * ftp://ftp.pombase.org/pombe/genome_sequence_and_features/genome_sequence/


* ***Saccharomyces kudriavzevii***
 * https://downloads.yeastgenome.org/sequence/fungi/S_kudriavzevii/IFO1802/
 * IFO1802_UColDMed_2011_SRX055455.fsa
```
mv IFO1802_UColDMed_2011_SRX055455.fsa kudriavzevii.fa
```

* Annotate the ncRNAs:
```
cd ~/projects/robustness-RNP/data/genomes-fungi
basename -a -s '.fa' `ls *fa ` | awk '{print "~/inst/infernal-1.1/src/cmsearch -o "$1"-ncRNAs-with-rfam130.cmsearch --tblout  "$1"-ncRNAs-with-rfam130.tbl  --cut_ga --rfam  ~/data/rfam/rfam14.0/Rfam.cm "$1".fa"}' | sh
basename -a -s '-ncRNAs-with-rfam130.tbl' `ls $tRNAscanOut/*-ncRNAs-with-rfam130.tbl` | awk '{print "~/bin/cmsearchtblout2gff.pl $tRNAscanOut/"$1"-ncRNAs-with-rfam130.tbl \t> $tRNAscanOut/"$1"-ncRNAs-with-rfam130.gff"}'  | sh

export tRNAscanOut=~/projects/robustness-RNP/data/genomes-fungi
cd ~/inst/tRNAscan-SE-1.3.1/
basename -a -s '.fa' `ls $tRNAscanOut/*fa` | awk '{print "tRNAscan-SE   --general --output  $tRNAscanOut/"$1"-tRNAscan131.out $tRNAscanOut/"$1".fa"}' | sh
cd -
basename -a -s '-tRNAscan131.out' `ls $tRNAscanOut/*-tRNAscan131.out` | awk '{print "~/bin/trnascan2gff.pl $tRNAscanOut/"$1"-tRNAscan131.out \t> $tRNAscanOut/"$1"-tRNAscan131.gff"}'  | sh
```

* Missing group I introns:
 * Get alignments from here: http://www.rna.whu.edu.cn/gissd/alignment.html (may need to use the waybackmachine)
 * Li, Z., & Zhang, Y. (2005). Predicting the secondary structures and tertiary interactions of 211 group I introns in IE subgroup. Nucleic acids research, 33(7), 2118-2128.
 
```
cd ~/projects/robustness-RNP/data/genomes-fungi/introns
cat ~/download/view-source_www.rna.whu.edu.cn_gissd_alignment.html | tr "<>" "\n" | grep ".sto.gz" | perl -lane 'if(/alignment\/(\S+.sto.gz)/){print "wget http://www.rna.whu.edu.cn/data/alignment/$1"}' | uniq | sh

ls -1 *sto | tr "." " " | awk '{print "/usr/local/bin/cmbuild "$1".cm "$1".sto"}' | sh     #fix IE3 alignment, missing IE stockholm alignment
ls -1 *cm               | awk '{print "/usr/local/bin/cmcalibrate "$1}'           | sh

ls *cm | egrep -v 'IA1.cm|IE3.cm' | tr "\n" " " | perl -lane 'print "cat $_"' | sh > gssi-grpi-introns.cm
cd ~/projects/robustness-RNP/data/genomes-fungi/
basename -a -s '.fa' `ls *fa ` | awk '{print "~/inst/infernal-1.1/src/cmsearch -o "$1"-gssi-grpi-introns.cmsearch --tblout  "$1"-gssi-grpi-introns.tbl  -E 0.001 --rfam introns/gssi-grpi-introns.cm  "$1".fa"}' | sh
basename -a -s '-gssi-grpi-introns.tbl' `ls $tRNAscanOut/*-gssi-grpi-introns.tbl` | awk '{print "~/bin/cmsearchtblout2gff.pl $tRNAscanOut/"$1"-gssi-grpi-introns.tbl \t> $tRNAscanOut/"$1"-gssi-grpi-introns.gff"}' | sh

basename -a -s '.fa' `ls $tRNAscanOut/*fa` | awk '{print "~/people/bird-genomes/scripts/compete_clans.pl -g "$1"-ncRNAs-with-rfam130.gff -g "$1"-tRNAscan131.gff -g "$1"-gssi-grpi-introns.gff -cl ~/people/bird-genomes/data/clan_info.txt -od ./ > "$1"-ncRNAs.gff"}' | sh 

#counts & join -- find deep/shallow families 
cat   cerevisiae-ncRNAs.gff | perl -lane '$id="unknown"; if(/rfam-id=(\S+);rfam-acc/){$id=$1;}elsif(/type=(\S+);anti/){$id="$1-tRNA"}  print $id;' | sort -d | uniq -c >   cerevisiae-ncRNAs-with-rfam130.counts
cat kudriavzevii-ncRNAs.gff | perl -lane '$id="unknown"; if(/rfam-id=(\S+);rfam-acc/){$id=$1;}elsif(/type=(\S+);anti/){$id="$1-tRNA"}  print $id;' | sort -d | uniq -c > kudriavzevii-ncRNAs-with-rfam130.counts
cat        pombe-ncRNAs.gff | perl -lane '$id="unknown"; if(/rfam-id=(\S+);rfam-acc/){$id=$1;}elsif(/type=(\S+);anti/){$id="$1-tRNA"}  print $id;' | sort -d | uniq -c >        pombe-ncRNAs-with-rfam130.counts
join -a 1 -e 0 -o 1.2 1.1 2.1     -1 2 -2 2 cerevisiae-ncRNAs-with-rfam130.counts  kudriavzevii-ncRNAs-with-rfam130.counts > blah1
join -a 1 -e 0 -o 1.1 1.2 1.3 2.1 -1 1 -2 2 blah1                                         pombe-ncRNAs-with-rfam130.counts         
```

* Annotate the Pfam domains:
```
hmmpress ~/data/pfam/pfam32/Pfam-A.hmm  ##speeds up file I/O

cd ~/projects/robustness-RNP/data/genomes-fungi
basename -a -s '.fa' `ls *fa ` | awk '{print "translate -q "$1".fa > "$1"-6frame-translation.fasta"}' | sh

basename -a -s '.fa' `ls *fa ` | awk '{print "hmmsearch  --cut_ga --noali --domtblout "$1".domtblout ~/data/pfam/pfam32/Pfam-A.hmm "$1"-6frame-translation.fasta > "$1"-pfam32.hmmsearch"}' | sh

#Then run this monstrous one-liner to convert the domtblout file to gff:
cat cerevisiae.domtblout      | perl -lane 'if(/nt\s+(\d+)\.\.(\d+)/){($f,$t)=($1,$2); if($f<$t){$str="+"; $fDNA=$f+$F[17]*3-3; $tDNA=$f+$F[18]*3-3;}else{$str="-"; ($fDNA,$tDNA)=($f-3*$F[18]+1,$f-3*$F[17]+3); } $sid=$F[0]; $sid=~s/\.\d+$//; print "$sid\thmmsearch\tpolypeptide_domain\t$fDNA\t$tDNA\t$F[7]\t$str\t.\tE-value=$F[6];ID=$F[3];ACC=$F[4];hmm-st=$F[15];hmm-en=$F[16];"; }' | sort -k4n > cerevisiae-pfam32.gff
cat kudriavzevii.domtblout    | perl -lane 'if(/nt\s+(\d+)\.\.(\d+)/){($f,$t)=($1,$2); if($f<$t){$str="+"; $fDNA=$f+$F[17]*3-3; $tDNA=$f+$F[18]*3-3;}else{$str="-"; ($fDNA,$tDNA)=($f-3*$F[18]+1,$f-3*$F[17]+3); } $sid=$F[0]; $sid=~s/\.\d+$//; print "$sid\thmmsearch\tpolypeptide_domain\t$fDNA\t$tDNA\t$F[7]\t$str\t.\tE-value=$F[6];ID=$F[3];ACC=$F[4];hmm-st=$F[15];hmm-en=$F[16];"; }' | sort -k4n > kudriavzevii-pfam32.gff
cat pombe.domtblout           | perl -lane 'if(/nt\s+(\d+)\.\.(\d+)/){($f,$t)=($1,$2); if($f<$t){$str="+"; $fDNA=$f+$F[17]*3-3; $tDNA=$f+$F[18]*3-3;}else{$str="-"; ($fDNA,$tDNA)=($f-3*$F[18]+1,$f-3*$F[17]+3); } $sid=$F[0]; $sid=~s/\.\d+$//; print "$sid\thmmsearch\tpolypeptide_domain\t$fDNA\t$tDNA\t$F[7]\t$str\t.\tE-value=$F[6];ID=$F[3];ACC=$F[4];hmm-st=$F[15];hmm-en=$F[16];"; }' | sort -k4n > pombe-pfam32.gff

cat   cerevisiae-pfam32.gff | perl -lane '$id="unknown"; if(/;ID=(\S+);ACC/){$id=$1;} print $id;' | sort -d | uniq -c >   cerevisiae-pfam32.counts
cat kudriavzevii-pfam32.gff | perl -lane '$id="unknown"; if(/;ID=(\S+);ACC/){$id=$1;} print $id;' | sort -d | uniq -c > kudriavzevii-pfam32.counts
cat        pombe-pfam32.gff | perl -lane '$id="unknown"; if(/;ID=(\S+);ACC/){$id=$1;} print $id;' | sort -d | uniq -c >        pombe-pfam32.counts
join -a 1 -e 0 -o 1.2 1.1 2.1     -1 2 -2 2 cerevisiae-pfam32.counts  kudriavzevii-pfam32.counts > blah1
join -a 1 -e 0 -o 1.1 1.2 1.3 2.1 -1 1 -2 2 blah1                            pombe-pfam32.counts         
```

* Generate reverse translated mRNA alignments: 
```
cat deep-families.tsv    | cut -f 2 | sort -d | uniq | awk '{print "~/projects/robustness-RNP/bin/fetchPfam.sh "$1" cerevisiae-pfam32.gff  cerevisiae.fa        pombe-pfam32.gff         pombe.fa "}' | sh
cat shallow-families.tsv | cut -f 2 | sort -d | uniq | awk '{print "~/projects/robustness-RNP/bin/fetchPfam.sh "$1" cerevisiae-pfam32.gff  cerevisiae.fa kudriavzevii-pfam32.gff  kudriavzevii.fa "}' | sh
ls -s mrna-seqs/*.nuc.stk
```

* Generate ncRNA alignments: 
```
cat deep-families.tsv    | cut -f 1 | sort -d | uniq | awk '{print "~/projects/robustness-RNP/bin/fetchRfam.sh "$1" cerevisiae-ncRNAs.gff  cerevisiae.fa        pombe-ncRNAs.gff         pombe.fa"}' | sh
cat shallow-families.tsv | cut -f 1 | sort -d | uniq | awk '{print "~/projects/robustness-RNP/bin/fetchRfam.sh "$1" cerevisiae-ncRNAs.gff  cerevisiae.fa kudriavzevii-ncRNAs.gff  kudriavzevii.fa"}' | sh
#The Rfam Group I Intron missed this:
~/inst/infernal-1.1.1/src/cmalign  introns/IC2.cm   ncrna-seqs/Intron_gpI-IC2.fasta > ncrna-seqs/blah && ~/inst/infernal-1.1.1/easel/miniapps/esl-reformat pfam ncrna-seqs/blah > ncrna-seqs/Intron_gpI-IC2.stk

ls -1 ncrna-seqs/*.stk     | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 1 -i "$1}' | sh > synonNonsynon-ncRNA-results.txt
ls -1  mrna-seqs/*.nuc.stk | awk '{print "echo "$1"&& ../../bin/computeSynonNonsynon.pl -r 0 -i "$1}' | sh > synonNonsynon-mRNA-results.txt

# Some hoop-jumping to get into tabular format:
echo -e "ID\tdistID\ttotVar\tlen\taaSynTot\taaSynProp\taaNonSynTot\taaNonSynProp\tbanpSynTot\tbanpSynProp\tbanpNonSynTot\tbanpNonSynProp\tblossSynTot\tblossSynProp\tblossNonSynTot\tblossNonSynProp" > synonNonsynon-mRNA-results.tsv
cat synonNonsynon-mRNA-results.txt | perl -lane 'if(  /^mrna-seqs\/(\S+)\.nuc/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/aaSynon:\s+(\d+) \((\S+)\)\s+aaNonsynon:\s+(\d+) \((\S+)\)/){($aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP)=($1,$2,$3,$4);}elsif(/banpSynon:\s+(\d+) \((\S+)\)\s+banpNonsynon:\s+(\d+) \((\S+)\)/){($banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP)=($1,$2,$3,$4);}elsif(/blossumSynon:\s+(\d+) \((\S+)\)\s+blossumNonsynon:\s+(\d+) \((\S+)\)/){($blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP)=($1,$2,$3,$4); $dist='deep'; $dist='shallow' if (`grep $id shallow-families.tsv`);   @a=($id,$dist,$tV,$len,$aaSynT,$aaSynP,$aaNonSynT,$aaNonSynP,$banpSynT,$banpSynP,$banpNonSynT,$banpNonSynP,$blossSynT,$blossSynP,$blossNonSynT,$blossNonSynP); $a=join("\t",@a); print $a}' >> synonNonsynon-mRNA-results.tsv

echo -e "ID\tdistID\ttotVar\tlen\trySynTot\trySynProp\tryNonSynTot\tryNonSynProp\tbpSynTot\tbpSynProp\tbpNonSynTot\tbpNonSynProp" > synonNonsynon-ncRNA-results.tsv
cat synonNonsynon-ncRNA-results.txt | perl -lane 'if(/^ncrna-seqs\/(\S+)\.stk/){$id=$1}elsif(/Total variation:\s+(\d+)\s+len=(\d+)/){($tV,$len)=($1,$2)}elsif(/rySynon:\s+(\d+) \((\S+)\)\s+ryNonsynon:\s+(\d+) \((\S+)\)/){($rySynT,$rySynP,$ryNonSynT,$ryNonSynP)=($1,$2,$3,$4);}elsif(/bpSynon:\s+(\d+) \((\S+)\)\s+bpNonsynon:\s+(\d+) \((\S+)\)/){($bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP)=($1,$2,$3,$4); $dist='deep'; $dist='shallow' if (`grep $id shallow-families.tsv`); @a=($id,$dist,$tV,$len,$rySynT,$rySynP,$ryNonSynT,$ryNonSynP,$bpSynT,$bpSynP,$bpNonSynT,$bpNonSynP); $a=join("\t",@a); print $a }' >> synonNonsynon-ncRNA-results.tsv
```

* Curate functional linkages between the protein and RNA families  
 * MUCH EASIER TO FIND LINKS WITH THIS:
 * https://www.yeastgenome.org/search?q=ribonucleoprotein&category=complex
 * e.g. https://www.yeastgenome.org/complex/CPX-32
 * Ribosome: https://www.annualreviews.org/doi/10.1146/annurev-biochem-060614-033917

```
Table 1

#Deep
5_8S_rRNA               Ribosomal_L7Ae
5_8S_rRNA               Ribosomal_L19e
5S_rRNA                 Ribosomal_L5e
5S_rRNA                 Brix
5S_rRNA                 SART-1
SSU_rRNA_eukarya        S4
SSU_rRNA_eukarya        Ribosomal_S5
SSU_rRNA_eukarya        Ribosomal_S6e
SSU_rRNA_eukarya        Ribosomal_S7e
SSU_rRNA_eukarya        Ribosomal_S8
SSU_rRNA_eukarya        Ribosomal_S8e
SSU_rRNA_eukarya        Ribosomal_S9
SSU_rRNA_eukarya        Ribosomal_S10
SSU_rRNA_eukarya        S10_plectin
SSU_rRNA_eukarya        Ribosomal_S11
SSU_rRNA_eukarya        Ribosom_S12_S23
SSU_rRNA_eukarya        Ribosomal_S13
SSU_rRNA_eukarya        Ribosomal_S14
SSU_rRNA_eukarya        Ribosomal_S15
SSU_rRNA_eukarya        Ribosomal_S17_N
SSU_rRNA_eukarya        Ribosomal_S17
SSU_rRNA_eukarya        Ribosomal_S17e
SSU_rRNA_eukarya        Ribosomal_S19
SSU_rRNA_eukarya        Ribosomal_S19e
SSU_rRNA_eukarya        Ribosomal_S21e
LSU_rRNA_eukarya        Ribosomal_L4
LSU_rRNA_eukarya        Ribosomal_L5_C
LSU_rRNA_eukarya        Ribosomal_L6
LSU_rRNA_eukarya        Ribosomal_L10
LSU_rRNA_eukarya        Ribosomal_L13
LSU_rRNA_eukarya        Ribosomal_L14
LSU_rRNA_eukarya        Ribosomal_L16
LSU_rRNA_eukarya        Ribosomal_L22
LSU_rRNA_eukarya        Ribosomal_L23
LSU_rRNA_eukarya        Ribosomal_L26
LSU_rRNA_eukarya        KOW
LSU_rRNA_eukarya        Ribosomal_L27A
LSU_rRNA_eukarya        Ribosomal_L29
LSU_rRNA_eukarya        Ribosomal_L30
LSU_rRNA_eukarya        Ribosomal_60s
RNase_MRP               POP1	#https://www.yeastgenome.org/complex/CPX-3284
RNase_MRP               POPLD
RNase_MRP               RNase_P_pop3
RNase_MRP               UPF0086
RNase_MRP               RNase_P_Rpp14
RNase_MRP               Rpp20
RNase_MRP               RNase_P_p30
RNase_MRP               Rpr2
RNaseP_nuc              POP1      #https://www.yeastgenome.org/complex/CPX-1294
RNaseP_nuc              Rpp20
RNaseP_nuc              RNase_P_Rpp14                
U2                      DUF382
U2                      PSP
U2                      RRM_5
U2                      LRR_9
U2                      FHA
U2                      zf-met
U2                      SF3A2
U2                      Surp
U2                      SF3A3
U2                      zf-C2H2_jaz
U2                      PHF5
U2                      MMS1_N
U2                      CPSF_A
U2                      LSM
U2                      SF3b10
U4                      LSM	
U4                      PRP3
U4                      DUF1115
U4                      Nop
U4                      Prp31_C
U5                      LSM
U6                      LSM
U6                      PRP3
U6                      DUF1115
U6                      Nop
U6                      Prp31_C
Fungi_U3                NOP5NT
Fungi_U3                Ribosomal_S4
Fungi_U3                S4
Fungi_U3                Brix
Asn-tRNA                tRNA_anti-codon  #tRNAs  -- low copy number
Asn-tRNA                tRNA-synt_2
Cys-tRNA                tRNA-synt_1e
Gln-tRNA                tRNA-synt_1c
His-tRNA                tRNA-synt_His
Met-tRNA                tRNA-synt_1g
Phe-tRNA                tRNA-synt_2d
Pro-tRNA                tRNA-synt_2b
Trp-tRNA                tRNA-synt_1b
Tyr-tRNA                tRNA-synt_1b


#Shallow
Sacc_telomerase		EST1
Sacc_telomerase		EST1_DNA_bind
Sacc_telomerase		Telomerase_RBD
Sacc_telomerase		TPP1
srg1                    2-Hacid_dh
srg1                    2-Hacid_dh_C
Intron_gpII             COX1
Intron_gpII             LAGLIDADG_2
Intron_gpI-IC2          Cytochrome_B
Intron_gpI-IC2          LAGLIDADG_1
Intron_gpI-IC2          Ydc2-catalyt
U1_yeast                RRM_1
U1_yeast                WW
U1_yeast                LUC7
U1_yeast                zf-U1
U1_yeast                LSM
Fungi_SRP               SRP54
Fungi_SRP               SRP54_N
Fungi_SRP               SRP_SPB
Fungi_SRP               SRP19
Fungi_SRP               SRP14
Fungi_SRP               SRP9-21
Fungi_SRP               SRP68
Fungi_SRP               SRP72
snR4                    Nop             #C/D box 
snR40                   Nop
snR41                   Nop
snR45                   Nop
snR50                   Nop
snR4                    NOP5NT
snR40                   NOP5NT
snR41                   NOP5NT
snR45                   NOP5NT
snR50                   NOP5NT
snR3                    Nop10p    #H/ACA
snR5                    Nop10p
snR8                    Nop10p
snR9                    Nop10p
snR10                   Nop10p
snR3                    SHQ1
snR5                    SHQ1
snR8                    SHQ1
snR9                    SHQ1
snR10                   SHQ1
snR3                    TruB_N
snR5                    TruB_N
snR8                    TruB_N
snR9                    TruB_N
snR10                   TruB_N
```

* Compute a phylogenetic tree

```
#The kudriavzevii genome is missing the rRNAs, hence fetching seqs from RNA central:
>S.cerevisiae            Saccharomyces cerevisiae rRNA SSU 18S   https://rnacentral.org/rna/URS00005F2C2D/4932
UAUCUGGUUGAUCCUGCCAGUAGUCAUAUGCUUGUCUCAAAGAUUAAGCCAUGCAUGUCUAAGUAUAAGCAAUUUAUACAGUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUCGUUUAUUUGAUAGUUCCUUUACUACAUGGUAUAACUGUGGUAAUUCUAGAGCUAAUACAUGCUUAAAAUCUCGACCCUUUGGAAGAGAUGUAUUUAUUAGAUAAAAAAUCAAUGUCUUCGGACUCUUUGAUGAUUCAUAAUAACUUUUCGAAUCGCAUGGCCUUGUGCUGGCGAUGGUUCAUUCAAAUUUCUGCCCUAUCAACUUUCGAUGGUAGGAUAGUGGCCUACCAUGGUUUCAACGGGUAACGGGGAAUAAGGGUUCGAUUCCGGAGAGGGAGCCUGAGAAACGGCUACCACAUCCAAGGAAGGCAGCAGGCGCGCAAAUUACCCAAUCCUAAUUCAGGGAGGUAGUGACAAUAAAUAACGAUACAGGGCCCAUUCGGGUCUUGUAAUUGGAAUGAGUACAAUGUAAAUACCUUAACGAGGAACAAUUGGAGGGCAAGUCUGGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGUUGCAGUUAAAAAGCUCGUAGUUGAACUUUGGGCCCGGUUGGCCGGUCCGAUUUUUUCGUGUACUGGAUUUCCAACGGGGCCUUUCCUUCUGGCUAACCUUGAGUCCUUGUGGCUCUUGGCGAACCAGGACUUUUACUUUGAAAAAAUUAGAGUGUUCAAAGCAGGCGUAUUGCUCGAAUAUAUUAGCAUGGAAUAAUAGAAUAGGACGUUUGGUUCUAUUUUGUUGGUUUCUAGGACCAUCGUAAUGAUUAAUAGGGACGGUCGGGGGCAUCAGUAUUCAAUUGUCAGAGGUGAAAUUCUUGGAUUUAUUGAAGACUAACUACUGCGAAAGCAUUUGCCAAGGACGUUUUCAUUAAUCAAGAACGAAAGUUAGGGGAUCGAAGAUGAUCAGAUACCGUCGUAGUCUUAACCAUAAACUAUGCCGACUAGGGAUCGGGUGGUGUUUUUUUAAUGACCCACUCGGCACCUUACGAGAAAUCAAAGUCUUUGGGUUCUGGGGGGAGUAUGGUCGCAAGGCUGAAACUUAAAGGAAUUGACGGAAGGGCACCACCAGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACGGGGAAACUCACCAGGUCCAGACACAAUAAGGAUUGACAGAUUGAGAGCUCUUUCUUGAUUUUGUGGGUGGUGGUGCAUGGCCGUUCUUAGUUGGUGGAGUGAUUUGUCUGCUUAAUUGCGAUAACGAACGAGACCUUAACCUACUAAAUAGUGGUGCUAGCAUUUGCUGGUUAUCCACUUCUUAGAGGGACUAUCGGUUUCAAGCCGAUGGAAGUUUGAGGCAAUAACAGGUCUGUGAUGCCCUUAGACGUUCUGGGCCGCACGCGCGCUACACUGACGGAGCCAGCGAGUCUAACCUUGGCCGAGAGGUCUUGGUAAUCUUGUGAAACUCCGUCGUGCUGGGGAUAGAGCAUUGUAAUUAUUGCUCUUCAACGAGGAAUUCCUAGUAAGCGCAAGUCAUCAGCUUGCGUUGAUUACGUCCCUGCCCUUUGUACACACCGCCCGUCGCUAGUACCGAUUGAAUGGCUUAGUGAGGCCUCAGGAUCUGCUUAGAGAAGGGGGCAACUCCAUCUCAGAGCGGAGAAUUUGGACAAACUUGGUCAUUUAGAGGAACUAAAAGUCGUAACAAGGUUUCCGUAGGUGAACCUGCGGAAGGAUCAUUA
>S.kudriavzevii          Saccharomyces kudriavzevii    https://rnacentral.org/rna/URS00001926EE/114524
CUGGUUGAUCCUGCCAGUAGUCAUAUGCUUGUCUCAAAGAUUAAGCCAUGCAUGUCUAAGUAUAAGCAAUUUAUACAGUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUCGUUUAUUUGAUAGUUCCUUUACUACAUGGUAUAACUGUGGUAAUUCUAGAGCUAAUACAUGCUUAAAAUCUCGACCCUUUGGAAGAGAUGUAUUUAUUAGAUAAAAAAUCAAUGUCUUCGGACUCUUUGAUGAUUCAUAAUAACUUUUCGAAUCGCAUGGCCUUGUGCUGGCGAUGGUUCAUUCAAAUUUCUGCCCUAUCAACUUUCGAUGGUAGGAUAGUGGCCUACCAUGGUUUCAACGGGUAACGGGGAAUAAGGGUUCGAUUCCGGAGAGGGAGCCUGAGAAACGGCUACCACAUCCAAGGAAGGCAGCAGGCGCGCAAAUUACCCAAUCCUAAUUCAGGGAGGUAGUGACAAUAAAUAACGAUACAGGGCCCAUUCGGGUCUUGUAAUUGGAAUGAGUACAAUGUAAAUACCUUAACGAGGAACAAUUGGAGGGCAAGUCUGGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGUUGCAGUUAAAAAGCUCGUAGUUGAACUUUGGGCUCGGUUGGCCGGUCCGAUUUUUUCGUGUACUGGAUUUCCAACGGGGCCUUUCCUUCUGGCUAACCUUGAGUCCUUGUGGCUCUUGGCGAACCAGGACUUUUACUUUGAAAAAAUUAGAGUGUUCAAAGCAGGCGUAUUGCUCGAAUAUAUUAGCAUGGAAUAAUAGAAUAGGACGUUUGGUUCUAUUUUGUUGGUUUCUAGGACCAUCGUAAUGAUUAAUAGGGACGGUCGGGGGCAUCAGUAUUCAAUUGUCAGAGGUGAAAUUCUUGGAUUUAUUGAAGACUAACUACUGCGAAAGCAUUUGCCAAGGACGUUUUCAUUAAUCAAGAACGAAAGUUAGGGGAUCGAAGAUGAUCAGAUACCGUCGUAGUCUUAACCAUAAACUAUGCCGACUAGGGAUCGGGUGGUGUUUUUUUAAUGACCCACUCGGCACCUUACGAGAAAUCAAAGUCUUUGGGUUCUGGGGGGAGUAUGGUCGCAAGGCUGAAACUUAAAGGAAUUGACGGAAGGGCACCACCAGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACGGGGAAACUCACCAGGUCCAGACACAAUAAGGAUUGACAGAUUGAGAGCUCUUUCUUGAUUUUGUGGGUGGUGGUGCAUGGCCGUUCUUAGUUGGUGGAGUGAUUUGUCUGCUUAAUUGCGAUAACGAACGAGACCUUAACCUACUAAAUAGUGGUGCUAGCAUUUGCUGGUUAUCCACUUCUUAGAGGGACUAUCGGUUUCAAGCCGAUGGAAGUUUGAGGCAAUAACAGGUCUGUGAUGCCCUUAGACGUUCUGGGCCGCACGCGCGCUACACUGACGGAGCCAGCGAGUCUAACCUUGGCCGAGAGGUCUUGGUAAUCUUGUGAAACUCCGUCGUGCUGGGGAUAGAGCAUUGUAAUUAUUGCUCUUCAACGAGGAAUUCCUAGUAAGCGCAAGUCAUCAGCUUGCGUUGAUUACGUCCCUGCCCUUUGUACACACCGCCCGUCGCUAGUACCGAUUGAAUGGCUUAGUGAGGCCUCAGGAUCUGCUUAGAGAAGGGGGCAACUCCAUCUCAGAGCGGAGAAUUUGGACAAACUUGGUCAUUUAGAGGAACUAAAAGUCGUAACAAGGUUUCC
>S.pombe                 Schizosaccharomyces pombe (fission yeast) 18S ribosomal RNA    https://rnacentral.org/rna/URS000055C986/4896
UACCUGGUUGAUCCUGCCAGUAGUCAUAUGCUUGUCUCAAAGAUUAAGCCAUGCAUGUCUAAGUAUAAGCAAUUUUGUACUGUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUCGUUUAUUUGAUAGUACCUCAACUACUUGGAUAACCGUGGUAAUUCUAGAGCUAAUACAUGCUAAAAAUCCCGACUUUUUUGGAAGGGAUGUAUUUAUUAGAUAAAAAACCAAUGCCUUCGGGCUUUUUUUGGUGAGUCAUAAUAACUUUUCGAAUCGCAUGGCCUUGCGCCGGCGAUGGUUCAUUCAAAUUUCUGCCCUAUCAACUUUCGAUGGUAGGAUAGAGGCCUACCAUGGUUUUAACGGGUAACGGGGAAUUAGGGUUCGAUUCCGGAGAGGGAGCCUGAGAAACGGCUACCACAUCCAAGGAAGGCAGCAGGCGCGCAAAUUACCCAAUCCCGACACGGGGAGGUAGUGACAAGAAAUAACAAUGCAGGGCCCUUUCGGGUCUUGUAAUUGGAAUGAGUACAAUGUAAAUACCUUAACGAGGAACAAUUGGAGGGCAAGUCUGGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGUUGCAGUUAAAAAGCUCGUAGUUGAACUUUGAGCCUGGUCGACUGGUCCGCCGCAAGGCGUGUUUACUGGUCAUGACCGGGGUCGUUAACCUUCUGGCAAACUACUCAUGUUCUUUAUUGAGCGUGGUAGGGAACCAGGACUUUUACCUUGAAAAAAUUAGAGUGUUCAAAGCAGGCAAGUUUUGCUCGAAUACAUUAGCAUGGAAUAAUAAAAUAGGACGUGUGGUUCUAUUUUGUUGGUUUCUAGGACCGCCGUAAUGAUUAAUAGGGAUAGUCGGGGGCAUUCGUAUUCAAUUGUCAGAGGUGAAAUUCUUGGAUUUAUUGAAGACGAACUACUGCGAAAGCAUUUGCCAAGGAUGUUUUCAUUAAUCAAGAACGAAAGUUAGGGGAUCGAAGACGAUCAGAUACCGUCGUAGUCUUAACCAUAAACUAUGCCGACUAGGGAUCGGGCAAUGUUUCAUUUAUCGACUUGCUCGGCACCUUACGAGAAAUCAAAGUCUUUGGGUUCCGGGGGGAGUAUGGUCGCAAGGCUGAAACUUAAAGGAAUUGACGGAAGGGCACCACAAUGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACGGGGAAACUCACCAGGUCCAGACAUAGUAAGGAUUGACAGAUUGAGAGCUCUUUCUUGAUUCUAUGGGUGGUGGUGCAUGGCCGUUCUUAGUUGGUGGAGUGAUUUGUCUGCUUAAUUGCGAUAACGAACGAGACCUUAACCUGCUAAAUAGCUGGAUCAGCCAUUUUGGCUGAUCAUUAGCUUCUUAGAGGGACUAUUGGCAUAAAGCCAAUGGAAGUUUGAGGCAAUAACAGGUCUGUGAUGCCCUUAGAUGUUCUGGGCCGCACGCGCGCUACACUGACGGAGCCAACGAGUUGAAAAAAAUCUUUUGAUUUUUUAUCCUUGGCCGGAAGGUCUGGGUAAUCUUGUUAAACUCCGUCGUGCUGGGGAUAGAGCAUUGCAAUUAUUGCUCUUCAACGAGGAAUUCCUAGUAAGCGCAAGUCAUCAGCUUGCGUUGAAUACGUCCCUGCCCUUUGUACACACCGCCCGUCGCUACUACCGAUUGAAUGGCUUAGUGAGGCCUCUGGAUUGGCUUGUUUCUGCUGGCAACGGCGGAAACAUUGCCGAGAAGUUGGACAAACUUGGUCAUUUAGAGGAAGUAAAAGUCGUAACAAGGUUUCCGUAGGUGAACCUGCGGAAGGAUCAUUA

cd data/phylogeny
~/inst/infernal-1.1/src/cmfetch ~/data/rfam/rfam14.0/Rfam.cm  SSU_rRNA_eukarya > SSU_rRNA_eukarya.cm

~/inst/infernal-1.1/src/cmalign --outformat PHYLIP -g -o fungi-ssu.phy  SSU_rRNA_eukarya.cm fungi-ssu.fasta
cp fungi-ssu.phy   infile
dnaml
mv outtree fungi-ssu.dnd
figtree   fungi-ssu.dnd &

#handy tab files for Rfam & Pfam IDs/ACCs:
cat ~/data/rfam/rfam14.1/Rfam.cm | perl -lane 'if(/^INFERNAL1\/a/){$p=1}elsif(/^HMMER3\/f/){$p=0}  if($p && /^NAME\s+(\S+)/){$i=$1}elsif($p && /^ACC\s+(RF0\d+)/){$a=$1}elsif($p && /DESC\s+(.*)/){print "$i\t$a\t$1"; $i=""; $a=""; $p=0;}' > ../rfam-id-acc-desc.tsv
cat ~/data/pfam/pfam32/Pfam-A.hmm |  perl -lane 'if(/^NAME\s+(\S+)/){$i=$1}elsif(/^ACC\s+(PF\d+)/){$a=$1}elsif(/DESC\s+(.*)/){print "$i\t$a\t$1"; $i=""; $a="";}' > ../pfam-id-acc-desc.tsv
```
