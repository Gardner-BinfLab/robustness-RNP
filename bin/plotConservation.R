#!/usr/bin/Rscript

mRNA  <- read.table("synonNonsynon-mRNA-results.tsv",  sep="\t", header=TRUE)
ncRNA <- read.table("synonNonsynon-ncRNA-results.tsv", sep="\t", header=TRUE)

mRNA.shallow  <- mRNA[mRNA$distID == "ecol2styp",]$totVar   / mRNA[mRNA$distID == "ecol2styp",]$len
mRNA.deep     <- mRNA[mRNA$distID == "ecol2nmen",]$totVar   / mRNA[mRNA$distID == "ecol2nmen",]$len
ncRNA.shallow <- ncRNA[ncRNA$distID == "ecol2styp",]$totVar / ncRNA[ncRNA$distID == "ecol2styp",]$len
ncRNA.deep    <- ncRNA[ncRNA$distID == "ecol2nmen",]$totVar / ncRNA[ncRNA$distID == "ecol2nmen",]$len

len.mRNA.shallow  <- length(mRNA.shallow)
len.mRNA.deep     <- length(mRNA.deep)
len.ncRNA.shallow <- length(ncRNA.shallow)
len.ncRNA.deep    <- length(ncRNA.deep)

maxRows <- max(len.mRNA.shallow,len.mRNA.deep,len.ncRNA.shallow,len.ncRNA.deep)
conservationRP <- matrix(NA,ncol=4,nrow=maxRows)
colnames(conservationRP)<-c("Shallow \n(RNA)", "Deep\n(RNA)", "Shallow\n(Protein)", "Deep\n(Protein)")

conservationRP[1:len.ncRNA.shallow,1] <- ncRNA.shallow
conservationRP[1:len.mRNA.shallow, 2] <- mRNA.shallow
conservationRP[1:len.ncRNA.deep,   3] <- ncRNA.deep
conservationRP[1:len.mRNA.deep,    4] <- mRNA.deep

pdf(file="../../manuscript/figures/figure1a.pdf", width=12, height=12)
par(cex=2.5,las=2)
col <- c("orchid1", "skyblue", "lightpink", "lightblue1" )
boxplot(conservationRP, col=col, ylab="Nucleotide variation", xlab="", main="Nucleotide conservation of RNA & protein", outline=FALSE, xaxt = "n")
axis(1, 1:4, colnames(conservationRP))
cols<-c("deeppink4","dodgerblue4","deeppink4","dodgerblue4")
for (i in 1:4){
    xp <-  rnorm(length(conservationRP[,i]), mean = 0, sd = 0.1)   
    points(i+xp, conservationRP[,i],col=cols[i], pch=20, cex=0.75)
}
legend("bottomright", c("RNA", "", "Protein", ""), col=col[c(1,3,2,4)], fil=col[c(1,3,2,4)],cex=0.6, ncol = 2)
dev.off()

#system("convert -flatten conservation-new-old-genes.pdf conservation-new-old-genes.png")

wilcox.test(conservationRP[,1], conservationRP[,2])
wilcox.test(conservationRP[,3], conservationRP[,4])

#for the tree:
col2rgb("orchid1")
col2rgb("lightpink")
col2rgb("skyblue")
col2rgb("lightblue1")

#SHALLOW
# > col2rgb("orchid1")
#       [,1]
# red    255
# green  131
# blue   250
# > col2rgb("skyblue")
#       [,1]
# red    135
# green  206
# blue   235

#DEEP
# > col2rgb("lightpink")
#       [,1]
# red    255
# green  182
# blue   193
# > col2rgb("lightblue1")
#       [,1]
# red    191
# green  239
# blue   255
######################################################################
#Plot fraction of synonymous mutations for each group shallow/deep & RNA/protein

synonRP <- matrix(NA,ncol=10,nrow=maxRows)
synonRP[1:len.mRNA.shallow,  1] <- mRNA[mRNA$distID   == "ecol2styp",]$aaSynProp
synonRP[1:len.mRNA.deep,     2] <- mRNA[mRNA$distID   == "ecol2nmen",]$aaSynProp
synonRP[1:len.ncRNA.shallow, 3] <- ncRNA[ncRNA$distID == "ecol2styp",]$rySynProp
synonRP[1:len.ncRNA.deep,    4] <- ncRNA[ncRNA$distID == "ecol2nmen",]$rySynProp
synonRP[1:len.mRNA.shallow,  5] <- mRNA[mRNA$distID   == "ecol2styp",]$banpSynProp
synonRP[1:len.mRNA.deep,     6] <- mRNA[mRNA$distID   == "ecol2nmen",]$banpSynProp
synonRP[1:len.mRNA.shallow,  7] <- mRNA[mRNA$distID   == "ecol2styp",]$blossSynProp
synonRP[1:len.mRNA.deep,     8] <- mRNA[mRNA$distID   == "ecol2nmen",]$blossSynProp
synonRP[1:len.ncRNA.shallow, 9] <- ncRNA[ncRNA$distID == "ecol2styp",]$bpSynProp
synonRP[1:len.ncRNA.deep,   10] <- ncRNA[ncRNA$distID == "ecol2nmen",]$bpSynProp

pdf(file="../../manuscript/figures/suppfigure1.pdf", width=16, height=10)
par(cex=2.5,las=2, mar = c(7,4,2,2) + .1, tck=-0.01) #c(bottom, left, top, right).  c(5, 4, 4, 2)
cols <- c(c("skyblue", "lightblue1", "orchid1", "lightpink", "skyblue", "lightblue1", "skyblue", "lightblue1", "orchid1", "lightpink" ))
boxplot(synonRP, col=cols, ylab="Proportion of (near) neutral mutations", xlab="", main="", outline=FALSE, xaxt = "n")
names <- c("Degeneracy", "Biochemistry\n(R<->Y)", "Biochemistry\n(BANP)", "BLOSUM", "Secondary\nstructure" )
axis(1, 2*(1:length(names))-0.5, names ) #
cols<-c("dodgerblue4","dodgerblue4","deeppink4","deeppink4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","deeppink4","deeppink4")
for (i in 1:length(synonRP[1,])){
    #ii <- order[i] #sortedMedians$ix[i]
    xp <-  rnorm(length(synonRP[,i]), mean = 0, sd = 0.1)   
    points(i+xp, synonRP[,i],col=cols[i], pch=20, cex=0.75)
}
legend("bottomright", c("Shallow RNA", "Deep RNA", "Shallow protein", "Deep protein"), col=c("orchid1", "lightpink", "skyblue", "lightblue1"), fil=c("orchid1", "lightpink", "skyblue", "lightblue1"),cex=0.65)
dev.off()

######################################################################

# 5S_rRNA	Ribosomal_L5_C
# 5S_rRNA	Ribosomal_L5
# 6S		Sigma70_ner
# Arg		Arg_tRNA_synt_N
# Arg		tRNA-synt_1d
# Bacteria_small_SRP	SRP54_N
# Bacteria_small_SRP	SRP54
# Bacteria_small_SRP	SRP_SPB
# Cobalamin		TonB_dep_Rec
# Cobalamin		Plug
# cspA			CSD
# CsrB			CsrA
# FMN			DHBP_synthase
# GcvB			Hfq
# GlmY_tke1		ATP_bind_2
# GlmZ_SraJ		RNase_E_G
# His_leader		His_leader
# Leu_leader		Leu_leader
# Leu			tRNA-synt_1
# Leu			tRNA-synt_1_2
# LSU_rRNA_bacteria	Ribosomal_L6
# LSU_rRNA_bacteria	ribosomal_L24
# LSU_rRNA_bacteria	Ribosomal_L31
# Mg_sensor		MGTL
# MicA			Hfq
# OxyS			Hfq
# Phe			tRNA_bind
# Phe			B3_4
# RNaseP_bact_a		Ribonuclease_P
# ROSE_2		HSP20
# rseX			Hfq
# S15			Ribosomal_S15
# Ser			tRNA-synt_2b
# Ser			Seryl_tRNA_N
# SgrS			SgrT
# SSU_rRNA_bacteria	Ribosomal_S15
# SSU_rRNA_bacteria	Ribosomal_S17
# SSU_rRNA_bacteria	Ribosomal_S21
# Thr_leader		Leader_Thr
# tmRNA			S1
# tmRNA			GTP_EFTU
# TPP			ThiC-associated
# TPP			ThiC_Rad_SAM
# Trp_leader		Leader_Trp

pairs  <- read.table("rna-protein-pairs.tsv",  sep="\t", header=TRUE)
percentVar           <- matrix(NA,ncol=2,nrow=length(pairs$RNA))
colnames(percentVar) <- c("RNA.var", "protein.var")
for(i in 1:length(pairs$RNA)){
      percentVar[i,1] <- 100*ncRNA[ncRNA$ID == as.character(pairs[i,]$RNA),    ]$totVar / ncRNA[ncRNA$ID == as.character(pairs[i,]$RNA),    ]$len      
      percentVar[i,2] <- 100*mRNA[  mRNA$ID == as.character(pairs[i,]$protein),]$totVar / mRNA[  mRNA$ID == as.character(pairs[i,]$protein),]$len            
}


pdf(file="../../manuscript/figures/suppfigure2.pdf", width=11, height=10)
par(cex=2.5,las=2, mar = c(5,4,0,2) + .1, tck=-0.01) #c(bottom, left, top, right).  c(5, 4, 4, 2)
plot(NA, NA, ylim=c(0,90), xlim=c(0,90), ylab="mRNA variation (%)", xlab="ncRNA variation (%)")
#plot lines first
for(i in 1:(length(pairs$RNA)-1)){
      for(j in (i+1):(length(pairs$RNA)-1)){
	    col <- "purple"; 
      	    if(pairs[i,]$phylogroup == 'ecol2nmen'){
      	 	col <- "violet"; 
      	    }
      	    #shared RNA family:
      	    if(pairs$RNA[i] == pairs$RNA[j]){
		lines(percentVar[c(i,j),1],percentVar[c(i,j),2],col=col,lwd=2)
	    }
	    #shared protein family:
      	    if(pairs$protein[i] == pairs$protein[j]){
		lines(percentVar[c(i,j),1],percentVar[c(i,j),2],col=col,lwd=2)
	    }
      }
}
#plot points
col <- "white"; 
for(i in 1:length(pairs$RNA)){
      pch <- 21; bg <- "purple"; 
      if(pairs[i,]$phylogroup == 'ecol2nmen'){
      	 pch <- 22; bg <- "violet"; 
      }
      points(percentVar[i,1], percentVar[i,2], pch=pch, col=col, bg = bg)
      pch2 <- switch(as.character(pairs[i,]$type),
        "RNP"             = 3,
     	"Cis-regulatory"  = 4,
      	"tRNA+synthetase" = 8, 
      	"Ribosome"        = 11,
      	"Dual-function"   = 14
      )
      points(percentVar[i,1], percentVar[i,2], pch=pch2, col="black", cex = 0.5)
}
d <- 2
rect(
min(percentVar[pairs$phylogroup == 'ecol2styp',1]) - d,
min(percentVar[pairs$phylogroup == 'ecol2styp',2]) - d,
max(percentVar[pairs$phylogroup == 'ecol2styp',1]) + d,
max(percentVar[pairs$phylogroup == 'ecol2styp',2]) + d
)
rect(
min(percentVar[pairs$phylogroup == 'ecol2nmen',1]) - d,
min(percentVar[pairs$phylogroup == 'ecol2nmen',2]) - d,
max(percentVar[pairs$phylogroup == 'ecol2nmen',1]) + d,
max(percentVar[pairs$phylogroup == 'ecol2nmen',2]) + d
)
x <- as.vector(percentVar[pairs$phylogroup == 'ecol2styp',1])
y <- as.vector(percentVar[pairs$phylogroup == 'ecol2styp',2])
lmObj <- lm( y ~ x )
xs <- range(x)
ys <- predict(lmObj, newdata = data.frame(x = xs))
lines(xs, ys, col = "red", lty = 2, lwd = 2)
cor.test(x, y, method="spearman")
x <- as.vector(percentVar[pairs$phylogroup == 'ecol2nmen',1])
y <- as.vector(percentVar[pairs$phylogroup == 'ecol2nmen',2])
lmObj <- lm( y ~ x )
xs <- range(x)
ys <- predict(lmObj, newdata = data.frame(x = xs))
lines(xs, ys, col = "red", lty = 2, lwd = 2)
cor.test(x, y, method="spearman")
legend("bottomright", c("Shallow", "Deep"), col=c("purple","violet"), fil=c("purple","violet"),cex=0.6)
legend(73,22,c("RNP","Cis-regulatory","tRNA+synthetase","Ribosome","Dual-function"), pch=c(3,4,8,11,14) ,cex=0.3)
dev.off()
