#!/usr/bin/Rscript

#BACTERIA
 mRNA.bac <- read.table("data/genomes/synonNonsynon-mRNA-results.tsv",  sep="\t", header=TRUE)
ncRNA.bac <- read.table("data/genomes/synonNonsynon-ncRNA-results.tsv", sep="\t", header=TRUE)

 mRNA.bac.shallow <-  mRNA.bac[ mRNA.bac$distID == "shallow",]$totVar /  mRNA.bac[ mRNA.bac$distID == "shallow",]$len
 mRNA.bac.deep    <-  mRNA.bac[ mRNA.bac$distID == "deep",   ]$totVar /  mRNA.bac[ mRNA.bac$distID == "deep",   ]$len
ncRNA.bac.shallow <- ncRNA.bac[ncRNA.bac$distID == "shallow",]$totVar / ncRNA.bac[ncRNA.bac$distID == "shallow",]$len
ncRNA.bac.deep    <- ncRNA.bac[ncRNA.bac$distID == "deep",   ]$totVar / ncRNA.bac[ncRNA.bac$distID == "deep",   ]$len

len.mRNA.bac.shallow  <- length( mRNA.bac.shallow)
len.mRNA.bac.deep     <- length( mRNA.bac.deep)
len.ncRNA.bac.shallow <- length(ncRNA.bac.shallow)
len.ncRNA.bac.deep    <- length(ncRNA.bac.deep)

maxRows.bac <- max(len.mRNA.bac.shallow,len.mRNA.bac.deep,len.ncRNA.bac.shallow,len.ncRNA.bac.deep)
conservationRP.bac <- matrix(NA,ncol=4,nrow=maxRows.bac)
colnames(conservationRP.bac)<-c("Shallow \n(RNA)", "Shallow\n(Protein)", "Deep\n(RNA)", "Deep\n(Protein)")

conservationRP.bac[1:len.ncRNA.bac.shallow,1] <- ncRNA.bac.shallow
conservationRP.bac[1:len.mRNA.bac.shallow, 2] <- mRNA.bac.shallow
conservationRP.bac[1:len.ncRNA.bac.deep,   3] <- ncRNA.bac.deep
conservationRP.bac[1:len.mRNA.bac.deep,    4] <- mRNA.bac.deep

pdf(file="manuscript/figures/figure1a.pdf", width=12, height=12)
par(cex=2.5,las=2)
col <- c("orchid1", "skyblue", "lightpink", "lightblue1" )
boxplot(conservationRP.bac, col=col, ylab="Total nucleotide variation per site", xlab="", main="Nucleotide variation in Bacteria", outline=FALSE, xaxt = "n")
axis(1, 1:4, colnames(conservationRP.bac))
cols<-c("deeppink4","dodgerblue4","deeppink4","dodgerblue4")
for (i in 1:4){
    xp <-  rnorm(length(conservationRP.bac[,i]), mean = 0, sd = 0.1)   
    points(i+xp, conservationRP.bac[,i],col = "white", bg=cols[i], pch=21, cex=0.75)
}
#legend("bottomright", c("RNA", "", "Protein", ""), col=col[c(1,3,2,4)], fil=col[c(1,3,2,4)],cex=0.6, ncol = 2)
dev.off()

#system("convert -flatten conservation-new-old-genes.pdf conservation-new-old-genes.png")

wilcox.test(conservationRP.bac[,1], conservationRP.bac[,2])
wilcox.test(conservationRP.bac[,3], conservationRP.bac[,4])

######################################################################
#FUNGI 
mRNA.fun <- read.table("data/genomes-fungi/synonNonsynon-mRNA-results.tsv", sep="\t", header=TRUE)
ncRNA.fun <- read.table("data/genomes-fungi/synonNonsynon-ncRNA-results.tsv", sep="\t", header=TRUE)

 mRNA.fun.shallow <-  mRNA.fun[ mRNA.fun$distID == "shallow",]$totVar /  mRNA.fun[ mRNA.fun$distID == "shallow",]$len
 mRNA.fun.deep    <-  mRNA.fun[ mRNA.fun$distID == "deep",   ]$totVar /  mRNA.fun[ mRNA.fun$distID == "deep",   ]$len
ncRNA.fun.shallow <- ncRNA.fun[ncRNA.fun$distID == "shallow",]$totVar / ncRNA.fun[ncRNA.fun$distID == "shallow",]$len
ncRNA.fun.deep    <- ncRNA.fun[ncRNA.fun$distID == "deep",   ]$totVar / ncRNA.fun[ncRNA.fun$distID == "deep",   ]$len

len.mRNA.fun.shallow  <- length( mRNA.fun.shallow)
len.mRNA.fun.deep     <- length( mRNA.fun.deep)
len.ncRNA.fun.shallow <- length(ncRNA.fun.shallow)
len.ncRNA.fun.deep    <- length(ncRNA.fun.deep)

maxRows.fun <- max(len.mRNA.fun.shallow,len.mRNA.fun.deep,len.ncRNA.fun.shallow,len.ncRNA.fun.deep)
conservationRP.fun <- matrix(NA,ncol=4,nrow=maxRows.fun)
colnames(conservationRP.fun)<-c("Shallow \n(RNA)", "Shallow\n(Protein)", "Deep\n(RNA)", "Deep\n(Protein)")

conservationRP.fun[1:len.ncRNA.fun.shallow,1] <- ncRNA.fun.shallow
conservationRP.fun[1:len.mRNA.fun.shallow, 2] <- mRNA.fun.shallow
conservationRP.fun[1:len.ncRNA.fun.deep,   3] <- ncRNA.fun.deep
conservationRP.fun[1:len.mRNA.fun.deep,    4] <- mRNA.fun.deep

pdf(file="manuscript/figures/figure1c.pdf", width=12, height=12)
par(cex=2.5,las=2)
col <- c("orchid1", "skyblue", "lightpink", "lightblue1" )
boxplot(conservationRP.fun, col=col, ylab="Total nucleotide variation per site", xlab="", main="Nucleotide variation in Fungi", outline=FALSE, xaxt = "n")
axis(1, 1:4, colnames(conservationRP.fun))
cols<-c("deeppink4","dodgerblue4","deeppink4","dodgerblue4")
for (i in 1:4){
    xp <-  rnorm(length(conservationRP.fun[,i]), mean = 0, sd = 0.1)   
    points(i+xp, conservationRP.fun[,i],col = "white", bg=cols[i], pch=21, cex=0.75)
}
#legend("bottomright", c("RNA", "", "Protein", ""), col=col[c(1,3,2,4)], fil=col[c(1,3,2,4)],cex=0.6, ncol = 2)
dev.off()

wilcox.test(conservationRP.fun[,1], conservationRP.fun[,2])
wilcox.test(conservationRP.fun[,3], conservationRP.fun[,4])

######################################################################
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

synonRP.bac <- matrix(NA,ncol=10,nrow=maxRows.bac)
synonRP.bac[1:len.mRNA.bac.shallow,  1] <- mRNA.bac[mRNA.bac$distID   == "shallow",]$aaSynProp
synonRP.bac[1:len.mRNA.bac.deep,     2] <- mRNA.bac[mRNA.bac$distID   == "deep",]$aaSynProp
synonRP.bac[1:len.ncRNA.bac.shallow, 3] <- ncRNA.bac[ncRNA.bac$distID == "shallow",]$rySynProp
synonRP.bac[1:len.ncRNA.bac.deep,    4] <- ncRNA.bac[ncRNA.bac$distID == "deep",]$rySynProp
synonRP.bac[1:len.mRNA.bac.shallow,  5] <- mRNA.bac[mRNA.bac$distID   == "shallow",]$banpSynProp
synonRP.bac[1:len.mRNA.bac.deep,     6] <- mRNA.bac[mRNA.bac$distID   == "deep",]$banpSynProp
synonRP.bac[1:len.mRNA.bac.shallow,  7] <- mRNA.bac[mRNA.bac$distID   == "shallow",]$blossSynProp
synonRP.bac[1:len.mRNA.bac.deep,     8] <- mRNA.bac[mRNA.bac$distID   == "deep",]$blossSynProp
synonRP.bac[1:len.ncRNA.bac.shallow, 9] <- ncRNA.bac[ncRNA.bac$distID == "shallow",]$bpSynProp
synonRP.bac[1:len.ncRNA.bac.deep,   10] <- ncRNA.bac[ncRNA.bac$distID == "deep",]$bpSynProp

synonRP.fun <- matrix(NA,ncol=10,nrow=maxRows.fun)
synonRP.fun[1:len.mRNA.fun.shallow,  1] <- mRNA.fun[mRNA.fun$distID   == "shallow",]$aaSynProp
synonRP.fun[1:len.mRNA.fun.deep,     2] <- mRNA.fun[mRNA.fun$distID   == "deep",]$aaSynProp
synonRP.fun[1:len.ncRNA.fun.shallow, 3] <- ncRNA.fun[ncRNA.fun$distID == "shallow",]$rySynProp
synonRP.fun[1:len.ncRNA.fun.deep,    4] <- ncRNA.fun[ncRNA.fun$distID == "deep",]$rySynProp
synonRP.fun[1:len.mRNA.fun.shallow,  5] <- mRNA.fun[mRNA.fun$distID   == "shallow",]$banpSynProp
synonRP.fun[1:len.mRNA.fun.deep,     6] <- mRNA.fun[mRNA.fun$distID   == "deep",]$banpSynProp
synonRP.fun[1:len.mRNA.fun.shallow,  7] <- mRNA.fun[mRNA.fun$distID   == "shallow",]$blossSynProp
synonRP.fun[1:len.mRNA.fun.deep,     8] <- mRNA.fun[mRNA.fun$distID   == "deep",]$blossSynProp
synonRP.fun[1:len.ncRNA.fun.shallow, 9] <- ncRNA.fun[ncRNA.fun$distID == "shallow",]$bpSynProp
synonRP.fun[1:len.ncRNA.fun.deep,   10] <- ncRNA.fun[ncRNA.fun$distID == "deep",]$bpSynProp


pdf(file="manuscript/figures/suppfigure1.pdf", width=16, height=20)
par(mfrow=c(2,1)) #c(bottom, left, top, right).  c(5, 4, 4, 2)
par(cex=2.5,las=2, mar = c(7,4,2,2) + .1, tck=-0.01)
cols <- c(c("skyblue", "lightblue1", "orchid1", "lightpink", "skyblue", "lightblue1", "skyblue", "lightblue1", "orchid1", "lightpink" ))
boxplot(synonRP.bac, col=cols, ylab="Proportion of (near) neutral mutations", xlab="", main="Bacteria", outline=FALSE, xaxt = "n")
names <- c("Degeneracy", "Biochemistry\n(R<->Y)", "Biochemistry\n(BANP)", "BLOSUM", "Secondary\nstructure" )
axis(1, 2*(1:length(names))-0.5, names ) #
cols<-c("dodgerblue4","dodgerblue4","deeppink4","deeppink4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","deeppink4","deeppink4")
for (i in 1:length(synonRP.bac[1,])){
    #ii <- order[i] #sortedMedians$ix[i]
    xp <-  rnorm(length(synonRP.bac[,i]), mean = 0, sd = 0.1)   
    points(i+xp, synonRP.bac[,i],col=cols[i], pch=20, cex=0.75)
}
legend("bottomright", c("Shallow RNA", "Deep RNA", "Shallow protein", "Deep protein"), col=c("orchid1", "lightpink", "skyblue", "lightblue1"), fil=c("orchid1", "lightpink", "skyblue", "lightblue1"),cex=0.65)
####
par(cex=2.5,las=2, mar = c(7,4,2,2) + .1, tck=-0.01)
cols <- c(c("skyblue", "lightblue1", "orchid1", "lightpink", "skyblue", "lightblue1", "skyblue", "lightblue1", "orchid1", "lightpink" ))
boxplot(synonRP.fun, col=cols, ylab="Proportion of (near) neutral mutations", xlab="", main="Fungi", outline=FALSE, xaxt = "n")
names <- c("Degeneracy", "Biochemistry\n(R<->Y)", "Biochemistry\n(BANP)", "BLOSUM", "Secondary\nstructure" )
axis(1, 2*(1:length(names))-0.5, names ) #
cols<-c("dodgerblue4","dodgerblue4","deeppink4","deeppink4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","deeppink4","deeppink4")
for (i in 1:length(synonRP.fun[1,])){
    #ii <- order[i] #sortedMedians$ix[i]
    xp <-  rnorm(length(synonRP.fun[,i]), mean = 0, sd = 0.1)   
    points(i+xp, synonRP.fun[,i],col=cols[i], pch=20, cex=0.75)
}
legend("bottomright", c("Shallow RNA", "Deep RNA", "Shallow protein", "Deep protein"), col=c("orchid1", "lightpink", "skyblue", "lightblue1"), fil=c("orchid1", "lightpink", "skyblue", "lightblue1"),cex=0.65)
dev.off()


######################################################################
######################################################################
#Supp figure 2:

rnaVprotDot <- function(pairs, ncRNA, mRNA) {
percentVar           <- matrix(NA,ncol=2,nrow=length(pairs$RNA))
colnames(percentVar) <- c("RNA.var", "protein.var")
for(i in 1:length(pairs$RNA)){
      percentVar[i,1] <- 100*ncRNA[ncRNA$ID == as.character(pairs[i,]$RNA),    ]$totVar / ncRNA[ncRNA$ID == as.character(pairs[i,]$RNA),    ]$len      
      percentVar[i,2] <- 100*mRNA[  mRNA$ID == as.character(pairs[i,]$protein),]$totVar / mRNA[  mRNA$ID == as.character(pairs[i,]$protein),]$len            
}
par(cex=2.5,las=2, mar = c(5,4,0,2) + .1, tck=-0.01) #c(bottom, left, top, right).  c(5, 4, 4, 2)
plot(NA, NA, ylim=c(0,90), xlim=c(0,90), ylab="mRNA variation (%)", xlab="ncRNA variation (%)")
#plot lines first
for(i in 1:(length(pairs$RNA)-1)){
      for(j in (i+1):(length(pairs$RNA)-1)){
	    col <- "purple"; 
      	    if(pairs[i,]$phylogroup == 'deep'){
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
      if(pairs[i,]$phylogroup == 'deep'){
      	 pch <- 22; bg <- "violet"; 
      }
      points(percentVar[i,1], percentVar[i,2], pch=pch, col=col, bg = bg)
      pch2 <- switch(as.character(pairs[i,]$type),
        "RNP"             = 3,
     	"Cis-regulatory"  = 4,
      	"tRNA+synthetase" = 8, 
      	"Ribosome"        = 11,
      	"Dual-function"   = 13,
        "Spliceosome"     = 14
      )
      points(percentVar[i,1], percentVar[i,2], pch=pch2, col="black", cex = 0.5)
}
## d <- 2
## rect(
## min(percentVar[pairs$phylogroup == 'shallow',1]) - d,
## min(percentVar[pairs$phylogroup == 'shallow',2]) - d,
## max(percentVar[pairs$phylogroup == 'shallow',1]) + d,
## max(percentVar[pairs$phylogroup == 'shallow',2]) + d
## )
## rect(
## min(percentVar[pairs$phylogroup == 'deep',1]) - d,
## min(percentVar[pairs$phylogroup == 'deep',2]) - d,
## max(percentVar[pairs$phylogroup == 'deep',1]) + d,
## max(percentVar[pairs$phylogroup == 'deep',2]) + d
## )
x <- as.vector(percentVar[pairs$phylogroup == 'shallow',1])
y <- as.vector(percentVar[pairs$phylogroup == 'shallow',2])
lmObj <- lm( y ~ x )
xs <- range(x)
ys <- predict(lmObj, newdata = data.frame(x = xs))
lines(xs, ys, col = "purple4", lty = 2, lwd = 4)
    s.cor <- cor.test(x, y, method="spearman")
    print("SHALLOW")
    print(s.cor)
x <- as.vector(percentVar[pairs$phylogroup == 'deep',1])
y <- as.vector(percentVar[pairs$phylogroup == 'deep',2])
lmObj <- lm( y ~ x )
xs <- range(x)
ys <- predict(lmObj, newdata = data.frame(x = xs))
lines(xs, ys, col = "hotpink4", lty = 2, lwd = 4)
    s.cor <- cor.test(x, y, method="spearman")
    print("DEEP")
    print(s.cor)
legend("bottomright", c("Shallow", "Deep"), col=c("purple","violet"), fil=c("purple","violet"),cex=0.6)
legend(73,22,c("RNP","Cis-regulatory","tRNA+synthetase","Ribosome","Dual-function", "Spliceosome"), pch=c(3,4,8,11,13,14) ,cex=0.3)
}


pairs   <- read.table("data/genomes/rna-protein-pairs.tsv",  sep="\t", header=TRUE)
pdf(file="manuscript/figures/suppfigure2a.pdf", width=11, height=10)
rnaVprotDot(pairs, ncRNA.bac, mRNA.bac)
dev.off()


pairs   <- read.table("data/genomes-fungi/rna-protein-pairs.tsv",  sep="\t", header=TRUE)
pdf(file="manuscript/figures/suppfigure2b.pdf", width=11, height=10)
rnaVprotDot(pairs, ncRNA.fun, mRNA.fun)
dev.off()

#convert  -flatten -density 300 -trim  manuscript/figures/suppfigure2a.pdf  -quality 100   manuscript/figures/suppfigure2a.png
#convert  -flatten -density 300 -trim  manuscript/figures/suppfigure2b.pdf  -quality 100   manuscript/figures/suppfigure2b.png




percentVar           <- matrix(NA,ncol=2,nrow=length(pairs$RNA))
colnames(percentVar) <- c("RNA.var", "protein.var")
for(i in 1:length(pairs$RNA)){
      percentVar[i,1] <- 100*ncRNA.bac[ncRNA.bac$ID == as.character(pairs[i,]$RNA),    ]$totVar / ncRNA.bac[ncRNA.bac$ID == as.character(pairs[i,]$RNA),    ]$len      
      percentVar[i,2] <- 100*mRNA.bac[  mRNA.bac$ID == as.character(pairs[i,]$protein),]$totVar / mRNA.bac[  mRNA.bac$ID == as.character(pairs[i,]$protein),]$len            
}






















