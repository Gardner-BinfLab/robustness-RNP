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
