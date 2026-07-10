#!/usr/bin/Rscript

fluoroRP <- read.table("data/fluoro/figure3b.tsv", sep="\t", header=TRUE)
fluoroRP<-fluoroRP[,c(3,4,1,2)]

pdf(file="manuscript/figures/figure3b.pdf", width=18, height=18)
par(cex=6.0,las=2, mar = c(5,4,0,2) + .1) #c(bottom, left, top, right).
col <- c(c("orchid1", "lightpink", "skyblue", "lightblue1" ))
boxplot(fluoroRP, col=col, ylab="Relative fluorescence", xlab="", main="", outline=FALSE, xaxt = "n", ylim=c(0,2))
axis(1, 1:4, c("Parental","Mutants", "Parental","Mutants")) #
cols<-c("deeppink4","deeppink4","dodgerblue4","dodgerblue4")
for (i in 1:length(fluoroRP[1,])){
    xp <-  rnorm(length(fluoroRP[,i]), mean = 0, sd = 0.1)   
    points(i+xp, fluoroRP[,i],col=cols[i], pch=20, cex=0.5)
}
text(4,1.825,"4.09\n mutations/kb",cex=0.5)
text(2,1.825,"3.94\n mutations/kb",cex=0.5)

legend("bottomleft", c("RNA", "", "Protein", ""), col=col, fil=col,cex=0.5, ncol = 1)
dev.off()

wilcox.test(fluoroRP$RNA.mutants, fluoroRP$protein.mutants)

######################################################################

fluoroR <- read.table("data/fluoro/figure3c-1.tsv", sep="\t", header=TRUE)
#get relative values:
mean.Fluor.R <- mean( apply(fluoroR[,4:6], 1, median)[fluoroR$NoMutations == 0] )
rel.Fluor.R <- apply(fluoroR[,4:6], 1, median)/mean.Fluor.R

fluoroP <- read.table("data/fluoro/figure3c-2.tsv", sep="\t", header=TRUE)
#get relative values:
mean.Fluor.P <- mean( apply(fluoroP[,4:5], 1, median)[fluoroP$NoMutations == 0] )
rel.Fluor.P <- apply(fluoroP[,4:5], 1, median)/mean.Fluor.P

pdf(file="manuscript/figures/figure3c.pdf", width=18, height=18)
par(cex=6.0,las=1, mar = c(5,4,0,2) + .1) #c(bottom, left, top, right).
plot(fluoroR$NoMutations, rel.Fluor.R, col="deeppink4", pch=20, ylim=c(0,1.5),xlim=c(0,6), xlab="Number of mutations", ylab="Relative fluorescence", cex=0.5  )
lines(0:4, predict(loess(rel.Fluor.R ~ fluoroR$NoMutations), 0:4), col = "deeppink4", lwd=8)
points(fluoroP$NoMutations, rel.Fluor.P, col="dodgerblue4", pch=20, cex=0.5)
lines( 0:6, predict(loess(rel.Fluor.P ~ fluoroP$NoMutations), 0:6), col = "dodgerblue4", lwd=8)
legend("topright", c("RNA", "Protein"), col=c("deeppink4", "dodgerblue4"), fil=c("deeppink4", "dodgerblue4"),cex=0.5)
dev.off()

######################################################################
#Considering:
#[mut.Fl - mean(single WT.Fl)] / mean(single WT.Fl)


mean.fluoroR0 <- mean( apply(fluoroR[,4:6], 1, median)[fluoroR$NoMutations == 0] )
rel.Delta.Fluor.R <- ( apply(fluoroR[,4:6], 1, median) - mean.fluoroR0 )/mean.fluoroR0

mean.fluoroP0 <- mean( apply(fluoroP[,4:5], 1, median)[fluoroP$NoMutations == 0] )
rel.Delta.Fluor.P <- ( apply(fluoroP[,4:5], 1, median) - mean.fluoroP0 )/mean.fluoroP0

pdf(file="manuscript/figures/figure3c-delta.pdf", width=18, height=18)
par(cex=6.0,las=1, mar = c(5,4,0,2) + .1) #c(bottom, left, top, right).
plot(fluoroR$NoMutations, rel.Delta.Fluor.R, col="deeppink4", pch=20, ylim=c(-1.0,0.5),xlim=c(0,6), xlab="Number of mutations per molecule", ylab="Change in fluorescence", cex=0.5  )
lines(0:4, predict(loess(rel.Delta.Fluor.R ~ fluoroR$NoMutations), 0:4), col = "deeppink4", lwd=8)
points(fluoroP$NoMutations, rel.Delta.Fluor.P, col="dodgerblue4", pch=20, cex=0.5)
lines( 0:6, predict(loess(rel.Delta.Fluor.P ~ fluoroP$NoMutations), 0:6), col = "dodgerblue4", lwd=8)
legend("topright", c("RNA", "Protein"), col=c("deeppink4", "dodgerblue4"), fil=c("deeppink4", "dodgerblue4"),cex=0.5)
dev.off()


