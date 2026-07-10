#!/usr/bin/Rscript

require(graphics)
library(hash)
library(RColorBrewer)
#library(gplots)
library(Manu)

my.cols.corr <- c("#000000", brewer.pal(11, "PiYG"), brewer.pal(9, "Greys") )
mutBroc <- read.table("sam2var.broccoli-118-318.txt", sep="\t", header=TRUE)
mutGFP <- read.table("sam2var.egfp-2410-2610.txt",    sep="\t", header=TRUE)
#Re-normalise counts for each bin to expectation from cell-sorting
#Barcodes groups were balance, this to be corrected for.
#                    low    lo.med med.hi  high              

samStatsBroc  <- read.table("broccoli.sam.stats", sep="\t", header=TRUE)
samStatsGFP   <- read.table("egfp.sam.stats", sep="\t", header=TRUE)

#ADJUST APPROACH TO TPMSs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###########
#Function to convert mutation tables from sam2var to a matrix for stacked barplots: 
#--altered to return TPMs, with a pseudocount added. 
mutTable2stackMatrixTPM <- function( mutTable, column ){
     x<-matrix(100/10^6, nrow = 201, ncol = 4, dimnames = list(1:201, c('A','C','G','T'))) ##initialise with a pseudocount
     totReads <- sum(mutTable[, column])
     x[mutTable[mutTable$querynuc == 'A', 1], 1] <- x[mutTable[mutTable$querynuc == 'A', 1], 1] + (10^6)*(mutTable[mutTable$querynuc == 'A', column])/totReads
     x[mutTable[mutTable$querynuc == 'C', 1], 2] <- x[mutTable[mutTable$querynuc == 'C', 1], 2] + (10^6)*(mutTable[mutTable$querynuc == 'C', column])/totReads
     x[mutTable[mutTable$querynuc == 'G', 1], 3] <- x[mutTable[mutTable$querynuc == 'G', 1], 3] + (10^6)*(mutTable[mutTable$querynuc == 'G', column])/totReads
     x[mutTable[mutTable$querynuc == 'T', 1], 4] <- x[mutTable[mutTable$querynuc == 'T', 1], 4] + (10^6)*(mutTable[mutTable$querynuc == 'T', column])/totReads	     
     return(x)
}
###########

paired.cols <- palette.colors(palette='Paired')
cols = c("lightgray", paired.cols[c(5,1,9,7)])    #c("lightgray","pink", "lightblue","orchid","orange") #neutral, GoF, LoF, Uniform, Other
brocCol <- paired.cols[3] #'seagreen1'
gfpCol  <- paired.cols[4] #'darkgreen'


######################################################################

pdf(file="plotMutagenesis.pdf", width=64, height=64)


########
#stats on number of variants

par(mfrow=c(2,2),las=2,cex=7.0)
barplot(t(as.matrix(samStatsBroc[,2])),xlim=c(0,20),col=brocCol, ylab="#reads", xlab="#variants", main="Broccoli number of variants per read", names.arg = samStatsBroc[,1])
barplot(t(as.matrix(samStatsBroc[,3])),xlim=c(0,20),col=brocCol, ylab="%reads", xlab="#variants", main="Broccoli proportion of variants per read", names.arg = samStatsBroc[,1])

barplot(t(as.matrix(samStatsGFP[,2])),xlim=c(0,20),col=gfpCol, ylab="#reads", xlab="#variants", main="GFP number of variants per read", names.arg = samStatsGFP[,1])
barplot(t(as.matrix(samStatsGFP[,3])),xlim=c(0,20),col=gfpCol, ylab="%reads", xlab="#variants", main="GFP proportion of variants per read", names.arg = samStatsGFP[,1])


#Presort population distribution:
presort.tpm.broc <- mutTable2stackMatrixTPM(mutBroc, 4)
presort.tpm.gfp  <- mutTable2stackMatrixTPM(mutGFP,  4)
broc.sum.sorted  <- mutTable2stackMatrixTPM(mutBroc, 5)+mutTable2stackMatrixTPM(mutBroc, 6)+mutTable2stackMatrixTPM(mutBroc, 7)+mutTable2stackMatrixTPM(mutBroc, 8)
gfp.sum.sorted   <- mutTable2stackMatrixTPM(mutGFP,  5)+mutTable2stackMatrixTPM(mutGFP,  6)+mutTable2stackMatrixTPM(mutGFP,  7)+mutTable2stackMatrixTPM(mutGFP,  8)

mn <- 0
mx <- 6
brks <- seq(mn,mx,length.out = 50)
broc.hist         <- hist(log10(presort.tpm.broc +1), breaks=brks, plot=F)
gfp.hist          <- hist(log10(presort.tpm.gfp  +1), breaks=brks, plot=F)
broc.mn.sort.hist <- hist(log10(broc.sum.sorted/4+1), breaks=brks, plot=F)
 gfp.mn.sort.hist <- hist(log10(gfp.sum.sorted/4 +1), breaks=brks, plot=F)

par(mfrow=c(1,1),las=2,cex=15.0)
plot(NA,NA, xlim=c(mn,mx), ylim=c(0,0.6), main="Presort distributions",xlab="log10(TPM per variant)", ylab="Normalised Freq.")
lines(        broc.hist$mids,         broc.hist$counts/max(        broc.hist$counts), lwd=18, col=brocCol)#"darkslateblue"        )
lines(broc.mn.sort.hist$mids, broc.mn.sort.hist$counts/max(broc.mn.sort.hist$counts), lwd=18, col=brocCol, lty=3)#"darkslateblue" )
lines(          gfp.hist$mids,         gfp.hist$counts/max(         gfp.hist$counts), lwd=18, col=gfpCol)#"green"                )
lines( gfp.mn.sort.hist$mids,  gfp.mn.sort.hist$counts/max( gfp.mn.sort.hist$counts), lwd=18, col=gfpCol, lty=3)#"green",         lty=3 )

legend("topright", c("Brocolli", "GFP"), fil=c(brocCol, gfpCol)) 
legend(4,0.425, c("Presort", "Sorted"), col=c("black", "black"), lty=c(1,3), lwd=c(18,18),cex=0.75)

###################################
#Renormalised TPMs:
#Summary stats of presort data (remove the zeros -- these are not useful for the analysis):
x<-log10(c(presort.tpm.broc, presort.tpm.gfp)+1)
x<-x[x>1]
mn.presort <- mean(x)
sd.presort <-   sd(x)

brks <- seq(-10,10,length.out = 50)
broc.z.hist         <- hist( (log10(presort.tpm.broc +1)-mn.presort)/sd.presort, breaks=brks, plot=F)
broc.z.mn.sort.hist <- hist( (log10(broc.sum.sorted/4+1)-mn.presort)/sd.presort, breaks=brks, plot=F)
 gfp.z.hist         <- hist( (log10(presort.tpm.gfp  +1)-mn.presort)/sd.presort, breaks=brks, plot=F)
 gfp.z.mn.sort.hist <- hist( (log10(gfp.sum.sorted/4 +1)-mn.presort)/sd.presort, breaks=brks, plot=F)

par(mfrow=c(1,1),las=2,cex=15.0)
plot(NA,NA, xlim=c(-5,5), ylim=c(0,0.8), main="Pre. & post sort distributions",xlab="Z-score (normalised TPM)", ylab="Normalised Freq.")
lines(      broc.z.hist$mids,           broc.z.hist$counts/max(        broc.z.hist$counts), lwd=18, col=brocCol        )
lines(broc.z.mn.sort.hist$mids, broc.z.mn.sort.hist$counts/max(broc.z.mn.sort.hist$counts), lwd=18, col=brocCol, lty=3 )
lines(         gfp.z.hist$mids,          gfp.z.hist$counts/max(         gfp.z.hist$counts), lwd=18, col=gfpCol         )
lines( gfp.z.mn.sort.hist$mids,  gfp.z.mn.sort.hist$counts/max( gfp.z.mn.sort.hist$counts), lwd=18, col=gfpCol,  lty=3 )

legend("topright", c("Brocolli", "GFP"), fil=c(brocCol, gfpCol)) 
legend(1.5,0.5, c("Presort", "Sorted"), col=c("black", "black"), lty=c(1,3), lwd=c(18,18),cex=0.75)

######################################################################
par(mfrow=c(6,1),las=2,cex=3.0)
#Rework data into a dataframe for stacked barplots:
stackMatrix.Broc.total <- mutTable2stackMatrixTPM(mutBroc, 9)
barplot(t(stackMatrix.Broc.total  ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Total: Broccoli")
legend("topleft", c('A', 'C', 'G', 'T'), fill=c('green', 'blue', 'black', 'red'),cex=1,ncol = 2)
stackMatrix.Broc.high <- mutTable2stackMatrixTPM(mutBroc, 8)
barplot(t(stackMatrix.Broc.high   ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="High")
stackMatrix.Broc.medHigh <- mutTable2stackMatrixTPM(mutBroc, 7)
barplot(t(stackMatrix.Broc.medHigh),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Median-High")
stackMatrix.Broc.lowMed <- mutTable2stackMatrixTPM(mutBroc, 6)
barplot(t(stackMatrix.Broc.lowMed ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Low-Median")
stackMatrix.Broc.low <- mutTable2stackMatrixTPM(mutBroc, 5)
barplot(t(stackMatrix.Broc.low    ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Low")
stackMatrix.Broc.presort <- mutTable2stackMatrixTPM(mutBroc, 4)
barplot(t(stackMatrix.Broc.presort),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="Position", main="Presort")
#dev.off()

pseudocount <- 50

log.Broc.total.presort   <- log10(t(stackMatrix.Broc.total + pseudocount  )/t(stackMatrix.Broc.presort + pseudocount))
log.Broc.high.presort    <- log10(t(stackMatrix.Broc.high + pseudocount   )/t(stackMatrix.Broc.presort + pseudocount))
log.Broc.medHigh.presort <- log10(t(stackMatrix.Broc.medHigh + pseudocount)/t(stackMatrix.Broc.presort + pseudocount))
log.Broc.lowMed.presort  <- log10(t(stackMatrix.Broc.lowMed + pseudocount )/t(stackMatrix.Broc.presort + pseudocount))
log.Broc.low.presort     <- log10(t(stackMatrix.Broc.low + pseudocount    )/t(stackMatrix.Broc.presort + pseudocount))

par(mfrow=c(5,1),las=2,cex=3.0)
barplot(log.Broc.total.presort,  ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.tot/tpm.presort)", xlab="", beside=TRUE, main="Total: Broccoli")
legend("topleft", c('A', 'C', 'G', 'T'), fill=c('green', 'blue', 'black', 'red'),cex=1,ncol = 2)
barplot(log.Broc.high.presort   ,ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.hi /tpm.presort)", xlab="", beside=TRUE, main="High"           )
barplot(log.Broc.medHigh.presort,ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.mhi/tpm.presort)", xlab="", beside=TRUE, main="Median-High"    )
barplot(log.Broc.lowMed.presort, ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.lom/tpm.presort)", xlab="", beside=TRUE, main="Low-Median"     )
barplot(log.Broc.low.presort,    ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.lo /tpm.presort)", xlab="", beside=TRUE, main="Low"            )

##############################

par(mfrow=c(6,1),las=2,cex=3.0)
stackMatrix.GFP.total <- mutTable2stackMatrixTPM(mutGFP, 9)
barplot(t(stackMatrix.GFP.total  ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Total: GFP")
legend("topright", c('A', 'C', 'G', 'T'), fill=c('green', 'blue', 'black', 'red'),cex=1,ncol = 2)
stackMatrix.GFP.high <- mutTable2stackMatrixTPM(mutGFP, 8)
barplot(t(stackMatrix.GFP.high   ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="High")
stackMatrix.GFP.medHigh <- mutTable2stackMatrixTPM(mutGFP, 7)
barplot(t(stackMatrix.GFP.medHigh),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Median-High")
stackMatrix.GFP.lowMed <- mutTable2stackMatrixTPM(mutGFP, 6)
barplot(t(stackMatrix.GFP.lowMed ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Low-Median")
stackMatrix.GFP.low <- mutTable2stackMatrixTPM(mutGFP, 5)
barplot(t(stackMatrix.GFP.low    ),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Low")
stackMatrix.GFP.presort <- mutTable2stackMatrixTPM(mutGFP, 4)
barplot(t(stackMatrix.GFP.presort),ylim=c(0,40000),col=c('green', 'blue', 'black', 'red'), ylab="TPM", xlab="", main="Presort")
#dev.off()

log.GFP.total.presort   <- log10(t(stackMatrix.GFP.total + pseudocount  )/t(stackMatrix.GFP.presort + pseudocount))
log.GFP.high.presort    <- log10(t(stackMatrix.GFP.high + pseudocount   )/t(stackMatrix.GFP.presort + pseudocount))
log.GFP.medHigh.presort <- log10(t(stackMatrix.GFP.medHigh + pseudocount)/t(stackMatrix.GFP.presort + pseudocount))
log.GFP.lowMed.presort  <- log10(t(stackMatrix.GFP.lowMed + pseudocount )/t(stackMatrix.GFP.presort + pseudocount))
log.GFP.low.presort     <- log10(t(stackMatrix.GFP.low + pseudocount    )/t(stackMatrix.GFP.presort + pseudocount))

par(mfrow=c(5,1),las=2,cex=3.0)
barplot(log.GFP.total.presort,  ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio((tpm.tot+p)/(tpm.presort+p))", xlab="", beside=TRUE, main="Total: GFP")
legend("topleft", c('A', 'C', 'G', 'T'), fill=c('green', 'blue', 'black', 'red'),cex=1,ncol = 2)
barplot(log.GFP.high.presort   ,ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.hi /tpm.presort)", xlab="", beside=TRUE, main="High"           )
barplot(log.GFP.medHigh.presort,ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.mhi/tpm.presort)", xlab="", beside=TRUE, main="Median-High"    )
barplot(log.GFP.lowMed.presort, ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.lom/tpm.presort)", xlab="", beside=TRUE, main="Low-Median"     )
barplot(log.GFP.low.presort,    ylim=c(-1,1),col=c('green', 'blue', 'black', 'red'), ylab="log-ratio(tpm.lo /tpm.presort)", xlab="", beside=TRUE, main="Low"            )

#####################

#From Leak's cell sorting experiments:
par(mfrow=c(2,1),cex=5.0)
expectedFreqsBroc <- c(23.06, 50.66, 23.77,  2.50, 100)/100
countsBroc <- matrix(NA, nrow = 3, ncol = 6, dimnames = list(c('WT','Variants','Expected'), c('Presort','Low','Med.-Low','High-Med.','High','Total')))
countsBroc[1,6] <- 106821 # WT total from broccoli.sam.stats
countsBroc[2,1] <- sum(mutBroc$presort )
countsBroc[2,2] <- sum(mutBroc$low     )
countsBroc[2,3] <- sum(mutBroc$low.med )
countsBroc[2,4] <- sum(mutBroc$med.high)
countsBroc[2,5] <- sum(mutBroc$high    )
countsBroc[2,6] <- sum(mutBroc$total   )
countsBroc[3,2:6]  <- expectedFreqsBroc * countsBroc[2,6]
barplot(countsBroc, col=c('black', 'blue', 'gold'), ylab="# reads", xlab="", beside=TRUE, main="Read counts for broccoli: expected (from gates) vs observed")
legend('topleft', c('WT','Variants','Expected'), col=c('black', 'blue', 'gold'), fil=c('black', 'blue', 'gold'))
lines(c(5.5,20.5), c(countsBroc[2,6], countsBroc[2,6])/4, lty=2, lwd=2)
text(5.5, countsBroc[2,6]/4, "Uniform freq.",pos=2)

expectedFreqsGFP  <- c(78.25, 15.16,  6.38,  0.20, 100)/100
countsGFP  <- matrix(NA, nrow = 3, ncol = 6, dimnames = list(c('WT','Variants','Expected'), c('Presort','Low','Med.-Low','High-Med.','High','Total')))
countsGFP[1,6] <- 48163 # WT total from GFP.sam.stats
countsGFP[2,1] <- sum(mutGFP$presort )
countsGFP[2,2] <- sum(mutGFP$low     )
countsGFP[2,3] <- sum(mutGFP$low.med )
countsGFP[2,4] <- sum(mutGFP$med.high)
countsGFP[2,5] <- sum(mutGFP$high    )
countsGFP[2,6] <- sum(mutGFP$total   )
countsGFP[3,2:6]  <- expectedFreqsGFP * countsGFP[2,6]
barplot(countsGFP, col=c('black', 'blue', 'gold'), ylab="# reads", xlab="", beside=TRUE, main="Read counts for GFP: expected (from gates) vs observed")
legend('topleft', c('WT','Variants','Expected'), col=c('black', 'blue', 'gold'), fil=c('black', 'blue', 'gold'))
lines(c(5.5,20.5), c(countsGFP[2,6], countsGFP[2,6])/4, lty=2, lwd=2)
text(5.5, countsGFP[2,6]/4, "Uniform freq.",pos=2)

######################################################################

par(mfrow=c(1,1),cex=7.0)
plot(NA,NA,xlim=c(0.95,4.7), ylim=c(-1.2,1.2), main="Candidate GoF variants", ylab="log(TPM/TPM.presort)", xlab="", xaxt = "n")
axis(1, c(4,3,2,1), c('High','High-Mid', 'Mid-Low','Low'))

thresh <- 0.095
nucs<-c("A","C","G","T")
x <- 1:4
for (i in 1:4  ){
    for (j in 1:201){
    	
    	y <- c(log.GFP.low.presort[i,j], log.GFP.lowMed.presort[i,j], log.GFP.medHigh.presort[i,j], log.GFP.high.presort[i,j])
	if(sum(abs(y)) > 0){
	gfp.spear <- cor.test(x, y, method = "spearman")
		  if(is.na( gfp.spear$p.value ) == FALSE && gfp.spear$p.value < thresh && gfp.spear$estimate>0){#Significant
		  		       lines(x,y, col=gfpCol, lwd=10)
				       if(gfp.spear$estimate>0 && y[4]>0.45){
					     ref<-mutGFP[mutGFP$pos == j, 2]
					     ref <- ref[1]					    
					     text(4,y[4],paste(ref, j, nucs[i],',N~', floor(stackMatrix.GFP.total[j,i]), sep=""),pos=4,col=gfpCol)
				 	}
		  }
	}
	
    	y <- c(log.Broc.low.presort[i,j], log.Broc.lowMed.presort[i,j], log.Broc.medHigh.presort[i,j], log.Broc.high.presort[i,j])
	if(sum(abs(y)) > 0){
		 broc.spear <- cor.test(x, y, method = "spearman")
		 if(is.na( broc.spear$p.value ) == FALSE && broc.spear$p.value < thresh  && broc.spear$estimate>0){#Significant
		       			     lines(x,y, col=brocCol, lwd=10)
					     if(broc.spear$estimate>0 && y[4]>0.45){
						ref<-mutBroc[mutBroc$pos == j, 2]
					     	ref <- ref[1]					    
					     	text(4,y[4],paste(ref, j, nucs[i],',N~', floor(stackMatrix.Broc.total[j,i]), sep=""),pos=4,col=brocCol)
					     }
		       }
		 }
	}    
}
legend(4.05,-0.8, c("Brocolli", "GFP"), fill=c(brocCol, gfpCol))



######################################################################
#SLICES:

#target:119-319
mutBrocSlice <- mutBroc[1 <= mutBroc[,1] & mutBroc[,1] <= 201, ]
nucs <- c('A', 'C', 'G', 'T')
nuc2num<-hash( keys=nucs, values=1:4 )
##
mutBrocSliceCorr <- matrix(NA,nrow=5, ncol=319-119+1)
rownames(mutBrocSliceCorr) <- c(nucs, 'sampleSize')
mutBrocSliceCorr[5,]<-0
##

#target:2410-2610
mutGFPSlice <- mutGFP; #mutGFP[2410 <= mutGFP[,1] & mutGFP[,1] <= 2610, ]


######################################################################
#PRESORT vs TOTALS
par(mfrow=c(2,2),cex=7)
a <- mutBrocSlice$presort / sum(mutBrocSlice$presort)
b <- mutBrocSlice$total   / sum(mutBrocSlice$total)
h <- mutBrocSlice$high    / sum(mutBrocSlice$high)
l <- mutBrocSlice$low     / sum(mutBrocSlice$low )
plot( a, b, xlab='Pre-sort Freq.', ylab='Post-sort Freq.', main='Broccoli'  )
broc.cor.pre.post<-cor.test( a, b   )
text(0,0.015, paste('cor=',     signif(broc.cor.pre.post$estimate, 2) , "\n",
                    'P.value=', signif(broc.cor.pre.post$p.value,  2)),
		    pos=4 )
abline(lm( b ~ a ), col='black', lwd=2)
points(a, h, pch=4, col='red')
abline(lm( h ~ a ), col='red', lwd=2)
points(a, l, pch=3, col='blue')
abline(lm( l ~ a ), col='blue', lwd=2)
legend('bottomright', c('total', 'high', 'low'), pch=c(1,4,3),col=c('black','red','blue'))
###
hist(log10((b+0.001)/(a+0.001)), breaks=seq(-0.65,0.65,length.out = 50), col='lightcyan', xlab='log10((post+d)/(pre+d))', main='Broccoli', xlim=c(-0.6,0.6))
###
a <- mutGFPSlice$presort / sum(mutGFPSlice$presort)
b <- mutGFPSlice$total   / sum(mutGFPSlice$total)
h <- mutGFPSlice$high    / sum(mutGFPSlice$high)
l <- mutGFPSlice$low     / sum(mutGFPSlice$low )
plot( a, b, xlab='Pre-sort Freq.', ylab='Post-sort Freq.', main='GFP'  )
gfp.cor.pre.post<-cor.test( a, b   )
text(0,0.020, paste('cor=',     signif(gfp.cor.pre.post$estimate, 2) , "\n",
                    'P.value=', signif(gfp.cor.pre.post$p.value,  2)),
		    pos=4 )
abline(lm( b ~ a ), col='red', lwd=2)
points(a, h, pch=4, col='red')
abline(lm( h ~ a ), col='red', lwd=2)
points(a, l, pch=3, col='blue')
abline(lm( l ~ a ), col='blue', lwd=2)
###
hist(log10((b+0.001)/(a+0.001)), breaks=seq(-0.65,0.65,length.out = 50), col='lightcyan', xlab='log10((post+d)/(pre+d))', main='GFP', xlim=c(-0.6,0.6))

######################################################################
#

par(mfrow=c(2,1),cex=7.0)
plot( c(log.Broc.high.presort    + log.Broc.medHigh.presort - log.Broc.lowMed.presort - log.Broc.low.presort),
      c(log.Broc.medHigh.presort + log.Broc.lowMed.presort  - log.Broc.high.presort   - log.Broc.low.presort),
      type='p',pch=20,col='red',cex=2, xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),xlab='high + med.high - low.med - low', ylab='med.high + low.med - high - low',main='Broccoli: log ratios of counts with presort')
plot( c(log.GFP.high.presort    + log.GFP.medHigh.presort - log.GFP.lowMed.presort - log.GFP.low.presort),
      c(log.GFP.medHigh.presort + log.GFP.lowMed.presort  - log.GFP.high.presort   - log.GFP.low.presort),
      type='p',pch=20,col='red',cex=2, xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),xlab='high + med.high - low.med - low', ylab='med.high + low.med - high - low',main='GFP: log ratios of counts with presort')

######################
#HISTOGRAMS
par(mfrow=c(2,5),las=3,cex=5.0)
numBreaks <- 50;
breaks <- c(-1001, seq(-10, 10, length.out=numBreaks), 1001)
hist(c(log.Broc.low.presort    ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=brocCol, breaks=breaks, main='Low')
hist(c(log.Broc.lowMed.presort ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=brocCol, breaks=breaks, main='Low-Med.')
hist(c(log.Broc.medHigh.presort), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=brocCol, breaks=breaks, main='Med.-High')
hist(c(log.Broc.high.presort   ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=brocCol, breaks=breaks, main='High')
hist(c(log.Broc.total.presort  ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=brocCol, breaks=breaks, main='Total')
mtext('Broccoli', side = 4, line = 1)

hist(c(log.GFP.low.presort    ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=gfpCol, breaks=breaks, main='Low')
hist(c(log.GFP.lowMed.presort ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=gfpCol, breaks=breaks, main='Low-Med.')
hist(c(log.GFP.medHigh.presort), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=gfpCol, breaks=breaks, main='Med.-High')
hist(c(log.GFP.high.presort   ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=gfpCol, breaks=breaks, main='High')
hist(c(log.GFP.total.presort  ), xlim=c(-4,4), xlab='log(tpm/tpm.presort)', col=gfpCol, breaks=breaks, main='Total')
mtext('GFP', side = 4, line = 1)



######################################################################


ncols <- 319-119+1
mutBrocSliceCorr <- matrix(NA,nrow=5, ncol=ncols, dimnames=list(c(nucs, 'sampleSize'), 1:ncols))
#rownames(mutBrocSliceCorr) <- c(nucs, 'sampleSize')
mutBrocSliceCorr[5,]<-0
##
minSigGoFBroc <-  1
maxSigLoFBroc <- -1

 for (i in 1:length(mutBrocSlice$pos)){
     pos <- mutBrocSlice$pos[i]  #- 119 + 1
     colnames(mutBrocSliceCorr)[pos] <- mutBrocSlice$pos[i]
     mutBrocSliceCorr[ nuc2num[[ mutBrocSlice$refnuc[i]  ]] , pos ]  <- -2 #Set reference nucleotides to "-2" 
     if(mutBrocSlice$total[i]>10){#If there are a reasonable number of reads, compute a spearman correlation 
 	spear <- cor.test(c(-2,-1,1,2), as.numeric(mutBrocSlice[i,5:8]), method = "spearman")
 	mutBrocSliceCorr[ nuc2num[[ mutBrocSlice$querynuc[i]  ]] , pos ]  <- as.numeric(spear$estimate)
 	}
	
	if(spear$p.value < 0.1 && as.numeric(spear$estimate)>0){#Significant & GoF 
		minSigGoFBroc <- min(minSigGoFBroc, as.numeric(spear$estimate))
	}
	else if (spear$p.value < 0.1 && as.numeric(spear$estimate)>0){#Significant & LoF 
		maxSigLoFBroc <- max(maxSigLoFBroc, as.numeric(spear$estimate))
	}
	
	
 	mutBrocSliceCorr[5,pos] <- mutBrocSliceCorr[5,mutBrocSlice$pos[i] ] + mutBrocSlice$total[i]
 }
#mutBrocSliceCorr

#target:2410-2610
ncols <- 2610-2410+1
mutGFPSliceCorr <- matrix(NA,nrow=5, ncol=ncols, dimnames=list(c(nucs, 'sampleSize'), 1:ncols))
#rownames(mutGFPSliceCorr) <- c(nucs, 'sampleSize')
mutGFPSliceCorr[5,]<-0
##

minSigGoFGFP <-  1
maxSigLoFGFP <- -1
for (i in 1:length(mutGFPSlice$pos)){
     pos <- mutGFPSlice$pos[i]  #- 2410 + 1
     colnames(mutGFPSliceCorr)[pos] <- mutGFPSlice$pos[i]
     mutGFPSliceCorr[ nuc2num[[ mutGFPSlice$refnuc[i]  ]] , pos ]  <- -2 #Set reference nucleotides to "-2" 
     if(mutGFPSlice$total[i]>10){#If there are a reasonable number of reads, compute a spearman correlation 
 	spear <- cor.test(c(-2,-1,1,2), as.numeric(mutGFPSlice[i,5:8]), method = "spearman")
 	mutGFPSliceCorr[ nuc2num[[ mutGFPSlice$querynuc[i]  ]] , pos ]  <- as.numeric(spear$estimate)
 	}	

	print(spear$p.value)
	if(spear$p.value < 0.1 && as.numeric(spear$estimate)>0){#Significant & GoF 
		minSigGoFGFP <- min(minSigGoFGFP, as.numeric(spear$estimate))
		}
	else if (spear$p.value < 0.1 && as.numeric(spear$estimate)>0){#Significant & LoF 
		maxSigLoFGFP <- max(maxSigLoFGFP, as.numeric(spear$estimate))
	}
	
 	mutGFPSliceCorr[5,pos] <- mutGFPSliceCorr[5,mutGFPSlice$pos[i] ] + mutGFPSlice$total[i]
 }

mutGFPSliceCorr[5,] <- log10( mutGFPSliceCorr[5,] +1 )
#mutGFPSliceCorr

#pdf(file="hists-spearman-gfp-broc.pdf", width=18, height=8)
par(mfrow=c(1,2),las=2, cex=7)
hist(mutBrocSliceCorr[1:4,], breaks=c(-2.1,-1.1,seq(-1,1,by=0.1)), xlim=c(-1,1),xlab='Spearman Correlation Coefficient', main='Broccoli',freq=TRUE, col=brocCol)
#minSigGoFBroc,minSigGoFBroc
lines(c(0.50,0.50),c(0,100), col='red',lwd=3)
text(0.5,70,'GoF?',pos=4)
lines(c(-0.50,-0.50),c(0,100), col='red',lwd=3)
text(-0.5,70,'LoF?',pos=2)
#lines(c(maxSigLoFBroc,maxSigLoFBroc),c(0,100), col='red',lwd=3)
hist( mutGFPSliceCorr[1:4,], breaks=c(-2.1,-1.1,seq(-1,1,by=0.1)), xlim=c(-1,1),xlab='Spearman Correlation Coefficient', main='GFP',freq=TRUE, col=gfpCol)
#minSigGoFGFP,minSigGoFGFP
lines(c(0.50,0.50),c(0,100), col='red',lwd=3)
text(0.5,70,'GoF?',pos=4)
lines(c(-0.50,-0.50),c(0,100), col='red',lwd=3)
text(-0.5,70,'LoF?',pos=2)
#lines(c(maxSigLoFGFP,maxSigLoFGFP),c(0,100), col='red',lwd=3)
#dev.off()

######################################################################
library('gplots')

#Broccoli 
start<-1   #118
end  <-201 #318
numRows <- length(mutBroc$pos)

pooledCountMatrix <- matrix(1,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('presort', 'low', 'low.med', 'med.high', 'high') ) )
correlationMatrix <- matrix(-1.1,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('','','','','Rho') ) )

for (i in 1:numRows){
    if(start <= mutBroc$pos[i] & mutBroc$pos[i] <= end){
    	  j <- mutBroc$pos[i] - start +1
    	  pooledCountMatrix[j, ] <- as.matrix(mutBroc[i,4:8]) + t(as.matrix(pooledCountMatrix[j, 1:5]))
	  
    }
}
#head(as.matrix(pooledCountMatrix))

for (i in 1:length(pooledCountMatrix[,2])){
    correlationMatrix[i,5] <- cor(1:4,  pooledCountMatrix[i,2:5], method='spearman')
}

totalsCounts <- apply(pooledCountMatrix, 1, sum)
norm.pooledCountMatrix <-  pooledCountMatrix/totalsCounts

logTots <- matrix(0.8,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('','','','','total')) )
logTots[,5] <- as.matrix(log10(totalsCounts))


#pdf(file=    "heatmap-readCounts.pdf", width = 5,  height = 15)
par(cex=5)
heatmap.2(norm.pooledCountMatrix,
	col=rev(brewer.pal(n = 11, name = "RdBu")), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=seq(0,0.54,length=11+1), scale="none",
	cexRow=7.0, cexCol=7.0, margins = c(30, 2), key.title = "Fraction reads", Colv=FALSE, Rowv=FALSE,main="Broccoli", notecex=5)
#dev.off()

#pdf(file=    "heatmap-logTotals.pdf", width = 5,  height = 15)
heatmap.2(logTots,
	col=c('white',rev(brewer.pal(n = 11, name = "Spectral"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(0.8,seq(1,4,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, margins = c(30, 2), key.title = "log10(Tot. Reads)", Colv=FALSE, Rowv=FALSE,main="Broccoli")
#dev.off()

#pdf(file=    "heatmap-spearmanRho.pdf", width = 5,  height = 15)
heatmap.2(correlationMatrix,
	col=c('white',rev(brewer.pal(n = 11, name = "PRGn"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(-1.1,seq(-1,1,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, margins = c(30, 2), key.title = "Rho", Colv=FALSE, Rowv=FALSE,main="Broccoli")
#dev.off()

###################################
#GFP

start<-1   #2410
end  <-201 #2610
numRows <- length(mutGFP$pos)

pooledCountMatrix <- matrix(1,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('presort', 'low', 'low.med', 'med.high', 'high') ) )
correlationMatrix <- matrix(-1.1,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('','','','','Rho') ) )

for (i in 1:numRows){
    if(start <= mutGFP$pos[i] & mutGFP$pos[i] <= end){
    	  j <- mutGFP$pos[i] - start +1
    	  pooledCountMatrix[j, ] <- as.matrix(mutGFP[i,4:8]) + t(as.matrix(pooledCountMatrix[j, 1:5]))
	  
    }
}
#head(as.matrix(pooledCountMatrix))

for (i in 1:length(pooledCountMatrix[,2])){
    correlationMatrix[i,5] <- cor(1:4,  pooledCountMatrix[i,2:5], method='spearman')
}

totalsCounts <- apply(pooledCountMatrix, 1, sum)
norm.pooledCountMatrix <-  pooledCountMatrix/totalsCounts

logTots <- matrix(0.8,nrow = end-start+1, ncol = 5,
		  dimnames = list(start:end,   c('','','','','total')) )
logTots[,5] <- as.matrix(log10(totalsCounts))


#pdf(file=    "heatmap-readCounts-gfp.pdf", width = 5,  height = 15)
heatmap.2(norm.pooledCountMatrix,
	col=rev(brewer.pal(n = 11, name = "RdBu")), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=seq(0,0.54,length=11+1), scale="none",
	cexRow=0.5, cexCol=7.0, margins = c(30, 2), key.title = "Fraction reads", Colv=FALSE, Rowv=FALSE,main="GFP")
#dev.off()

#pdf(file=    "heatmap-logTotals-gfp.pdf", width = 5,  height = 15)
heatmap.2(logTots,
	col=c('white',rev(brewer.pal(n = 11, name = "Spectral"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(0.8,seq(1,4,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, margins = c(30, 2), key.title = "log10(Tot. Reads)", Colv=FALSE, Rowv=FALSE,main="GFP")
#dev.off()

#pdf(file=    "heatmap-spearmanRho-gfp.pdf", width = 5,  height = 15)
heatmap.2(correlationMatrix,
	col=c('white',rev(brewer.pal(n = 11, name = "PRGn"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(-1.1,seq(-1,1,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, margins = c(30, 2), key.title = "Rho", Colv=FALSE, Rowv=FALSE,main="GFP")
#dev.off()


######################################################################

######################################################################
#library('gplots')

cocBroc <- read.table("broccoli.sam.co-occurance", sep="\t", header=TRUE, row.names=1)
colnames(cocBroc) <- 1:201

#pdf(file=    "heatmap-mutation-co-occurence-broc.pdf", width = 15,  height = 15)
par(cex=5.0)
heatmap.2(as.matrix(log10(cocBroc + 0.1)),
	col=c('white',rev(brewer.pal(n = 11, name = "Spectral"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(-1,seq(0,3.6,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, notecex=5.0, margins = c(30, 2), key.title = "log10(mutation co-occurance)", Colv=FALSE, Rowv=FALSE,
	main = 'Brocolli mutation co-occurrence matrix')
#dev.off()


cocGFP <- read.table("egfp.sam.co-occurance", sep="\t", header=TRUE, row.names=1)
colnames(cocGFP) <- 1:201


#pdf(file=    "heatmap-mutation-co-occurence-GFP.pdf", width = 15,  height = 15)
par(cex=5.0)
heatmap.2(as.matrix(log10(cocGFP + 0.1)),
	col=c('white',rev(brewer.pal(n = 11, name = "Spectral"))), density.info="none", trace="none", dendrogram="none",
	symm=F,symkey=F,symbreaks=T, breaks=c(-1,seq(0,3.6,length=11+1)), scale="none",
	cexRow=0.5, cexCol=7.0, notecex=5.0, margins = c(30, 2), key.title = "log10(mutation co-occurance)", Colv=FALSE, Rowv=FALSE,
	main = 'GFP mutation co-occurrence matrix'
	)
dev.off()
######################################################################


#⫸(primer: GFP1_091) — 311 bp — (primer: GFPR2)⫷         Q = 940.6  (moderate amplification)

#Amplified sequence: 
#

# gaacaccgtggtgagaatccaagctagcgTCGACATAAAAAATTTATTTGCTTTGTGAGCGGATAACAATTATAATAGATTCAATTGTGAGCGGATAACAATTTCACACAGAATTCATTAAAGAGGAGAAATTAACTCAT
# ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCAC
# AAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGtgttcaatgctttgcgagataccca

#                                                                                                                                             |START<--                                                         -->MUTNS|
# gaacaccgtggtgagaatccaagctagcgTCGACATAAAAAATTTATTTGCTTTGTGAGCGGATAACAATTATAATAGATTCAATTGTGAGCGGATAACAATTTCACACAGAATTCATTAAAGAGGAGAAATTAACTCATATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGtgttcaatgctttgcgagataccca
#                                                                                                                                                                                                                                                        -->|
# 																														         250
#                                                                                                                    |<--
#                                                                                                                   250



# > predictions.GFP[ predictions.GFP == 'GoF' ]
#  10  14  22  34  79  96 218 219 220 221 223 224 226 227 
# GoF GoF GoF GoF GoF GoF GoF GoF GoF GoF GoF GoF GoF GoF 

#1  4  7  10 13 16 19 22 25
#ATGCGTAAAGGAGAAGAACTTTTCACT
#M  R  K  G  E  E  L  F  T
#GoF:
#   S  T  G              T
#   P

#GAAGAACTTTTCACT



#MRKGEELFT
