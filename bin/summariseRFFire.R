#!/usr/bin/Rscript

#R CMD BATCH summariseRFFire.R
#Read results of multiple random forest classifications of GoF/LoF predictions for Brocolli and GFP. 
#Plot results in a *nice* way. 
#
require(graphics)
library(scales)
######################################################################
gofBrocFire <- read.table(file="gofBroc-appended2.Rdata",sep=",",header=T)
lofBrocFire <- read.table(file="lofBroc-appended2.Rdata",sep=",",header=T)
gofGFPFire  <- read.table(file= "gofGFP-appended2.Rdata",sep=",",header=T)
lofGFPFire  <- read.table(file= "lofGFP-appended2.Rdata",sep=",",header=T)
######################################################################
brocProbs <- read.table(file="brocProbs-appended.Rdata",sep=",",header=T)
 gfpProbs <- read.table(file="gfpProbs-appended.Rdata",sep=",",header=T)
######################################################################
brocCnts <- read.table(file="lofBroc-counts.tsv",sep="\t",header=F)
gfpCnts  <- read.table(file= "lofGFP-counts.tsv",sep="\t",header=F)
######################################################################



######################################################################
paired.cols <- palette.colors(palette='Paired')
cols = c("lightgray", paired.cols[c(5,1,9,7)])    #c("lightgray","pink", "lightblue","orchid","orange") #neutral, GoF, LoF, Uniform, Other
#SEQUENCE PLOTS:
nuc2num <- function( nuc ){
    nuc <- toupper(nuc);
    if(nuc == 'A'){
       return(1);
    }else if(nuc == 'C'){
       return(2);
    }else if(nuc == 'G'){
       return(3);
    }else if(nuc == 'T'){
       return(4);
    } else {
      return(NA);
    }
}
######################################################################
pdf(file="summariseRFFire-plots.pdf", width=20, height=10)

######################################################################
#4x motif
brocSeq <- 'GAGACGGUCGGGUCCAUCUGAGACGGUCGGGUCCAGAUAUUCGUAUCUGUCGAGUAGAGUGUGGGCUCAGAUGUCGAGUAGAGUGUGGGCUCCCACAUACUCUGAUGAUCCAGACGGUCGGGUCCAUCUGAGACGGUCGGGUCCAGAUAUUCGUAUCUGUCGAGUAGAGUGUGGGCUCAGAUGUCGAGUAGAGUGUGGGCU'
struct  <- '(((...........((((((((...........((((((....))))))................))))))))................)))............(((((((((...........((((((((...........((((((....))))))................))))))))................)))))))))....'
#struct  <- '(((...........((((((((...........((((((....))))))................))))))))................)))...................((...........((((((((...........((((((....))))))................))))))))................))'

structA <- unlist(strsplit(struct,""))
brocSeqA<- unlist(strsplit(brocSeq,""))
prev    <- 0
numRecClassifications <- 5
broc.p.lof.mean <- c()
###################################
GFPSeq <- 'AUGCGUAAAGGAGAAGAACUUUUCACUGGAGUUGUCCCAAUUCUUGUUGAAUUAGAUGGUGAUGUUAAUGGGCACAAAUUUUCUGUCAGUGGAGAGGGUGAAGGUGAUGCAACAUACGGAAAACUUACCCUUAAAUUUAUUUGCACUACUGGAAAACUACCUGUUCCAUGGCCAACACUUGUCACUACUUUCGGUUAUGGU'
GFPSeqA<- unlist(strsplit(GFPSeq,""))
prot <- 'MRKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYG'
protA   <- unlist(strsplit(prot,""))


######################################################################
#LoF plots
##Brocolli:
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1, cex=2) #mar: bottom, left, top, right
maxP    <- 1.0 #max(brocProbs[ brocProbs[,5] == 2, 6])
lineSpace <- 0.03
bottomLine<- -0.05
plot(NA,NA, ylim=c(0,maxP+4*lineSpace),xlim=c(8,length(brocSeqA)+1),xlab="", ylab="", main="Broccoli loss of function variants", xaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
for (i in 0:10){
    lines(c(-10,length(brocSeqA)+1), c( i/10, i/10), col='lightgray')
}
lofBroc <- brocProbs[ brocProbs[,5] == 2, ]
uniqLoFBroc <- rle(sort(lofBroc[,2])) #uniq -c: return unique'ed values, and their counts
for (i in 1:length(uniqLoFBroc$lengths) ){
    if( uniqLoFBroc$lengths[i] > numRecClassifications ){
    	#Enough RFs classed this site as LoF:
	lofIslice <- lofBroc[ lofBroc[,2] == uniqLoFBroc$values[i], ]
	for (j in 1:length(lofIslice[,1]) ){	    
	    #Find x-value:	
	    offSet <- nuc2num( lofIslice[j,4] )
    	    xPos   <- lofIslice[j,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	    points(  xPos,               lofIslice[j,6] , pch=20,col=cols[2])	
	}
	mn <- mean( lofIslice[,6])
    	text(    xPos, mn, lofIslice[1,4], pos=3, cex=0.80)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
text(180, maxP-5*lineSpace, paste('N=',length(uniqLoFBroc$lengths),"/603", sep=''), pos=3, cex=1.2, font=2)
for (i in 1:length(brocSeqA)){
    text( i,    bottomLine+lineSpace,brocSeqA[i], pos=3, cex=0.35) #real sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( i,    bottomLine,          structA[i],  pos=3, cex=0.40) #structure
    text( i,    maxP+1*lineSpace,    brocSeqA[i], pos=3, cex=0.35) #structure
    text( i,    maxP+0*lineSpace,    structA[i],  pos=3, cex=0.40) #structure
    if(i %% 10 == 0){
    	 text( i, bottomLine+2*lineSpace, i, pos=3, cex=0.60, col='black')
    	 text( i,       maxP+2*lineSpace, i, pos=3, cex=0.60, col='black')
	 lines(c(i,i), c( -2*lineSpace, 1.00), col='lightgray')
    }
}

##GFP
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1, cex=2) #mar: bottom, left, top, right
maxP    <- 1.0 #max(gfpProbs[ gfpProbs[,5] == 2, 6])
lineSpace <- 0.015
plot(NA,NA, ylim=c(0,maxP+4*lineSpace),xlim=c(8,length(GFPSeqA)+1),xlab="", ylab="", main="GFP loss of function variants", xaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
for (i in 0:10){
    lines(c(-10,length(GFPSeqA)+1), c( i/10, i/10), col='lightgray')
}
lofGFP <- gfpProbs[ gfpProbs[,5] == 2, ]
uniqLoFGFP <- rle(sort(lofGFP[,2])) #uniq -c: return unique'ed values, and their counts
for (i in 1:length(uniqLoFGFP$lengths) ){
    if( uniqLoFGFP$lengths[i] > numRecClassifications ){
    	#Enough RFs classed this site as LoF:
	lofIslice <- lofGFP[ lofGFP[,2] == uniqLoFGFP$values[i], ]
	for (j in 1:length(lofIslice[,1]) ){	    
	    #Find x-value:	
	    offSet <- nuc2num( lofIslice[j,4] )
    	    xPos   <- lofIslice[j,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	    points(  xPos,               lofIslice[j,6] , pch=20,col=cols[2])	
	}
	mn <- mean( lofIslice[,6])
    	text(    xPos, mn, lofIslice[1,4], pos=3, cex=0.80)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
text(180, maxP-5*lineSpace, paste('N=',length(uniqLoFGFP$lengths),"/603", sep=''), pos=3, cex=1.2, font=2)
j <- 1
for (i in 1:length(GFPSeqA)){
    text( i,  bottomLine+lineSpace,   GFPSeqA[i],  pos=3, cex=0.35) #sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( i,        maxP+lineSpace,   GFPSeqA[i],  pos=3, cex=0.35) #sequence
    if(i %% 10 == 0){
    	 text( i, bottomLine+2*lineSpace, i, pos=3, cex=0.60, col='black')
    	 text( i,       maxP+2*lineSpace, i, pos=3, cex=0.60, col='black')
	 lines(c(i,i), c(   -2*lineSpace, 1.00), col='lightgray')
    }    
    if(i %% 3 == 1){
    	 text( i,               maxP,    protA[j],  pos=3, cex=0.60) #structure
    	 text( i, bottomLine - 1*lineSpace,    protA[j],  pos=3, cex=0.60) #structure
	 j <- j+1
    }
}



######################################################################
#GoF/LoF plots
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1) #mar: bottom, left, top, right
maxP    <- max(c(gofBrocFire$GoF,lofBrocFire$LoF))
plot(NA,NA, ylim=c(-1*maxP,maxP),xlim=c(8,length(brocSeqA)+1),xlab="", ylab="", main="Broccoli sequence", xaxt = "n", yaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
axis(2, at=seq(-0.6,0.6,by=0.2), labels=abs(seq(-0.6,0.6,by=0.2)))
lines(c(-10,length(brocSeqA)+1), c( 0.00, 0.00), col='lightgray')
lines(c(-10,length(brocSeqA)+1), c( 0.50, 0.50), col='lightgray')
lines(c(-10,length(brocSeqA)+1), c(-0.50,-0.50), col='lightgray')
for (i in 1:length(gofBrocFire[,1])){
    if( length(gofBrocFire[gofBrocFire$pos == gofBrocFire[i,2],1]) > numRecClassifications ){
    	offSet <- nuc2num( gofBrocFire[i,4] )
    	xPos   <- gofBrocFire[i,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	points(  xPos,               gofBrocFire[i,12] , pch=20,col=cols[2])
    	mn <- mean(gofBrocFire[gofBrocFire$pos == gofBrocFire[i,2],12])
    	text(    xPos, mn, gofBrocFire[i,4], pos=3, cex=0.60)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
for (i in 1:length(lofBrocFire[,1])){
    if( length(lofBrocFire[lofBrocFire$pos == lofBrocFire[i,2],1]) > numRecClassifications ){
    	offSet <- nuc2num( lofBrocFire[i,4] )
    	xPos   <- lofBrocFire[i,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	points(  xPos,               -1*lofBrocFire[i,12] , pch=20,col=cols[3])
    	mn <- -1*mean(lofBrocFire[lofBrocFire$pos == lofBrocFire[i,2],12])
	broc.p.lof.mean <- c(broc.p.lof.mean, -1*mn)
    	text(    xPos, mn, lofBrocFire[i,4], pos=3, cex=0.60)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
text(200, 0.9*maxP,'GoF',col=cols[2], font=2, cex=3)
text(200,-0.9*maxP,'LoF',col=cols[3], font=2, cex=3)
for (i in 1:length(brocSeqA)){   
    text( i,    0.00,    brocSeqA[i], pos=3, cex=0.60) #real sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( i,    maxP,    structA[i],  pos=3, cex=0.60) #structure
    text( i,   -0.05,    structA[i],  pos=3, cex=0.60) #structure
    text( i, -1*maxP,    structA[i],  pos=3, cex=0.60) #structure
    if(i %% 10 == 0){
    	 text( i, -0.02, i, pos=3, cex=0.60, col='black')
    	 points(i,0,pch='|',col='lightgray')
	 lines(c(i,i), c( -1.00, 1.00), col='lightgray')
    }
}
######################################################################

gfp.p.lof.mean <- c()
maxP    <- max(c(gofGFPFire$GoF,lofGFPFire$LoF))
prev    <- 0
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1) #mar: bottom, left, top, right
plot(NA,NA, ylim=c(-1*maxP,maxP),xlim=c(8,length(GFPSeqA)+1),xlab="", ylab="", main="GFP sequence", xaxt = "n", yaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
axis(2, at=seq(-0.6,0.6,by=0.2), labels=abs(seq(-0.6,0.6,by=0.2)))
lines(c(-10,length(GFPSeqA)+1), c( 0.00, 0.00), col='lightgray')
lines(c(-10,length(GFPSeqA)+1), c( 0.50, 0.50), col='lightgray')
lines(c(-10,length(GFPSeqA)+1), c(-0.50,-0.50), col='lightgray')
for (i in 1:length(gofGFPFire[,1])){
    if( length(gofGFPFire[gofGFPFire$pos == gofGFPFire[i,2],1]) > numRecClassifications ){
    	offSet <- nuc2num( gofGFPFire[i,4] )
    	xPos   <- gofGFPFire[i,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	points(  xPos,               gofGFPFire[i,12] , pch=20,col=cols[2])
    	mn <- mean(gofGFPFire[gofGFPFire$pos == gofGFPFire[i,2],12])
    	text(    xPos, mn, gofGFPFire[i,4], pos=3, cex=0.60)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
for (i in 1:length(lofGFPFire[,1])){
    if( length(lofGFPFire[lofGFPFire$pos == lofGFPFire[i,2],1]) > numRecClassifications ){
    	offSet <- nuc2num( lofGFPFire[i,4] )
    	xPos   <- lofGFPFire[i,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	points(  xPos,               -1*lofGFPFire[i,12] , pch=20,col=cols[3])
    	mn <- -1*mean(lofGFPFire[lofGFPFire$pos == lofGFPFire[i,2],12])
	gfp.p.lof.mean <- c(gfp.p.lof.mean, -1*mn)
    	text(    xPos, mn, lofGFPFire[i,4], pos=3, cex=0.60)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')
    }
}
text(200, 0.9*maxP,'GoF',col=cols[2], font=2, cex=3)
text(200,-0.9*maxP,'LoF',col=cols[3], font=2, cex=3)
j=1
for (i in 1:length(GFPSeqA)){   
    text( i,    0.00,    GFPSeqA[i], pos=3, cex=0.60) #real sequence #pos: 1-below, 2-left, 3-above, 4-right
    
    if(i %% 3 == 1){
    	 text( i,    maxP,    protA[j],  pos=3, cex=0.60) #structure
    	 text( i,   -0.05,    protA[j],  pos=3, cex=0.60) #structure
    	 text( i, -1*maxP,    protA[j],  pos=3, cex=0.60) #structure
	 j <- j+1
    }
    
    if(i %% 10 == 0){
    	 text( i, -0.02, i, pos=3, cex=0.60, col='black')
    	 points(i,0,pch='|',col='lightgray')
	 lines(c(i,i), c( -1.00, 1.00), col='lightgray')
    }
}
######################################################################

require(graphics)
paired.cols <- palette.colors(palette='Paired')
brocCol <- paired.cols[3] #'seagreen1'
gfpCol  <- paired.cols[4] #'darkgreen'


par(las=0,bty='l',mfrow=c(1,2),mar=c(4, 5, 4, 1)+0.1) #mar: bottom, left, top, right
breaks <- seq(0,1,length.out=40)
hist(broc.p.lof.mean, freq=TRUE, breaks = breaks, main="Broccoli: recurrent LoF classification", xlab="Mean LoF classification probability", ,col=brocCol)
hist( gfp.p.lof.mean, freq=TRUE, breaks = breaks, main="GFP: recurrent LoF classification"     , xlab="Mean LoF classification probability",col= gfpCol)


#######################
#REPORT SUMMARY STATISTICS
#MEDIAN, IQR, confidence interval???
par(las=0,bty='l',mfrow=c(1,2),mar=c(4, 5, 4, 1)+0.1,cex=2.5,las=2) #mar: bottom, left, top, right
col <- c("orchid1", "skyblue", "lightpink", "lightblue1" )
lofRP <- cbind(brocCnts, gfpCnts) 
boxplot(lofRP, col=col, ylab="LoF variants", xlab="", main="LoF variants GFP vs Brocolli", outline=FALSE, xaxt = "n")
axis(1, 1:2, c('Brocolli','GFP'))
cols<-c("deeppink4","dodgerblue4")
for (i in 1:2){
    xp <-  rnorm(length(lofRP[,i]), mean = 0, sd = 0.1)   
    points(i+xp, lofRP[,i],col = "white", bg=cols[i], pch=21, cex=0.75)
}

#CHECK, are these the mean of all 600 possible variants?
#ADD A SIGNIFICANCE TEST!!!!!!!!!
lofpRP <- cbind(broc.p.lof.mean, gfp.p.lof.mean) 
boxplot(lofpRP, col=col, ylab="Mean LoF probabilities", xlab="", main="Mean Probability of LoF", outline=FALSE, xaxt = "n", ylim=c(0.25,0.9))
axis(1, 1:2, c('Brocolli','GFP'))
cols<-c("deeppink4","dodgerblue4")
for (i in 1:2){
    xp <-  rnorm(length(lofpRP[,i]), mean = 0, sd = 0.1)   
    points(i+xp, lofpRP[,i],col = "white", bg=cols[i], pch=21, cex=0.75)
}

dev.off()
############



#ADD A10T notation...
#Add box-whisker plots...


#cols<-c("deeppink4","dodgerblue4","deeppink4","dodgerblue4")
cols <- c("skyblue", "lightblue1", "orchid1", "lightpink", "skyblue", "lightblue1", "skyblue", "lightblue1", "orchid1", "lightpink" )

pdf(file="summariseRFFire-LoF-main.pdf", width=20, height=20)
par(las=0,bty='l',mfrow=c(2,1),mar=c(2, 7, 4, 1)+0.1, cex=2) #mar: bottom, left, top, right
maxP    <- 1.0 #max(brocProbs[ brocProbs[,5] == 2, 6])
minP    <- 0.2
lineSpace <- 0.03
bottomLine<- minP-0.05
################
plot(NA,NA, ylim=c(minP,maxP+4*lineSpace),xlim=c(8,length(GFPSeqA)+1),xlab="", ylab="", main="GFP loss of function variants", xaxt = "n")
axis(1,at=10*(1:20),labels=10*(1:20))
mtext("LoF probability", side = 2, line = 3, cex=2)
for (i in 0:10){
    lines(c(-10,length(GFPSeqA)+1), c( i/10, i/10), col='lightgray')
}
lofGFP <- gfpProbs[ gfpProbs[,5] == 2, ]
uniqLoFGFP <- rle(sort(lofGFP[,2])) #uniq -c: return unique'ed values, and their counts
for (i in 1:length(uniqLoFGFP$lengths) ){
    if( uniqLoFGFP$lengths[i] > numRecClassifications ){
    	#Enough RFs classed this site as LoF:
	lofIslice <- lofGFP[ lofGFP[,2] == uniqLoFGFP$values[i], ]
	for (j in 1:length(lofIslice[,1]) ){	    
	    #Find x-value:	
	    offSet <- nuc2num( lofIslice[j,4] )
    	    xPos   <- lofIslice[j,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	    points(  xPos,               lofIslice[j,6] , pch=20,col=cols[2])	
	}
	mn <- mean( lofIslice[,6])
    	text(    xPos, mn, lofIslice[1,4], pos=3, cex=0.80)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
text(180, maxP-5*lineSpace, paste('N=',length(uniqLoFGFP$lengths),"/603", sep=''), pos=3, cex=1.2, font=2)
j <- 1
for (i in 1:length(GFPSeqA)){
    text( i,  bottomLine+lineSpace,   GFPSeqA[i],  pos=3, cex=0.35) #sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( i,        maxP+lineSpace,   GFPSeqA[i],  pos=3, cex=0.35) #sequence
    if(i %% 10 == 0){
    	 text( i, bottomLine+2*lineSpace, i, pos=3, cex=0.60, col='black')
    	 text( i,       maxP+2*lineSpace, i, pos=3, cex=0.60, col='black')
	 lines(c(i,i), c(   -2*lineSpace, 1.00), col='lightgray')
    }    
    if(i %% 3 == 1){
    	 text( i,       maxP,    protA[j],  pos=3, cex=0.60) #structure
    	 text( i, bottomLine,    protA[j],  pos=3, cex=0.60) #structure
	 j <- j+1
    }
}
################
plot(NA,NA, ylim=c(minP,maxP+4*lineSpace),xlim=c(8,length(brocSeqA)+1),xlab="", ylab="", main="Broccoli loss of function variants", xaxt = "n")
axis(1,at=10*(1:20),labels=10*(1:20))
mtext("LoF probability", side = 2, line = 3, cex=2)
for (i in 0:10){
    lines(c(-10,length(brocSeqA)+1), c( i/10, i/10), col='lightgray')
}
lofBroc <- brocProbs[ brocProbs[,5] == 2, ]
uniqLoFBroc <- rle(sort(lofBroc[,2])) #uniq -c: return unique'ed values, and their counts
for (i in 1:length(uniqLoFBroc$lengths) ){
    if( uniqLoFBroc$lengths[i] > numRecClassifications ){
    	#Enough RFs classed this site as LoF:
	lofIslice <- lofBroc[ lofBroc[,2] == uniqLoFBroc$values[i], ]
	for (j in 1:length(lofIslice[,1]) ){	    
	    #Find x-value:	
	    offSet <- nuc2num( lofIslice[j,4] )
    	    xPos   <- lofIslice[j,2] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    	    points(  xPos,               lofIslice[j,6] , pch=20,col=cols[4])	
	}
	mn <- mean( lofIslice[,6])
    	text(    xPos, mn, lofIslice[1,4], pos=3, cex=0.80)
	lines( c(xPos,xPos), c(0,mn), col=alpha('lightgray',0.5))
	points(  xPos, mn, pch='-',col='black')    
    }
}
text(180, maxP-5*lineSpace, paste('N=',length(uniqLoFBroc$lengths),"/603", sep=''), pos=3, cex=1.2, font=2)
for (i in 1:length(brocSeqA)){
    text( i,    bottomLine+lineSpace,brocSeqA[i], pos=3, cex=0.35) #real sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( i,    bottomLine,          structA[i],  pos=3, cex=0.40) #structure
    text( i,    maxP+1*lineSpace,    brocSeqA[i], pos=3, cex=0.35) #structure
    text( i,    maxP+0*lineSpace,    structA[i],  pos=3, cex=0.40) #structure
    if(i %% 10 == 0){
    	 text( i, bottomLine+2*lineSpace, i, pos=3, cex=0.60, col='black')
    	 text( i,       maxP+2*lineSpace, i, pos=3, cex=0.60, col='black')
	 lines(c(i,i), c( -2*lineSpace, 1.00), col='lightgray')
    }
}
dev.off()
