#!/usr/bin/Rscript

#R CMD BATCH trainClassifier2.R
# Simulate count data sorted by expression level sample from
# distributions corresponding to neutral, GoF (higher fluorescence) or
# LoF (lower fluorescence) and a mix of GoF+LoF (double mutants).

######################################################################
require(graphics)
library(scales)
library(plot.matrix)
library(randomForest)
library('EnvStats') #for the 'revd' function (extreme value/Gumbel distribution)
#Colours: 
#https://www.zeileis.org/news/coloring/

###########
#Function to convert mutation tables from sam2var to a matrix for stacked barplots: 
#--altered to return TPMs, with a pseudocount added. 
mutTable2stackMatrix <- function( mutTable, column ){
     x<-matrix(100/10^6, nrow = 201, ncol = 4, dimnames = list(1:201, c('A','C','G','T'))) ##initialise with a pseudocount
     totReads <- sum(mutTable[, column])
     x[mutTable[mutTable$querynuc == 'A', 1], 1] <- x[mutTable[mutTable$querynuc == 'A', 1], 1] + (10^6)*(mutTable[mutTable$querynuc == 'A', column])/totReads
     x[mutTable[mutTable$querynuc == 'C', 1], 2] <- x[mutTable[mutTable$querynuc == 'C', 1], 2] + (10^6)*(mutTable[mutTable$querynuc == 'C', column])/totReads
     x[mutTable[mutTable$querynuc == 'G', 1], 3] <- x[mutTable[mutTable$querynuc == 'G', 1], 3] + (10^6)*(mutTable[mutTable$querynuc == 'G', column])/totReads
     x[mutTable[mutTable$querynuc == 'T', 1], 4] <- x[mutTable[mutTable$querynuc == 'T', 1], 4] + (10^6)*(mutTable[mutTable$querynuc == 'T', column])/totReads	     
     return(x)
}
######################################################################
# * Four types of variant: high, low, neutral or other (these may be double mutants -- enriched in both the high and low fractions)
# * there are 600 unique (singleton) variants, each from one of the above types
# * presorting, each variant has a different frequency, this distribution is exponential/Gumbel
# * sorting is based on both the proportion in the presort -- and which type the variant is -- sample roughly the same number of variants from each bin! -- this is equivalent to how PCR+sequencing kept the counts roughly the same
# * convert all counts to TPMs within each bin, compute log-odds scores relative to presorting for each variant in each bin
######################################################################
#Read in Broccoli & GFP data, construct dataframes for RF:
mutBroc <- read.table("sam2var.broccoli-118-318.txt",   sep="\t", header=TRUE)
mutGFP  <- read.table("sam2var.egfp-2410-2610.txt", sep="\t", header=TRUE)
#Filter any sites with a low sample size
mutBroc <- mutBroc[mutBroc$total>9,] 
mutGFP  <-  mutGFP[ mutGFP$total>9,]
presort.tpm.broc <- mutTable2stackMatrix(mutBroc, 4)
presort.tpm.gfp  <- mutTable2stackMatrix(mutGFP,  4)
#Summary stats of presort data (remove the zeros -- these are not useful for the analysis):
pdf(file="trainClassifier-plots.pdf", width=20, height=10)
combined.tpm <- c(presort.tpm.broc, presort.tpm.gfp)
hist( combined.tpm, xlab='TPM', xlim=c(0,20000), col='green', breaks=80, main='Presorted GFP & Brocolli'     )
text(15000,500, paste("Min=", round(min(combined.tpm),0),"\nMedian=",round(median(combined.tpm),0),"\nMean=",round(mean(combined.tpm),0),"\nMax=", round(max(combined.tpm),0), sep=""),pos=2)

x<-log10(c(presort.tpm.broc, presort.tpm.gfp)+1)
x<-x[x>1]
(mn.presort <- mean(x))
(sd.presort <-   sd(x))
######################################################################
#PARAMETERS:
reps   <- 20 #100 #0                            # number of times to repeat the experiment 
N      <- 600                                   # number of "variable" sites, 120 per class (GoF, LoF, Neutral, Uniform, Other)
R.tot  <- 30000                                 # number of reads 
n      <- c(15, 100)                            # range of total read counts for broccoli & GFP
pcount <- 1.50                                  # regularise the TPM logodds values 
presortParams <- c(  mn.presort,  sd.presort  ) # generate presort log10 TPMs: values determine the probability of variant being
                                                # sampled from presort population 
neutralParams <- c( -1.15,  1.15, 0.75, 1.25  ) # ranges for mean & standard deviation
    gofParams <- c(  1.20,  3.00, 0.80, 1.20  ) 
    lofParams <- c( -3.00, -1.20, 0.80, 1.20  )
    uniParams <- c( -5.00, -4.00, 4.00, 5.00  ) # Max & Min
  otherParams <- c( -4.50, -2.50, 0.20, 0.70  ) # location, scale (evd)
######################################################################
# Workflow:
#  1. sample variant/read from presort population 
#  2. determine expression level, conditioned on variant type
#  3. determine which bin variant is sorted into
#  4. Re-size the totals for each bin to be roughly the constant
#  5. Compute TPMs & log-odds values for each variant in each bin.....       ADD FEATURES?: DROP TPMs, ADD spearman(bins vs log-odds), add TOTAL reads 
#  6. REPEAT 1-5 rep times
######################################################################
#sample2features: input expression levels (0 - 1), sort values into bins lo, lo-med, med-hi, hi
sample2features <- function( y, label ){
    features      <- matrix(as.numeric(0.00), nrow = 1, ncol = 5)
    features[1,1] <- length( y[  y <= 0.025             ] ) 
    features[1,2] <- length( y[ 0.075 <= y & y <= 0.475 ] ) 
    features[1,3] <- length( y[ 0.525 <= y & y <= 0.925 ] ) 
    features[1,4] <- length( y[ 0.975 <= y              ] ) 
    features[1,5] <- label
    return(features)
}
######################################################################
#
corBin <- function( y ){
       x <- 1:4 
       return(cor(x,y, method="spearman"))
}
######################################################################
#counts2rf.features: convert read counts to a features data.frame 
counts2rf.features <- function( bin.counts, background.tpm ){
     #For each bin (Low, Low.Med, Med.High, High) compute TPMs and log-ratios with the presort 
     write.csv(bin.counts, file="bin.counts.Rdata")
     #write.csv(background.tpm, file="background.tpm.Rdata")
     #write.csv(c(length(bin.counts$low), length(bin.counts$low.med), length(bin.counts$med.high), length(bin.counts$high), length(background.tpm), length(bin.counts$class), pcount), file="lengths.Rdata" )
     low.tpm          <- (10^6)*bin.counts$low      /sum(bin.counts$low     )
     low.med.tpm      <- (10^6)*bin.counts$low.med  /sum(bin.counts$low.med )
     med.high.tpm     <- (10^6)*bin.counts$med.high /sum(bin.counts$med.high)
     high.tpm         <- (10^6)*bin.counts$high     /sum(bin.counts$high    )
     low.logodds      <- log10( ( low.tpm     +pcount) / ( background.tpm+pcount) )
     low.med.logodds  <- log10( ( low.med.tpm +pcount) / ( background.tpm+pcount) )
     med.high.logodds <- log10( ( med.high.tpm+pcount) / ( background.tpm+pcount) )
     high.logodds     <- log10( ( high.tpm    +pcount) / ( background.tpm+pcount) )

     corBins          <- apply( cbind(low.logodds, cbind(low.med.logodds, cbind(med.high.logodds, high.logodds))),
	 	1,
		corBin)

     write.csv(corBins, file="corBins.Rdata")
     rf.sim.features <- data.frame(
         low.tpm          = low.tpm,
         low.med.tpm      = low.med.tpm,
         med.high.tpm     = med.high.tpm,
         high.tpm         = high.tpm,
	 low.logodds      = low.logodds,
	 low.med.logodds  = low.med.logodds,
	 med.high.logodds = med.high.logodds,
	 high.logodds     = high.logodds,
	 correlation      = corBins,
	 total            = bin.counts$low + bin.counts$low.med + bin.counts$med.high + bin.counts$high
     )
     columnNames <- colnames(rf.sim.features)
     if(length(bin.counts$class) == length(bin.counts$low)){
	rf.sim.features<-cbind(rf.sim.features, bin.counts$class)	
     	colnames(rf.sim.features) <- c(columnNames, 'class')
     }
     #write.csv(rf.sim.features, file="rf.sim.features.Rdata")
     return(rf.sim.features)
}
######################################################################
##PLOT INITIALISING
par(cex=2)
plot(NA,NA,xlim=c(-5,5), ylim=c(0,0.60),xlab="Normalised Expression",ylab="Freq.",main="Simulated training and test data")
numBreaks <- 60
x <- seq(-7.3, 7.3, length=numBreaks)
sum.neutral <- 0*x; sum.gof <- 0*x; sum.lof <- 0*x; sum.uni <- 0*x; sum.other <- 0*x; 
paired.cols <- palette.colors(palette='Paired')
cols = c("lightgray", paired.cols[c(5,1,9,7)])    #c("lightgray","pink", "lightblue","orchid","orange") #neutral, GoF, LoF, Uniform, Other
brocCol <- paired.cols[3] #'seagreen1'
gfpCol  <- paired.cols[4] #'darkgreen'
######################################################################
##REPLICATE EXPERIMENT "rep" times!!!!
simulation.features  <- c()
for (j in 1:reps){
  ###################################
  #Probabilities of  variant being sampled
  presort.sim.logTPMs  <- rnorm(N, mn.presort, sd.presort)
  # presort.sim.dens     <- dnorm(presort.sim.logTPMs, mean = mn.presort, sd = sd.presort)
  # convert simulated logTPMs to # of reads & verify the original dataset roughly matches, and get the TPMs:
  #presort.sim.reads <-   ceiling(  10^(presort.sim.logTPMs)  - 100/10^6 ) 
  #sumReads <- sum(presort.sim.reads)
  presort.sim.TPMs <- 10^(presort.sim.logTPMs) 
  #sum(presort.sim.reads)
  #write.csv(presort.sim.reads, "presort.sim.reads.Rdata")
  #Sample proportional number of reads from presort, allocate to sorted bins across neutral, LoF, GoF & other classes
  presort.sim.probs <- presort.sim.TPMs / sum(presort.sim.TPMs)    #(presort.sim.reads+1) / sum(c(presort.sim.reads, N))
  write.csv(presort.sim.probs, "presort.sim.probs.Rdata")
  #add some noise to R.tot
  numReads <- round(rnorm(1, mean = R.tot, sd = 1000),0)
  presort.sample <- rle(sort(c(1:N,  sample(1:(N+1),4*numReads,replace=TRUE,prob=c(presort.sim.probs, 10^(-3)) )) ) )
  write.csv(presort.sample$lengths, "presort.sample.Rdata")

  sim.bins <- data.frame( low = integer(N), low.med = integer(N), med.high = integer(N),  high = integer(N), class = rep('',N)  )
  ###################################
  for (i in 1:N ){
    #split "reads" across bins, based on distribution Params
    if(is.na(presort.sample$lengths[i]) == TRUE){
        next;
    }
    if( i <= floor(N/5)){
    	#NEUTRAL    
    	mn.neutral <- runif(1, min = neutralParams[1], max = neutralParams[2])
    	sd.neutral <- runif(1, min = neutralParams[3], max = neutralParams[4])
    	numReads.neutral <- presort.sample$lengths[i]
	expression.probs.neutral <- pnorm(rnorm(numReads.neutral, mean = mn.neutral, sd = sd.neutral), mean = 0, sd = 1)
	sim.bins[i,1:5] <- sample2features( expression.probs.neutral, 'Neutral'  )	
	##plot:
	y.neutral  <- dnorm(x,mean=mn.neutral,sd=sd.neutral)
	sum.neutral <- sum.neutral + y.neutral
	lines(x,      y.neutral           ,col=alpha(cols[1],0.1))
    } else if( floor(N/5) < i && i <= 2*floor(N/5)){
    	#GoF    
    	mn.gof <- runif(1, min = gofParams[1], max = gofParams[2])
    	sd.gof <- runif(1, min = gofParams[3], max = gofParams[4])
    	numReads.gof <- presort.sample$lengths[i]
	expression.probs.gof <- pnorm(rnorm(numReads.gof, mean = mn.gof, sd = sd.gof), mean = 0, sd = 1)
	sim.bins[i,1:5] <- sample2features( expression.probs.gof, 'GoF'  )	
	##plot:
	y.gof   <- dnorm(x,mean=mn.gof,sd=sd.gof)
    	sum.gof <- sum.gof + y.gof
	lines(x,      y.gof            ,col=alpha(cols[2],0.1))
    } else if( 2*floor(N/5) < i && i <= 3*floor(N/5)){
    	#LoF    
    	mn.lof <- runif(1, min = lofParams[1], max = lofParams[2])
    	sd.lof <- runif(1, min = lofParams[3], max = lofParams[4])
    	numReads.lof <- presort.sample$lengths[i]
	sample.lof <- rnorm(numReads.lof, mean = mn.lof, sd = sd.lof)
	expression.probs.lof <- pnorm(sample.lof, mean = 0, sd = 1)
	sim.bins[i,1:5] <- sample2features( expression.probs.lof, 'LoF'  )	
	##plot:
	y.lof   <- dnorm(x,mean=mn.lof,sd=sd.lof)
    	sum.lof <- sum.lof + y.lof
	lines(x,      y.lof            ,col=alpha(cols[3],0.1))
    } else if( 3*floor(N/5) < i && i <= 4*floor(N/5)){
    	#Uniform    
    	numReads.uni         <- presort.sample$lengths[i]
	mn.uni               <- runif(           1, min = uniParams[1], max = uniParams[2])
	mx.uni               <- runif(           1, min = uniParams[3], max = uniParams[4])
	sample.uni           <- runif(numReads.uni, min = mn.uni,       max = mx.uni)
	expression.probs.uni <- punif(   sample.uni,min = mn.uni,       max = mx.uni)
	sim.bins[i,1:5] <- sample2features( expression.probs.uni, 'Uniform'  )	
	##plot:
	y.uni   <- dunif(x,min = mn.uni, max = mx.uni)
    	sum.uni <- sum.uni + y.uni
	#lines(x,      4*y.uni            ,col=alpha(cols[4],0.1))
    } else if( 4*floor(N/5) < i ){
    	#Other (bimodal -- High AND Low enriched)
    	numReads.other <- floor(presort.sample$lengths[i]/2)
	loc.other <- runif(1, min = otherParams[1], max = otherParams[2])
	scl.other <- runif(1, min = otherParams[3], max = otherParams[4])
	#write.csv(c(i, presort.sample$lengths[i], numReads.other, loc.other, scl.other), file="revd.other.params.Rdata")
	r.other <- c(   revd(numReads.other,  location = loc.other, scale = scl.other),
		     -1*revd(numReads.other,  location = loc.other, scale = scl.other))
	expression.probs.other <- pnorm( r.other, mean = 0, sd = 1 )
	sim.bins[i,1:5] <- sample2features( expression.probs.other, 'Other'  )
	##plot:
	y.other   <- devd(   x,  location = loc.other, scale = scl.other) +
		     devd(-1*x,  location = loc.other, scale = scl.other)
    	sum.other <- sum.other + y.other
	scale.other <- 2*max(y.other)
	lines(x, y.other/scale.other     ,col=alpha(cols[5],0.025))
    }
  }
  #Save sim.bins features here:
  #First fix F'ing R datatypes 
  for (i in 1:4){
    sim.bins[,i] <- as.numeric( sim.bins[,i] )
  }
  ##APPEND NEW FEATURES TO "simulation.features"
  simulation.features  <- rbind(simulation.features,
  		                counts2rf.features(sim.bins, presort.sim.TPMs))
  write.csv(simulation.features, file="simulation.features.Rdata")
  
}#CLOSE REPS LOOP

#####PLOT AVERAGES,
#paired.cols[c(5,1,9,7)]
lines(x,  0.4*sum.neutral/max(sum.neutral), lwd=2, col='black'        )
lines(x,  0.4*sum.gof    /max(sum.gof    ), lwd=2, col=paired.cols[6] )
lines(x,  0.4*sum.lof    /max(sum.lof    ), lwd=2, col=paired.cols[2] )
lines(x,  0.4*sum.uni    /max(sum.uni    ), lwd=2, col=paired.cols[10])
lines(x,  0.4*sum.other  /max(sum.other  ), lwd=2, col=paired.cols[8] )
legend(4.0,0.60, c("Neutral", "GoF", "LoF", "Uniform", "Other"), col=cols, fill=cols)
##PLOT CELL SORTING BIN RECTANGLES: 
pVals    <- c(0.00001, 0.02500, 0.07500, 0.47500, 0.52500, 0.92500, 0.97500, 0.99999)
binNames <- c('Low', 'Low-Mid.', 'Mid.High', 'High')
text(-4,0.5, 'Bins:',      pos=3, font=2)
for(i in 1:4){
  x1 <- qnorm(pVals[2*i-1], mean = 0, sd = 1)
  x2 <- qnorm(pVals[2*i  ], mean = 0, sd = 1)
  rect(x1, 0, x2, 0.5); xName <- x1+(x2-x1)/2
  text(xName,  0.5, binNames[i], pos=3)
}
######################################################################
#CONSTRUCT BROCOLLI & GFP DATAFRAMES FOR RF:
tpmBrocPresort <- (10^6)*(mutBroc$presort) / sum(mutBroc$presort)
mutBroc.features <- counts2rf.features( mutBroc, tpmBrocPresort )
tpmGFPPresort <- (10^6)*(mutGFP$presort) / sum(mutGFP$presort)
mutGFP.features <- counts2rf.features( mutGFP, tpmGFPPresort )
######################################################################
#PLOT HISTOGRAMS FOR EACH FEATURE AND CLASS:
features2hists <- function(features.df, tpm.breaks, logodds.breaks, col, label){
   tpm.xlim <- c(0,20000)
   hist(  features.df$low.tpm     , xlim=tpm.xlim, xlab='TPM', col=col, breaks=tpm.breaks, main='Low'     )
   mtext(label, side = 2, line = 2, col=col)
   hist(  features.df$low.med.tpm , xlim=tpm.xlim, xlab='TPM', col=col, breaks=tpm.breaks, main='Low-Mid' )
   hist(  features.df$med.high.tpm, xlim=tpm.xlim, xlab='TPM', col=col, breaks=tpm.breaks, main='Mid-High')
   hist(  features.df$high.tpm    , xlim=tpm.xlim, xlab='TPM', col=col, breaks=tpm.breaks, main='High'    )
   logodds.xlim <- c(min(logodds.breaks),max(logodds.breaks))
   hist(  features.df$low.logodds     , xlim=logodds.xlim, xlab='log odds sort/presort', col=col, breaks=logodds.breaks, main='Low'     )
   hist(  features.df$low.med.logodds , xlim=logodds.xlim, xlab='log odds sort/presort', col=col, breaks=logodds.breaks, main='Low-Mid' )
   hist(  features.df$med.high.logodds, xlim=logodds.xlim, xlab='log odds sort/presort', col=col, breaks=logodds.breaks, main='Mid-High')
   hist(  features.df$high.logodds    , xlim=logodds.xlim, xlab='log odds sort/presort', col=col, breaks=logodds.breaks, main='High'    )
   mtext(label, side = 4, line = 1, col=col)
}
#BREAKS:
tpm.all      <- c(simulation.features$low.tpm,     simulation.features$low.med.tpm,      simulation.features$med.high.tpm,     simulation.features$high.tpm,     mutBroc.features$low.tpm,     mutBroc.features$high.tpm,     mutGFP.features$low.tpm,     mutGFP.features$high.tpm     )
logodds.all  <- c(simulation.features$low.logodds, simulation.features$low.med.logodds,  simulation.features$med.high.logodds, simulation.features$high.logodds, mutBroc.features$low.logodds, mutBroc.features$high.logodds, mutGFP.features$low.logodds, mutGFP.features$high.logodds )
tpm.breaks     <- seq(min(tpm.all    )-1, max(tpm.all    )+1, length.out=200)
logodds.breaks <- seq(min(logodds.all)-1, max(logodds.all)+1, length.out=40)
#
classes <- c("Neutral", "GoF", "LoF", "Uniform", "Other")
par(mfrow=c(length(classes),8))
for (i in 1:length(classes)){
   features2hists(simulation.features[simulation.features$class==classes[i],], tpm.breaks, logodds.breaks, cols[i], classes[i])
}
par(mfrow=c(length(classes),8))
features2hists(mutBroc.features, tpm.breaks, logodds.breaks, brocCol, 'Brocolli')
features2hists( mutGFP.features, tpm.breaks, logodds.breaks,  gfpCol,      'GFP')
######################################################################
#PCA: Principle Coordinate Analysis
class <- c(rep('Brocolli', length(mutBroc.features[,1])),rep('GFP', length(mutGFP.features[,1])))
ordination.data <- cbind(rbind(mutBroc.features, mutGFP.features),class)
ordination.data <- rbind(simulation.features,ordination.data)
simulation.pca <- prcomp(ordination.data[,1:8], scale = TRUE)
ordination.cols <- c(rep(c(rep(cols[1],N/length(classes)),
			   rep(cols[2],N/length(classes)),
			   rep(cols[3],N/length(classes)),
			   rep(cols[4],N/length(classes)),
			   rep(cols[5],N/length(classes))),reps),
			   rep(brocCol, length(mutBroc.features[,1])),
			   rep(gfpCol, length(mutGFP.features[,1])))
ordination.pch  <- c(rep('o',length(simulation.features[,1])),rep('+', length(mutBroc.features[,1])),rep('x', length(mutGFP.features[,1])))
par(mfrow=c(1,1),cex=2.0)
plot(simulation.pca$x[,1],simulation.pca$x[,2],col=ordination.cols,pch=ordination.pch,xlab='PCA1',ylab='PCA2',main='Principal component analysis',xlim=c(-5,7),ylim=c(-5,5))
legend('bottomright', c('Neutral', 'GoF', 'LoF', 'Uniform', 'Other','Brocolli','GFP'), col=c(cols,brocCol,gfpCol), pch=c('o','o','o','o','o','+','x'))
######################################################################
#TRAIN A RANDOM FOREST TO DISTINGUISH BETWEEN THE CLASSES (Neutral, GoF, LoF, Uniform or Other). 
ind <- sample(2,nrow(simulation.features),replace=TRUE,prob=c(0.8,0.2))
trainData <- simulation.features[ind==1,]
testData  <- simulation.features[ind==2,]
mutagenesis.rf <- randomForest(as.factor(class) ~ ., data=trainData, importance=TRUE, ntree=150, proximity=TRUE)
save(mutagenesis.rf,file = "trainClassifier2-simulation-rf.RData")
round(importance(mutagenesis.rf), 2)
###################################
#EVALUATE THE PERFORMANCE OF THE RF ON THE TESTDATA THAT WAS EXCLUDED FROM TRAINING:
mutagenesis.pred <- predict(mutagenesis.rf, testData)  
# 3 - Show Crosstab results:
(crossValTable <- table(observed = testData$class, predicted = mutagenesis.pred))
tp <- crossValTable[1,1] + crossValTable[2,2] + crossValTable[3,3]
tn <- crossValTable[4,4] ##OTHER
fp <- crossValTable[1,2] + crossValTable[1,3] + crossValTable[1,4] + crossValTable[2,3] + crossValTable[2,4]  + crossValTable[3,4] 
fn <- crossValTable[2,1] + crossValTable[3,1] + crossValTable[4,1] + crossValTable[3,2] + crossValTable[4,2]  + crossValTable[4,3] 
(sens <- tp/(tp+fn))
(ppv  <- tp/(tp+fp))
(spec <- tn/(tn+fp))
(acc  <- (tp+tn)/(tp+tn+fp+fn))
(f1   <- (2*sens*ppv)/(sens+ppv))
numer.1 <- (tp+fp)*(tp+fn)
numer.2 <- (tn+fp)*(tn+fn)
numer.3 <- (log10(numer.1) + log10(numer.2))/2 #Using logspace to get around a R integer overflow problem
(mcc  <- (tp*tn - fp*fn) / (10^numer.3) )
###################################
#PLOT RF PERFORMANCE
par(mfrow=c(2,1),mar=c(5.1, 4.1, 4.1, 4.1))
nms <- row.names(crossValTable)
res <- plot(matrix(crossValTable, nrow = length(classes), ncol = length(classes), dimnames = list(nms, nms)), cex=0.5, digits=0, xlab='Predicted', ylab='Category', xaxt = "n", yaxt = "n", col=rev(hcl.colors(12, palette = "Sunset")), main = 'Contingency Table: RF predictions vs sample category') #col=heat.colors
for (i in 1:nrow(crossValTable)) {
  for (j in 1:ncol(crossValTable)) {  
    args <- res$cell.text[[i,j]]
    args$cex    <- 1.5
    args$pos    <- 1
    do.call(text, args)
  }
}
par(las=1,mar=c(5.1, 8.1, 4.1, 2.1))
accuracies <- c(sens,ppv,spec,acc,f1,mcc)
acc.bp<-barplot(accuracies,  names.arg=c("Sensitivity","PPV","Specificity", "Accuracy","F1 score","MCC"), horiz=TRUE, main="Accuracy Measures on 10% of witheld simulated data", xlim=c(0,1) )
for (i in 1:length(accuracies)){
    text(0.8,acc.bp[i],round(accuracies[i], digits = 3),pos=4)
}
par(mfrow=c(1,1))
plot(mutagenesis.rf)
varImpPlot(mutagenesis.rf)
######################################################################
######################################################################
#RF PREDICTIONS ON SORTED SEQUENCING DATA: 	     
#BROCCOLI
(predictions.Broc       <- predict(mutagenesis.rf, mutBroc.features        )  )
(predictions.Broc.probs <- predict(mutagenesis.rf, mutBroc.features, "prob")  )
(outcomes.Broc <- table(predictions.Broc))
(gofBroc <- na.omit(cbind(mutBroc[predictions.Broc=='GoF',], predictions.Broc.probs[predictions.Broc == 'GoF',1])))
 lofBroc <- na.omit(cbind(mutBroc[predictions.Broc=='LoF',], predictions.Broc.probs[predictions.Broc == 'LoF',2]))
write.csv(gofBroc, file="gofBroc.Rdata",append=TRUE)
write.csv(lofBroc, file="lofBroc.Rdata",append=TRUE)
#GFP	     
(predictions.GFP        <- predict(mutagenesis.rf, mutGFP.features         )  )
(predictions.GFP.probs  <- predict(mutagenesis.rf, mutGFP.features,  "prob")  )
(outcomes.GFP <- table(predictions.GFP))
(gofGFP <- na.omit(cbind(mutGFP[predictions.GFP=='GoF',], predictions.GFP.probs[predictions.GFP=='GoF',1])))
 lofGFP <- na.omit(cbind(mutGFP[predictions.GFP=='LoF',], predictions.GFP.probs[predictions.GFP=='LoF',2]))
write.csv(gofGFP, file="gofGFP.Rdata",append=TRUE)
write.csv(lofGFP, file="lofGFP.Rdata",append=TRUE)
#######
#PLOT PREDICTION RESULTS:
par(mfrow=c(2,4), las=1, cex=0.75)
#Broccoli
#BARPLOT: PROPORTION OF OUTCOMES
(broc.bp<-barplot(outcomes.Broc, horiz=TRUE, main="Brocolli classifications", xlim=c(0,max(as.numeric(c(outcomes.Broc,outcomes.GFP)))),col=brocCol ))
for (i in 1:length(outcomes.Broc)){
    text(100,broc.bp[i],paste(round(100*outcomes.Broc[i]/sum(outcomes.Broc), digits = 1),"%, N=",outcomes.Broc[i],sep=""),pos=4)
}
##
if(length(gofBroc[,8]) > 0){
    #GoF #READS
    plot(NA,NA, ylim=c(0,max(gofBroc[,5:8], na.rm = TRUE) ),xlim=c(0.8,4),xlab="", ylab="# reads", main="Broccoli GoFs", xaxt = "n")
    axis(1, 4:1, c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofBroc[,8])){
        points(1:4,gofBroc[i,5:8], pch=20,col=cols[2],cex=3) 
    	lines( 1:4,gofBroc[i,5:8], pch=20,col=cols[2],lwd=3,lty=1) 
        points(1:4,gofBroc[i,5:8], pch=20,col='black',cex=1) 
    	lines( 1:4,gofBroc[i,5:8], pch=20,col='black',lwd=1,lty=1) 
    }
    #
    gofBrocFeatures<-mutBroc.features[predictions.Broc=='GoF',]
    #COLUMNS: low.tpm low.med.tpm med.high.tpm high.tpm low.logodds low.med.logodds med.high.logodds high.logodds
    #GoF TPM FEATURES
    plot(NA,NA, ylim=c(min(gofBrocFeatures[,1:4], na.rm = TRUE),max(gofBrocFeatures[,1:4], na.rm = TRUE) ),xlim=c(0.8,4),xlab="", ylab="TPM values", main="Broccoli GoFs", xaxt = "n")
    axis(1, 4:1, c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofBrocFeatures[,4])){
        points(1:4,gofBrocFeatures[i,1:4], pch=20,col=cols[2],cex=3) 
    	lines( 1:4,gofBrocFeatures[i,1:4], pch=20,col=cols[2],lwd=3,lty=1) 
        points(1:4,gofBrocFeatures[i,1:4], pch=20,col='black',cex=1) 
    	lines( 1:4,gofBrocFeatures[i,1:4], pch=20,col='black',lwd=1,lty=1) 
    }
    #GoF LOGODDS FEATURES log10(TPMs/PRESORT TPMs)
    plot(NA,NA, ylim=c(min(gofBrocFeatures[,5:8], na.rm = TRUE),max(gofBrocFeatures[,5:8], na.rm = TRUE) ),xlim=c(0.8,4),xlab="", ylab="log10(sort/presort)", main="Broccoli GoFs", xaxt = "n")
    axis(1, 4:1, c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofBrocFeatures[,5])){
        points(1:4,gofBrocFeatures[i,5:8], pch=20,col=cols[2],cex=3) 
    	lines( 1:4,gofBrocFeatures[i,5:8], pch=20,col=cols[2],lwd=3,lty=1) 
        points(1:4,gofBrocFeatures[i,5:8], pch=20,col='black',cex=1) 
    	lines( 1:4,gofBrocFeatures[i,5:8], pch=20,col='black',lwd=1,lty=1) 
    }
} else {
    plot(NA,NA)
    plot(NA,NA)
    plot(NA,NA)     
}


######
#GFP
(gfp.bp<-barplot(outcomes.GFP,  horiz=TRUE, main=     "GFP classifications", xlim=c(0,max(as.numeric(c(outcomes.Broc,outcomes.GFP)))),col=gfpCol ))
for (i in 1:length(gfp.bp)){
    text(100,gfp.bp[i],paste(round( 100*outcomes.GFP[i]/sum( outcomes.GFP), digits = 1),"%, N=",outcomes.GFP[i],sep=""),pos=4)
}
##
if(length(gofGFP[,8]) > 0){
    #GoF #READS
    plot(NA,NA, ylim=c(0,max(gofGFP[,5:8], na.rm = TRUE)),xlim=c(0.8,4),xlab="", ylab="# reads", main="GFP GoFs", xaxt = "n")
    axis(1, c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofGFP[,8])){
      points(1:4,gofGFP[i,5:8], pch=20,col=cols[2],cex=3) 
      lines( 1:4,gofGFP[i,5:8], pch=20,col=cols[2],lwd=3,lty=1)
      points(1:4,gofGFP[i,5:8], pch=20,col='black',cex=1) 
      lines( 1:4,gofGFP[i,5:8], pch=20,col='black',lwd=1,lty=1) 
      
    }
    #
    gofGFPFeatures<-mutGFP.features[predictions.GFP=='GoF',]
    #COLUMNS: low.tpm low.med.tpm med.high.tpm high.tpm low.logodds low.med.logodds med.high.logodds high.logodds
    #GoF TPM FEATURES
    plot(NA,NA, ylim=c(min(gofGFPFeatures[,1:4], na.rm = TRUE),max(gofGFPFeatures[,1:4], na.rm = TRUE) ),xlim=c(0.8,4),xlab="", ylab="TPM values", main="GFP GoFs", xaxt = "n")
    axis(1, 4:1, c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofGFPFeatures[,4])){
        points(1:4,gofGFPFeatures[i,1:4], pch=20,col=cols[2],cex=3) 
    	lines( 1:4,gofGFPFeatures[i,1:4], pch=20,col=cols[2],lwd=3,lty=1) 
        points(1:4,gofGFPFeatures[i,1:4], pch=20,col='black',cex=1) 
    	lines( 1:4,gofGFPFeatures[i,1:4], pch=20,col='black',lwd=1,lty=1) 
    }
    #GoF LOGODDS FEATURES log10(TPMs/PRESORT TPMs)
    plot(NA,NA, ylim=c(min(gofGFPFeatures[,5:8], na.rm = TRUE),max(gofGFPFeatures[,5:8], na.rm = TRUE) ),xlim=c(0.8,4),xlab="", ylab="log10(sort/presort)", main="GFP GoFs", xaxt = "n")
    axis(1, 4:1, c('High','High-Mid.', 'Mid.-Low','Low'))
    for (i in 1:length(gofGFPFeatures[,5])){
        points(1:4,gofGFPFeatures[i,5:8], pch=20,col=cols[2],cex=3) 
    	lines( 1:4,gofGFPFeatures[i,5:8], pch=20,col=cols[2],lwd=3,lty=1) 
        points(1:4,gofGFPFeatures[i,5:8], pch=20,col='black',cex=1) 
    	lines( 1:4,gofGFPFeatures[i,5:8], pch=20,col='black',lwd=1,lty=1) 
    }
} else {
    plot(NA,NA)
    plot(NA,NA)
    plot(NA,NA)     
}
######################################################################
#######
#A JITTERPLOT FUNCTION:
jitterPlot <- function(x, pos, col) {
    xp <-  rnorm(length(x), mean = pos, sd = 0.12)
    points(x,        xp,col=col, bg=col, pch=21)
    q  <-  quantile(x, probs<-(0:4)/4, na.rm=TRUE)
    ########!!!!!! https://www.r-bloggers.com/whisker-of-boxplot/
    upperWhisker <- min(max(x, na.rm=TRUE), q[4] + 1.5 * abs(q[4] - q[1]), na.rm=TRUE ) 
    lowerWhisker <- max(min(x, na.rm=TRUE), q[1] - 1.5 * abs(q[4] - q[1]), na.rm=TRUE )
    q  <-  c( lowerWhisker, q[2], median(x, na.rm=TRUE), q[4], upperWhisker)
    w  <-  c(          0.1, 0.25,      0.35, 0.25,          0.1)
    lwd <- c(            1,     2,        3,    2,            1)+1
    for(i in 1:length(q)){
        lines(c(q[i],q[i]),c(pos-w[i],pos+w[i]),lwd=lwd[i])
    }
#    lines(c(0.65,1.35),c(median(x),median(x)),lwd=3)
}
#######

par(cex=1.5,cex.axis=1.5,cex.lab=1.5,cex.main=2,las=1,mfrow=c(1,3),mar=c(6, 10, 3, 2)+0.1,mgp=c(3,2,0), bty='l')
#PROPORTION OF READS
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(0,1),xlab="", ylab="", main="", yaxt = "n")
par(mgp=c(3,1,0))
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('Proportion of reads', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
        cnts <- mutBroc[predictions.Broc==classes[j], 4+i]
    	tots <- mutBroc[predictions.Broc==classes[j],   9]
    	#pos <- length(classes)*(i-1)+j
	pos <- pos + 1
    	jitterPlot( cnts/tots , pos,  cols[j])
    }
    lines(c(0,1), c(pos+0.5, pos+0.5), col='lightgrey')    
}
legend('bottomright',classes, col=cols, fill=cols)
#TPMs
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(0,3000),xlab="", ylab="", main="Broccoli", yaxt = "n")
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('TPM', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
	pos <- pos + 1
    	jitterPlot( mutBroc.features[predictions.Broc==classes[j],i], pos,  cols[j]) #4*i-(j-1),  cols[j])
	#text( 3000, 4*i-(j-1), paste(i, j, classes[j], length(mutBroc.features[predictions.Broc==classes[j],j])), sep=":", cex=0.75)
	#text( 3000, 4*(i-1)+j, paste(i, j, classes[j], length(mutBroc.features[predictions.Broc==classes[j],j])), sep=":", cex=0.75)
    }
    lines(c(0,3000), c(pos+0.5, pos+0.5), col='lightgrey')
}
legend('bottomright',classes, col=cols, fill=cols)
#LOGODDS
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(-1.0,1.0),xlab="", ylab="", main="", yaxt = "n")
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('log10(sort/presort)', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
	pos <- pos + 1
    	jitterPlot( mutBroc.features[predictions.Broc==classes[j],4+i], pos ,  cols[j])
    }
    lines(c(-1,1), c(pos+0.5, pos+0.5), col='lightgrey')
}
legend('bottomright',classes, col=cols, fill=cols)


#######
par(cex=2.0,cex.axis=2,cex.lab=2,cex.main=2,las=1,mfrow=c(1,3),mar=c(6, 15, 3, 2)+0.1,mgp=c(3,2,0), bty='l')
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(0,1),xlab="", ylab="", main="", yaxt = "n")
par(mgp=c(3,1,0))
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('Proportion of reads', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
        cnts <- mutGFP[predictions.GFP==classes[j], 4+i]
    	tots <- mutGFP[predictions.GFP==classes[j],   9]
	pos <- pos + 1
    	jitterPlot( cnts/tots , pos,  cols[j])
    }
    lines(c(0,1), c(pos+0.5, pos+0.5), col='lightgrey')
}
legend('bottomright',classes, col=cols, fill=cols)
#TPMs
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(0,3000),xlab="", ylab="", main="GFP", yaxt = "n")
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('TPM', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
	pos <- pos + 1
    	jitterPlot( mutGFP.features[predictions.GFP==classes[j],i], pos,  cols[j])
    }
    lines(c(0,3000), c(pos+0.5, pos+0.5), col='lightgrey')
}
legend('bottomright',classes, col=cols, fill=cols)
#LOGODDS
plot(NA,NA, ylim=c(0,4*length(classes)+1),xlim=c(-1.0,1.0),xlab="", ylab="", main="", yaxt = "n")
axis(2, length(classes)*c(4,3,2,1), c('High','High-Mid.', 'Mid.-Low','Low'))
mtext('log10(sort/presort)', side = 1, line = 5, cex=3)
pos <- 0
for (i in 1:4 ){#Bin loop
    for (j in 1:length(classes)){#Classes loop
	pos <- pos + 1
	jitterPlot( mutGFP.features[predictions.GFP==classes[j],4+i], pos,  cols[j])
    }
    lines(c(-1,1), c(pos+0.5, pos+0.5), col='lightgrey')
}
legend('bottomright',classes, col=cols, fill=cols)

######################################################################
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
            
#Fold 1, 3x motif            
struct  <- '..(((.(((((((((((((((((((.((.(.((((((((....))))))..)).).)).)))...)))))))).....(((((((((((...)).)))))))))..))))).))).))).(((((((((((((((.((.(.((((((((....))))))..)).).)).)))...))))))..............))))))'
structA <- unlist(strsplit(struct,""))
maxP    <- max(predictions.Broc.probs[,1:2])
prev    <- 0
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1) #mar: bottom, left, top, right
plot(NA,NA, ylim=c(-1*maxP,maxP),xlim=c(8,202),xlab="", ylab="", main="Brocolli sequence", xaxt = "n", yaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
axis(2, at=seq(-0.6,0.6,by=0.2), labels=abs(seq(-0.6,0.6,by=0.2)))
lines(c(-10,202), c( 0.00, 0.00), col='lightgray')
lines(c(-10,202), c( 0.50, 0.50), col='lightgray')
lines(c(-10,202), c(-0.50,-0.50), col='lightgray')
brocClassesProbsPrint=c()
for (i in 1:length(mutBroc[,1])){    
    offSet <- nuc2num( mutBroc[i,3] )
    xPos   <- mutBroc[i,1] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    text( mutBroc[i,1],    0.00,            mutBroc[i,2] , pos=3, cex=0.60) #real sequence #pos: 1-below, 2-left, 3-above, 4-right
    text( mutBroc[i,1],    maxP,    structA[mutBroc[i,1]], pos=3, cex=0.60) #structure
    text( mutBroc[i,1],   -0.05,    structA[mutBroc[i,1]], pos=3, cex=0.60) #structure
    text( mutBroc[i,1], -1*maxP,    structA[mutBroc[i,1]], pos=3, cex=0.60) #structure
    if(predictions.Broc[i]=='GoF'){
        points(  xPos,               predictions.Broc.probs[i,1] , pch=20,col=cols[2])
    	lines( c(xPos,xPos), c(0,    predictions.Broc.probs[i,1]),        col=cols[2],lwd=3,lty=1)
	text(    xPos,               predictions.Broc.probs[i,1], mutBroc[i,3], pos=3, cex=0.60)
    } else if(predictions.Broc[i]=='LoF'){
        points(  xPos,            -1*predictions.Broc.probs[i,2] , pch=20,col=cols[3])
    	lines( c(xPos,xPos), c(0, -1*predictions.Broc.probs[i,2]),        col=cols[3],lwd=3,lty=1)
	text(    xPos,            -1*predictions.Broc.probs[i,2], mutBroc[i,3], pos=1, cex=0.60)
    }
    
    if ((mutBroc[i,1] %% 5)==1 && prev != mutBroc[i,1]){
       lines(c(mutBroc[i,1],mutBroc[i,1]), c( -1.00, 1.00), col='lightgray')
       prev    <- mutBroc[i,1]
    }
    brocClassesProbsPrint <- rbind(brocClassesProbsPrint, c(mutBroc[i,1], mutBroc[i,2], mutBroc[i,3], predictions.Broc[i],  predictions.Broc.probs[i,predictions.Broc[i]]))
}
write.csv(brocClassesProbsPrint, file="brocProbs.Rdata")
text(200, 0.9*maxP,'GoF',col=cols[2], font=2, cex=3)
text(200,-0.9*maxP,'LoF',col=cols[3], font=2, cex=3)
for (i in seq(0,200,by=10)){
    text( i, -0.02, i, pos=3, cex=0.60, col='black')
    points(i,0,pch='|',col='lightgray')
}

############

prot <- 'MRKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYG'
protA <- unlist(strsplit(prot,""))
protPos <- 1
prev    <- 0
maxP <- max(predictions.GFP.probs[,1:2])
par(las=0,bty='l',mfrow=c(1,1),mar=c(2, 7, 4, 1)+0.1) #mar: bottom, left, top, right
plot(NA,NA, ylim=c(-1*maxP,maxP),xlim=c(8,202),xlab="", ylab="", main="GFP sequence", xaxt = "n", yaxt = "n")
mtext("Classification probability", side = 2, line = 3, cex=2)
axis(2, at=seq(-0.6,0.6,by=0.2), labels=abs(seq(-0.6,0.6,by=0.2)))
lines(c(-10,202), c( 0.00, 0.00), col='lightgray')
lines(c(-10,202), c( 0.50, 0.50), col='lightgray')
lines(c(-10,202), c(-0.50,-0.50), col='lightgray')
gfpClassesProbsPrint=c()
for (i in 1:length(mutGFP[,1])){    
    offSet <- nuc2num( mutGFP[i,3] )
    xPos   <- mutGFP[i,1] + (offSet-1)/5 - 0.3   # each nt appears at i + (-0.3, -0.1, 0.1, 0.3)
    text( mutGFP[i,1],       0.00,         mutGFP[i,2] , pos=3, cex=0.60) #real sequence
    if( (mutGFP[i,1] %% 3)==1 && prev != mutGFP[i,1] ){
    	text( mutGFP[i,1],    maxP,    protA[protPos], pos=3, cex=0.60) #protein sequence
    	text( mutGFP[i,1],  -0.05,     protA[protPos], pos=3, cex=0.60) #protein sequence
    	text( mutGFP[i,1], -1*maxP,    protA[protPos], pos=3, cex=0.60) #protein sequence
       	lines(c(mutGFP[i,1],mutGFP[i,1]), c( -1.00, 1.00), col='lightgray')
	protPos <- protPos+1
	prev    <- mutGFP[i,1]
    }
    if(predictions.GFP[i]=='GoF'){
        points(  xPos,               predictions.GFP.probs[i,1] , pch=20,col=cols[2])
    	lines( c(xPos,xPos), c(0,    predictions.GFP.probs[i,1]),        col=cols[2],lwd=3,lty=1)
	text(    xPos,               predictions.GFP.probs[i,1], mutGFP[i,3], pos=3, cex=0.60)
    } else if(predictions.GFP[i]=='LoF'){
        points(  xPos,            -1*predictions.GFP.probs[i,2] , pch=20,col=cols[3])
    	lines( c(xPos,xPos), c(0, -1*predictions.GFP.probs[i,2]),        col=cols[3],lwd=3,lty=1)
	text(    xPos,            -1*predictions.GFP.probs[i,2], mutGFP[i,3], pos=1, cex=0.60)
    }
    gfpClassesProbsPrint <- rbind(gfpClassesProbsPrint, c(mutGFP[i,1], mutGFP[i,2], mutGFP[i,3], predictions.GFP[i],  predictions.GFP.probs[i,predictions.GFP[i]]))
}
write.csv(gfpClassesProbsPrint, file="gfpProbs.Rdata")
text(200, 0.9*maxP,'GoF',col=cols[2], font=2, cex=3)
text(200,-0.9*maxP,'LoF',col=cols[3], font=2, cex=3)
for (i in seq(0,200,by=10)){
    text( i, -0.02, i, pos=3, cex=0.60, col='black')
    points(i,0,pch='|',col='lightgray')
}

dev.off()
############
