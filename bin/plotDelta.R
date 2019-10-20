#!/usr/bin/Rscript

files<-dir(path="./")

delta1 <- read.table("data/delta-delta/computeDeltaDelta-1mut.out", header=T, sep="\t")
delta4 <- read.table("data/delta-delta/computeDeltaDelta-4mut.dat", header=T, sep="\t")

cor.test(delta1$dBS[ delta1$type == 'rna'    ], delta1$ddG[ delta1$type == 'rna'     ], method="spearman")
cor.test(delta1$dBS[ delta1$type == 'protein'], delta1$ddG[ delta1$type == 'protein' ], method="spearman")

cor.test(delta4$dBS[ delta4$type == 'rna'    ], delta4$ddG[ delta4$type == 'rna'     ], method="spearman")
cor.test(delta4$dBS[ delta4$type == 'protein'], delta4$ddG[ delta4$type == 'protein' ], method="spearman")

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

bitscoreAnnotation <- function(){
    arrows(0,3,3.0,3.0, lwd=2)
    text(  6.5,3,        "convergent", pos=3, srt=270,cex=2.5)
    w<-0.1 
    lines(c(0,0),c(3-w,3+w),lwd=3)
    arrows(0,3.0,-28,3.0,lwd=3)
    text(  -31,3.0,      "divergent", pos=3, srt=90,cex=2.5)
}

deltaGAnnotation <- function(){
    arrows(0,3.0,11,3.0,lwd=2)
    text(  13, 3.0,      "destabilising", pos=3, srt=270,cex=2.5)
    w<-0.1 
    lines(c(0,0),c(3-w,3+w),lwd=3)
    arrows(0,3.0,-4.5,3.0,lwd=3)
    text(  -6.0,3.0,      "stabilising", pos=3, srt=90,cex=2.5)
}


cairo_pdf(file="manuscript/figures/figure2ab.pdf", width=24, height=12, family="DejaVu Sans")  ##Switched to cairo -- normal function failed to plot Delta!?
par(cex=2.5,cex.axis=3,cex.lab=3,cex.main=3,las=1,mfrow=c(1,2),mar=c(6, 11, 3, 2)+0.1,mgp=c(3,2,0), bty='l') #c(bottom, left, top, right)
p.col <- adjustcolor("lightskyblue1",  alpha.f = 0.4)
r.col <- adjustcolor("pink",           alpha.f = 0.4)
######
plot(NA,NA, ylim=c(0.5,5.5),xlim=c(-7,15),xlab="", ylab="", main="", yaxt = "n")
lines(c(0,0),c(0,6),lwd=7, col=adjustcolor("grey", alpha.f = 0.7))
par(mgp=c(3,1,0))
axis(2, c(1,2,4,5), c('RNA','Protein', 'RNA','Protein'))
mtext(expression(paste(Delta, Delta, "G (kcal/mol)")), side = 1, line = 5, cex=3)
jitterPlot(delta1$ddG[ delta1$type == 'rna'    ], 4,  r.col)
###
jitterPlot(delta1$ddG[ delta1$type == 'protein'], 5,  p.col)
###
jitterPlot(delta4$ddG[ delta4$type == 'rna'    ], 1,  r.col)
###
jitterPlot(delta4$ddG[ delta4$type == 'protein'], 2,  p.col)
###
deltaGAnnotation()
text(6,5.5,"One mutation" ,font=2, cex=3)
text(6,2.5,"Four mutations",font=2, cex=3)
######
par(mgp=c(3,2,0))
plot(NA,NA, ylim=c(0.5,5.5),xlim=c(-34,8.5),xlab="", ylab="", main="", yaxt = "n")
lines(c(0,0),c(0,6),lwd=7, col=adjustcolor("grey", alpha.f = 0.7))
par(mgp=c(3,1,0))
axis(2, c(1,2,4,5), c('RNA','Protein', 'RNA','Protein'))
#axis(1, seq(-30, 10, 10), seq(-30, 10, 10), line=0)
mtext(expression(paste(Delta, "bitscore")), side = 1, line = 5, cex=3)
jitterPlot(delta1$dBS[ delta1$type == 'rna'    ], 4,  r.col)
###
jitterPlot(delta1$dBS[ delta1$type == 'protein'], 5,  p.col)
###
jitterPlot(delta4$dBS[ delta4$type == 'rna'    ], 1,  r.col)
###
jitterPlot(delta4$dBS[ delta4$type == 'protein'], 2,  p.col)
###
bitscoreAnnotation()
text(-13,5.5,"One mutation" ,font=2, cex=3)
text(-13,2.5,"Four mutations",font=2, cex=3)
###############
dev.off()

#convert  -flatten -density 200 -trim  manuscript/figures/figure2ab.pdf -quality  100   manuscript/figures/figure2ab.png

######################################################################
#Proportions of divergent or destabilising mutations:

##ADD THE NA COUNTS TOO!

#1 mutation
100*length(delta1$dBS[ delta1$type == 'rna'    & (delta1$dBS < 0 | delta1$ddG >  0) ] ) / length( delta1$dBS[ delta1$type == 'rna'    ] )
100*(length(delta1$dBS[ delta1$type == 'protein'& (delta1$dBS < 0 | delta1$ddG >  0) ] ) + length(delta1$ddG[ delta1$type == 'protein' & delta1$ddG == 'NA']) )  / length( delta1$dBS[ delta1$type == 'protein'] )
#66% RNA, 60% protein

#4 mutations
100*length(delta4$dBS[ delta4$type == 'rna'    & (delta4$dBS < 0 | delta4$ddG >  0) ] ) / length( delta4$dBS[ delta4$type == 'rna'    ] )
100*length(delta4$dBS[ delta4$type == 'protein'& (delta4$dBS < 0 | delta4$ddG >  0) ] ) / length( delta4$dBS[ delta4$type == 'protein'] )
#98% RNA, 93% protein

#
#Proportions of divergent mutations:  
#1 mutation
100*length(delta1$dBS[ delta1$type == 'rna'    & (delta1$dBS < 0) ] ) / length( delta1$dBS[ delta1$type == 'rna'    ] )
100*length(delta1$dBS[ delta1$type == 'protein'& (delta1$dBS < 0) ] ) / length( delta1$dBS[ delta1$type == 'protein'] )
#61% RNA, 50% protein

#4 mutations
100*length(delta4$dBS[ delta4$type == 'rna'    & (delta4$dBS < 0) ] ) / length( delta4$dBS[ delta4$type == 'rna'    ] )
100*length(delta4$dBS[ delta4$type == 'protein'& (delta4$dBS < 0) ] ) / length( delta4$dBS[ delta4$type == 'protein'] )
#96% RNA, 89% protein

#Proportions of destabilising mutations:  
#1 mutation
100*length(delta1$dBS[ delta1$type == 'rna'    & (delta1$ddG >  0) ] ) / length( delta1$dBS[ delta1$type == 'rna'    ] )
100*length(delta1$dBS[ delta1$type == 'protein'& (delta1$ddG >  0) ] ) / length( delta1$dBS[ delta1$type == 'protein'] )
#41% RNA, 30% protein

#4 mutations
100*length(delta4$dBS[ delta4$type == 'rna'    & (delta4$ddG >  0) ] ) / length( delta4$dBS[ delta4$type == 'rna'    ] )
100*length(delta4$dBS[ delta4$type == 'protein'& (delta4$ddG >  0) ] ) / length( delta4$dBS[ delta4$type == 'protein'] )
#73% RNA, 58% protein


## library(gplots)

## par(cex=1.5,cex.axis=2,cex.lab=2,mfrow=c(2,2),mar=c(7, 5, 2, 1)+0.1)
## input  <-list(
##     Divergent=    delta1$ID[ delta1$type == 'rna'     & (delta1$dBS < 0) ],
##     Destabilising=delta1$ID[ delta1$type == 'rna'     & (delta1$ddG > 0) ])
## #venn(input, names=c("Divergent", "Destabilising"))
## plot.new()
## vp <- venn.diagram( input, fill = 2:3, alpha = 0.3, filename = NULL); grid.draw(vp)
## input  <-list(
##     Divergent=    delta1$ID[ delta1$type == 'protein' & (delta1$dBS < 0) ],
##     Destabilising=delta1$ID[ delta1$type == 'protein' & (delta1$ddG > 0) ])
## #venn(input, names=c("Divergent", "Destabilising"))
## plot.new()
## vp <- venn.diagram( input, fill = 2:3, alpha = 0.3, filename = NULL); grid.draw(vp)

## input  <-list(
##     Divergent=    delta4$ID[ delta4$type == 'rna'    & (delta4$dBS < 0) ],
##     Destabilising=delta4$ID[ delta4$type == 'rna'    & (delta4$ddG > 0) ])
## #venn(input, names=c("Divergent", "Destabilising"))
## plot.new()
## vp <- venn.diagram( input, fill = 2:3, alpha = 0.3, filename = NULL); grid.draw(vp)

## input  <-list(
##     Divergent=    delta4$ID[ delta4$type == 'protein' & (delta4$dBS < 0) ],
##     Destabilising=delta4$ID[ delta4$type == 'protein' & (delta4$ddG > 0) ])
## #venn(input, names=c("Divergent", "Destabilising"))
## plot.new()
## vp <- venn.diagram( input, fill = 2:3, alpha = 0.3, filename = NULL); grid.draw(vp)


######################################################################
#Supplementary figures:

cairo_pdf(file="manuscript/figures/suppfigure3ab-delta.pdf", width=24, height=12, family="DejaVu Sans")
par(cex=1.5,cex.axis=2,cex.lab=2,mfrow=c(1,2),mar=c(7, 5, 2, 1)+0.1)
plot(NA,NA, main="1 point mutn", ylab=expression(paste(Delta, Delta, "G (kcal/mol)")), xlab=expression(paste(Delta, "bitscore (HMM/CM)")), xlim=c(-45,6), ylim=c(-10,20)  )  #xlim=c(-10,5), ylim=c(-3,5)  )  
lines(c(0,0),c(-100,100), lwd=1)
lines(c(-100,100),c(0,0), lwd=1)
points(delta1$dBS[ delta1$type == 'rna'     ], delta1$ddG[ delta1$type == 'rna'     ], col="deeppink4",   pch=20, cex=0.5)
DF <- na.omit(  data.frame(x = delta1$dBS[ delta1$type == 'rna'     ], y = delta1$ddG[ delta1$type == 'rna'     ])  )
lines(lowess( DF, f=.2), col="deeppink4",lwd=4)
points(delta1$dBS[ delta1$type == 'protein' ], delta1$ddG[ delta1$type == 'protein' ], col="dodgerblue4", pch=20, cex=0.5)
DF <- na.omit(  data.frame(x = delta1$dBS[ delta1$type == 'protein'     ], y = delta1$ddG[ delta1$type == 'protein'     ])  )
lines(lowess( DF, f=.2), col="dodgerblue4",lwd=4)
###
plot(NA,NA, main="4 point mutns", ylab=expression(paste(Delta, Delta, "G (kcal/mol)")), xlab=expression(paste(Delta, "bitscore (HMM/CM)")), xlim=c(-45,6), ylim=c(-10,20)  )  #xlim=c(-10,5), ylim=c(-3,5)  )  
lines(c(0,0),c(-100,100), lwd=1)
lines(c(-100,100),c(0,0), lwd=1)
points(delta4$dBS[ delta4$type == 'rna'     ], delta4$ddG[ delta4$type == 'rna'     ], col="deeppink4",   pch=20, cex=0.5)
DF <- na.omit(  data.frame(x = delta4$dBS[ delta4$type == 'rna'     ], y = delta4$ddG[ delta4$type == 'rna'     ])  )
lines(lowess( DF, f=.2), col="deeppink4",lwd=4)
points(delta4$dBS[ delta4$type == 'protein' ], delta4$ddG[ delta4$type == 'protein' ], col="dodgerblue4", pch=20, cex=0.5)
DF <- na.omit(  data.frame(x = delta4$dBS[ delta4$type == 'protein'     ], y = delta4$ddG[ delta4$type == 'protein'     ])  )
lines(lowess( DF, f=.2), col="dodgerblue4",lwd=4)
dev.off()

#convert  -flatten -density 200 -trim  manuscript/figures/suppfigure3ab-delta.pdf -quality  100   manuscript/figures/suppfigure3ab-delta.png

cairo_pdf(file="manuscript/figures/suppfigure3cdef-delta.pdf", width=12, height=12, family="DejaVu Sans")
par(cex=1.5,cex.axis=2,cex.lab=2, las=1,mfrow=c(2,2),mar=c(7, 5, 2, 1)+0.1, mgp=c(3,1,0))
col <- c("orchid1", "skyblue", "lightpink", "lightblue1" )
mn<-min(delta1$dBS, na.rm=TRUE)
mx<-max(delta1$dBS, na.rm=TRUE)
breaks<-seq(mn, mx, length.out=750)
h.r<-hist(delta1$dBS[ delta1$type == 'rna'    ], breaks=breaks, plot=F)
h.p<-hist(delta1$dBS[ delta1$type == 'protein'], breaks=breaks, plot=F)
plot( h.p$mids, log10(h.p$counts+1), lwd=5, col=col[2],xlim=c(-70,10),ylim=c(0,3.5),main="1 point mutn",xlab=expression(paste(Delta, "bitscore (HMM/CM)")), ylab="Frequency", type='l', yaxt = "n")
axis(2, 0:3, c(0,10^(1:3)))
lines(h.r$mids, log10(h.r$counts+1), lwd=5, col=col[1])
legend("topleft", c("RNA", "Protein"), fill=c(col[1], col[2]), cex=2.0)
###
mn<-min(delta1$ddG, na.rm=TRUE)
mx<-max(delta1$ddG, na.rm=TRUE)
breaks<-seq(mn, mx, length.out=30)
h.r<-hist(delta1$ddG[ delta1$type == 'rna'    ], breaks=breaks, plot=F)
h.p<-hist(delta1$ddG[ delta1$type == 'protein'], breaks=breaks, plot=F)
plot( h.p$mids, log10(h.p$counts+1), lwd=5, col=col[2],xlim=c(-10,20),ylim=c(0,3.5),main="1 point mutn",xlab=expression(paste(Delta, Delta, "G (kcal/mol)")), ylab="Frequency", type='l', yaxt = "n")
axis(2, 0:3, c(0,10^(1:3)))
lines(h.r$mids, log10(h.r$counts+1), lwd=5, col=col[1])
#####
mn<-min(delta4$dBS, na.rm=TRUE)
mx<-max(delta4$dBS, na.rm=TRUE)
breaks<-seq(mn, mx, length.out=750)
h.r<-hist(delta4$dBS[ delta4$type == 'rna'    ], breaks=breaks, plot=F)
h.p<-hist(delta4$dBS[ delta4$type == 'protein'], breaks=breaks, plot=F)
plot( h.p$mids, log10(h.p$counts+1), lwd=5, col=col[2],xlim=c(-70,10),ylim=c(0,3.5),main="4 point mutns",xlab=expression(paste(Delta, "bitscore (HMM/CM)")), ylab="Frequency", type='l', yaxt = "n")
axis(2, 0:3, c(0,10^(1:3)))
lines(h.r$mids, log10(h.r$counts+1), lwd=5, col=col[1])
###
mn<-min(delta4$ddG, na.rm=TRUE)
mx<-max(delta4$ddG, na.rm=TRUE)
breaks<-seq(mn, mx, length.out=30)
h.r<-hist(delta4$ddG[ delta4$type == 'rna'    ], breaks=breaks, plot=F)
h.p<-hist(delta4$ddG[ delta4$type == 'protein'], breaks=breaks, plot=F)
plot( h.p$mids, log10(h.p$counts+1), lwd=5, col=col[2],xlim=c(-10,20),ylim=c(0,3.5),main="4 point mutns",xlab=expression(paste(Delta, Delta, "G (kcal/mol)")), ylab="Frequency", type='l', yaxt = "n")
axis(2, 0:3, c(0,10^(1:3)))
lines(h.r$mids, log10(h.r$counts+1), lwd=5, col=col[1])
dev.off()

#convert  -flatten -density 200 -trim  manuscript/figures/suppfigure3cdef-delta.pdf -quality  100   manuscript/figures/suppfigure3cdef-delta.png







