#!/usr/bin/Rscript
#cd data/sgr-structural/
files<-dir(path="./")

mRNA.length <- 102
sRNA.length <- 227



######################################################################
#SUBSTITUTIONS
mrates<-c(0.001, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
if( length( files[grep("sgrT-mRNA-mrate", files)] ) == length( files[grep("sgrS-sRNA-mrate", files)] ) ){
    p.files<-files[grep("sgrT-mRNA-mrate", files)]
    r.files<-files[grep("sgrS-sRNA-mrate", files)]
    len <- length(p.files)
    substits              <- matrix(0, nrow = 100, ncol = 2*len)
    substitsP             <- matrix(0, nrow = 100, ncol =   len)
    substitsR             <- matrix(0, nrow = 100, ncol =   len)
    substitsP.num.substits<- matrix(0, nrow = 100, ncol =   len)
    substitsR.num.substits<- matrix(0, nrow = 100, ncol =   len)
    sgrT.len              <- matrix(0, nrow = 100, ncol =   len)
    substitsP.seqs        <- matrix("", nrow = 100, ncol =   len)
    substitsR.seqs        <- matrix("", nrow = 100, ncol =   len)
    substitsP.truncated   <- matrix("", nrow = 100, ncol =   len)
    substitsR.truncated   <- matrix("", nrow = 100, ncol =   len)
    for(i in 1:len){
    	  sgrT              <- read.table(p.files[i], sep="\t")$V2
    	  sgrS              <- read.table(r.files[i], sep="\t")$V2
	  sgrT.num.substits <- read.table(p.files[i], sep="\t")$V4
	  sgrS.num.substits <- read.table(r.files[i], sep="\t")$V4
    	  sgrT.seqs         <- read.table(p.files[i], sep="\t")$V10
    	  sgrS.seqs         <- read.table(r.files[i], sep="\t")$V8
	  
	  substits[ 1:length(sgrT),            2*i-1] <- sgrT
	  substits[ 1:length(sgrS),            2*i  ] <- sgrS
	  substitsP[1:length(sgrT),              i  ] <- sgrT
	  substitsR[1:length(sgrS),              i  ] <- sgrS
	  substitsP.num.substits[1:length(sgrT), i  ] <- sgrT.num.substits
	  substitsR.num.substits[1:length(sgrS), i  ] <- sgrS.num.substits
	  substitsP.seqs[        1:length(sgrT), i  ] <- sgrT.seqs
	  substitsR.seqs[        1:length(sgrS), i  ] <- sgrS.seqs
	  sgrT.len[ 1:length(sgrT),              i  ] <- nchar(as.character( read.table(p.files[i], sep="\t")$V10 ))

	  substitsP.truncated[1:length(sgrT),   i  ] <- nchar(                 as.vector(read.table(p.files[i], sep="\t")$V10) )   < ( 0.75 * mRNA.length/3 )
	  substitsR.truncated[1:length(sgrS),   i  ] <- nchar( gsub( "-", "",  as.vector(read.table(r.files[i], sep="\t")$V8 ) ) ) < ( 0.75 * sRNA.length   )

	  }
}

xp.subs <- 1000*c(substitsP.num.substits)/mRNA.length
xr.subs <- 1000*substitsR.num.substits/sRNA.length  
substitsP <- c(substitsP)     #[x$ix]
substitsR <- c(substitsR)     #[x$ix]

########
#INDELS
indelrates<-c(0.0001, 0.001, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
colnames.indels <- c()
if( length( files[grep("sgrT-mRNA-indelrate", files)] ) == length( files[grep("sgrS-sRNA-indelrate", files)] ) ){
    p.files<-files[grep("sgrT-mRNA-indelrate", files)]
    r.files<-files[grep("sgrS-sRNA-indelrate", files)]
    len <- length(p.files)
    indels            <- matrix(0, nrow = 100, ncol = 2*len)
    indelsP           <- matrix(0, nrow = 100, ncol =   len)
    indelsR           <- matrix(0, nrow = 100, ncol =   len)
    indelsP.num.indels<- matrix(0, nrow = 100, ncol =   len)
    indelsR.num.indels<- matrix(0, nrow = 100, ncol =   len)
    indelsP.truncated <- matrix(0, nrow = 100, ncol =   len)
    indelsR.truncated <- matrix(0, nrow = 100, ncol =   len)
    for(i in 1:len){
    	  sgrT <- read.table(p.files[i], sep="\t")$V2
    	  sgrS <- read.table(r.files[i], sep="\t")$V2
	  f1<-read.table(p.files[i], sep="\t")$V5
	  f2<-read.table(p.files[i], sep="\t")$V6
	  sgrT.num.indels <- as.numeric(f1) + as.numeric(f2)
	  f3<-read.table(r.files[i], sep="\t")$V5
	  f4<-read.table(r.files[i], sep="\t")$V6
	  sgrS.num.indels <-  as.numeric(f3) + as.numeric(f4)
	  indels[ 1:length(sgrT),            2*i-1] <- sgrT
	  indels[ 1:length(sgrS),            2*i  ] <- sgrS
	  indelsP[1:length(sgrT),              i  ] <- sgrT
	  indelsR[1:length(sgrS),              i  ] <- sgrS
	  indelsP.num.indels[1:length(sgrT),   i  ] <- sgrT.num.indels
	  indelsR.num.indels[1:length(sgrS),   i  ] <- sgrS.num.indels

	  indelsP.truncated[1:length(sgrT),   i  ] <- nchar(                 as.vector(read.table(p.files[i], sep="\t")$V10) )   < ( 0.75 * mRNA.length/3 )
	  indelsR.truncated[1:length(sgrS),   i  ] <- nchar( gsub( "-", "",  as.vector(read.table(r.files[i], sep="\t")$V8 ) ) ) < ( 0.75 * sRNA.length   )

	  colnames.indels <- c(colnames.indels, paste('sgrT.', indelrates[i], sep=""), paste('sgrS.', indelrates[i], sep=""))
    }
}
colnames(indels) <- colnames.indels
rownames.indels<-c(); for(i in 1:100){ rownames.indels <- c(rownames.indels, paste("sample", i, sep="")) }
rownames(indels) <- rownames.indels
write.table(indels, file="plotStructureMutageneratorData.tsv", sep="\t")

xp.in <- 1000*c(indelsP.num.indels)/mRNA.length
xr.in <- 1000*c(indelsR.num.indels)/sRNA.length
indelsP <- c(indelsP)
indelsR <- c(indelsR)

scatter.indels <- cbind(xp.in, indelsP, xr.in, indelsR)

write.table(indels, file="plotStructureMutageneratorData.tsv", sep="\t")



###################################
p.bootcol <- adjustcolor("lightskyblue1",  alpha.f = 0.01)
r.bootcol <- adjustcolor("pink",           alpha.f = 0.01)

cairo_pdf(file="../../manuscript/figures/figure2cd.pdf", width=24, height=12, family="DejaVu Sans")
#par(cex=2.5)
par(cex=2.5,cex.axis=3,cex.lab=3,cex.main=3,las=0,mfrow=c(1,2),mar=c(6, 11, 3, 2)+0.1,mgp=c(3,2,0), bty='l') #c(bottom, left, top, right)
#SUBSITUTIONS
plot(NA, NA, , col="dodgerblue4", pch=20, ylab="", xlab="", xlim=c(0,500), cex=0.5, ylim=c(-1,1) )
mtext("substitutions/kb", side = 1, line = 5, cex=3)
mtext("Structure Rho",    side = 2, line = 5, cex=3)
par(mgp=c(3,1,0))
minmax <- min(max(xp.subs), max(xr.subs))
xx<-seq(0,minmax,length.out=2000)
for(i in 1:500){
      ind<-sample.int(length(xp.subs), size = length(xp.subs), replace = TRUE)
      lines(xx, predict(loess( (substitsP[ind]) ~ c(xp.subs)[ind], span=0.5, degree=1), xx), col = p.bootcol, lwd=4)
      ind<-sample.int(length(c(xr.subs)), size = length(c(xr.subs)), replace = TRUE)
      lines(xx, predict(loess( (c(substitsR)[ind]) ~ c(xr.subs)[ind], span=0.5, degree=1), xx), col = r.bootcol, lwd=4)
}
bbp <- predict(loess(  substitsP  ~   xp.subs,  span=0.5, degree=1), xx)
bbr <- predict(loess(c(substitsR) ~ c(xr.subs), span=0.5, degree=1), xx)
points(xp.subs[substitsP.truncated == FALSE], substitsP[substitsP.truncated == FALSE], col="dodgerblue4", pch=20, cex=1.5)
points(xr.subs[substitsR.truncated == FALSE], substitsR[substitsR.truncated == FALSE], col="deeppink4",   pch=20, cex=1.5)
points(xp.subs[substitsP.truncated == TRUE ], substitsP[substitsP.truncated == TRUE ], col="dodgerblue4", pch=17, cex=1.5)
points(xr.subs[substitsR.truncated == TRUE ], substitsR[substitsR.truncated == TRUE ], col="deeppink4",   pch=17, cex=1.5)
lines(xx, bbp, col = "dodgerblue4", lwd=4)
lines(xx, bbr, col = "deeppink4",   lwd=4)
legend("bottomleft", c("Protein (sgrT)", "RNA (sgrS)"), col=c("dodgerblue4", "deeppink4"), fil=c("dodgerblue4", "deeppink4"),cex=2.5,ncol=2)
#INDELS##############
par(mgp=c(3,2,0))
plot(NA,NA, col="dodgerblue4", pch=20, ylab="", xlab="", xlim=c(0,500), cex=0.5, ylim=c(-1,1)  )  
mtext("INDELs/kb",     side = 1, line = 5, cex=3)
mtext("Structure Rho", side = 2, line = 5, cex=3)
par(mgp=c(3,1,0))
minmax <- min(max(xp.in), max(xr.in))
xx<-seq(0,minmax,length.out=2000)
for(i in 1:500){
      ind<-sample.int(length(xp.in), size = length(xp.in), replace = TRUE)
      lines(xx, predict(loess(as.vector(indelsP[ind]) ~ as.vector(xp.in[ind]), span=0.5, degree=1), xx), col = p.bootcol, lwd=4)
      ind<-sample.int(length(xr.in), size = length(xr.in), replace = TRUE)
      lines(xx, predict(loess( (indelsR[ind]) ~ xr.in[ind], span=0.5, degree=1), xx), col = r.bootcol, lwd=4)
}
bbp <- predict(loess(as.vector(indelsP) ~ as.vector(xp.in), span=0.5, degree=1), xx)
bbr <- predict(loess(indelsR ~ xr.in, span=0.5, degree=1), xx)
points(xp.in[indelsP.truncated == FALSE], indelsP[indelsP.truncated == FALSE], col="dodgerblue4", pch=20, cex=1.5)
points(xr.in[indelsR.truncated == FALSE], indelsR[indelsR.truncated == FALSE], col="deeppink4",   pch=20, cex=1.5)
points(xp.in[indelsP.truncated == TRUE ], indelsP[indelsP.truncated == TRUE ], col="dodgerblue4", pch=17, cex=1.5)
points(xr.in[indelsR.truncated == TRUE ], indelsR[indelsR.truncated == TRUE ], col="deeppink4",   pch=17, cex=1.5)
lines(xx, bbp, col = "dodgerblue4", lwd=4)
lines(xx, bbr, col = "deeppink4", lwd=4)
legend("bottomleft", c("Protein (sgrT)", "RNA (sgrS)"), col=c("dodgerblue4", "deeppink4"), fil=c("dodgerblue4", "deeppink4"),cex=2.5,ncol=2)
dev.off()

#convert  -flatten -density 200 -trim  manuscript/figures/figure2cd.pdf -quality  100   manuscript/figures/figure2cd.png

######################################################################
