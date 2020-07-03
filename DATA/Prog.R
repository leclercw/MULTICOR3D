### Cas depl

rm(list=ls())
don1 <- read.table(file="datadepl.txt",head=TRUE,sep="")

iter=don1$ITE
depl=don1$U

png(file="Evol_depl.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
plot(depl~iter,pch=-1,col=0,lwd=2,xlim=c(min(iter),max(iter)),ylim=c(min(depl),max(depl)),xlab="Iterations",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)

mtext(expression(paste("Displacement")), side=2, line=5.1,cex=1.8) 

points(depl~iter, col= 1,pch=1, lwd=1.5)
lines(depl~iter, col= 1,lty=1,lwd=2.5)

legend("bottomright",legend=c("Dis.") 
,col=c(1),lty=c(1),pch=c(1),lwd=2.5,cex=1.8,inset=0.01, border = "black",  merge = TRUE, bg = "gray92")

dev.off()

### Cas dila

rm(list=ls())
don1 <- read.table(file="datadila.txt",head=TRUE,sep="")

iter=don1$ITE
deplx=don1$DEPX
deply=don1$DEPY
deplz=don1$DEPZ

png(file="Evol_dila.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
plot(deplx~iter,pch=-1,col=0,lwd=2,xlim=c(min(iter),max(iter)),ylim=c(min(deplx),max(deplx)),xlab="Iterations",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)

mtext(expression(paste("Displacement")), side=2, line=5.1,cex=1.8) 

lines(deplx~iter, col= 1,lty=1,lwd=2.5)
lines(deply~iter, col= 2,lty=2,lwd=2.5)
lines(deplz~iter, col= 3,lty=3,lwd=2.5)


legend("bottomright",legend=c("Dis. x","Dis. y","Dis. z") 
,col=c(1,2,3),lty=c(1,2,3),pch=c(-1,-1,-1),lwd=2.5,cex=1.8,inset=0.01, border = "black",  merge = TRUE, bg = "gray92")

dev.off()

### Cas nrj

rm(list=ls())
don1 <- read.table(file="dataenergy.txt",head=TRUE,sep="")

iter=don1$ITE
enerp=don1$EP
enerc=don1$EC

png(file="Evol_nrj.png", width = 600, height = 480)
par(mar=c(5.0, 8.0, 0.5, 7.0)) 
plot(enerp~iter,pch=0,col=0,lwd=0,lty=1,xlim=c(min(iter),max(iter)),ylim=c(0.,max(enerp)),xlab="Iterations",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)

mtext("Potential energy (J)", side=2, line=5.6,cex=1.8) 
lines(enerp~iter, col= 1,lty=1, lwd=2.5)

par(new=T,yaxs="r",xpd = NA) 

plot(enerc~iter,pch=0,col=0,lwd=0,lty=2,xlim=c(min(iter),max(iter)),ylim=c(0.,max(enerp)),xlab="",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
axis(side = 4, las = 1,cex.axis=1.8)
mtext("Elastic energy (J)", side=4, line=5.6,cex=1.8) 

lines(enerc~iter, col= 2,lty=2, lwd=2.5)

legend("topleft",legend=c("Potential energy","Elastic energy") 
,col=c(1,2),lty=c(1,2),pch=c(-1,-1),lwd=2.5,cex=1.6,inset=0.01, border = "black",  merge = TRUE, bg = "gray92")

dev.off()


### Cas reac

rm(list=ls())
don1 <- read.table(file="datareac.txt",head=TRUE,sep="")

iter=don1$ITE
force=abs(don1$F)

png(file="Evol_force.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
plot(force~iter,pch=-1,col=0,lwd=2,xlim=c(min(iter),max(iter)),ylim=c(min(force),max(force)),xlab="Iterations",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)

mtext(expression(paste("Load (N)")), side=2, line=5.1,cex=1.8) 

points(force~iter, col= 1,pch=1, lwd=1.5)
lines(force~iter, col= 1,lty=1,lwd=2.5)

legend("bottomright",legend=c("Load") 
,col=c(1),lty=c(1),pch=c(-1),lwd=2.5,cex=1.8,inset=0.01, border = "black",  merge = TRUE, bg = "gray92")

dev.off()

### Cas rupt

rm(list=ls())
don1 <- read.table(file="datarupt.txt",head=TRUE,sep="")

iter=don1$ITE
nrc=don1$NRC
nrcis=don1$NRCIS
nrt=don1$NRT
nrtot=don1$NRTOT

png(file="Evol_rupt.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
plot(nrtot~iter,pch=-1,col=0,lwd=2,xlim=c(min(iter),max(iter)),ylim=c(min(nrtot),max(nrtot)),xlab="Iterations",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)

mtext(expression(paste("Number of delete elements")), side=2, line=5.1,cex=1.8) 

lines(nrtot~iter, col= 1,lty=1,lwd=2.5)
lines(nrt~iter, col= 2,lty=2,lwd=2.5)
lines(nrc~iter, col= 3,lty=3,lwd=2.5)
lines(nrcis~iter, col= 4,lty=4,lwd=2.5)

legend("topleft",legend=c("Total","Tensile","Comp.","Shear") 
,col=c(1,2,3,4),lty=c(1,2,3,4),pch=c(-1,-1,-1,-1),lwd=2.5,cex=1.8,inset=0.01, border = "black",  merge = TRUE, bg = "gray92")

dev.off()


