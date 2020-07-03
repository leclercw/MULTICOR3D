### Cas def

rm(list=ls())
don1 <- read.table(file="datadef.txt",head=TRUE,sep="")

num=don1$NUM
def=don1$DEF

png(file="Evol_def.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
hist(def,col="green",main="",xlab="Strain",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F, freq=FALSE)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)
mtext(expression(paste("Density")), side=2, line=5.1,cex=1.8) 

densite <- density(def) # estimer la densité que représente ces différentes valeurs
lines(densite, col = "red",lwd=3) # Superposer une ligne de densité à l'histogramme

dev.off()

cove=sd(def)/mean(def)

### Cas sig

don2 <- read.table(file="datasig.txt",head=TRUE,sep="")

num=don2$NUM
sig=don2$SIG

png(file="Evol_sig.png")
par(mar=c(4.5, 7.2, 0.5, 0.5)) 
hist(sig,col="green",main="",xlab="Stress",ylab="",cex.lab=1.8,cex.axis=1.8,axes=F, freq=FALSE)
box()
axis(side = 1, cex.axis=1.8)
axis(side = 2, las = 1,cex.axis=1.8)
mtext(expression(paste("Density")), side=2, line=5.1,cex=1.8) 

densite <- density(sig) # estimer la densité que représente ces différentes valeurs
lines(densite, col = "red",lwd=3) # Superposer une ligne de densité à l'histogramme

dev.off()

covs=sd(sig)/mean(sig)

