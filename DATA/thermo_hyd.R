#### Récupération de la table
donnees <- read.table(file='result_thyd',head=TRUE)
attach(donnees)

tps=donnees$Temps
x=donnees$X*100
conw=donnees$Cmilieu
temp=donnees$Tmilieu

nd=nrow(donnees)

grad_conw=c((conw[2]-conw[1])/(tps[2]-tps[1]))
for (i in 2:nd-1){
grad=(conw[i+1]-conw[i-1])/(tps[i+1]-tps[i-1])	
grad
grad_conw=c(grad_conw,grad)
}
grad_conw=c(grad_conw,(conw[nd]-conw[nd-1])/(tps[nd]-tps[nd-1]))
grad_conw2=abs(grad_conw)

grad_conw=c()
for (i in 1:nd){
grad=0
num=0
for (j in 1:8){
if(i+j<=nd){
grad=grad+grad_conw2[i+j]
num=num+1
}
if(i-j>=1){
grad=grad+grad_conw2[i-j]
num=num+1
}
}

grad=grad/num
grad_conw=c(grad_conw,grad)
}

tps=tps/3600


png(file="X_t.png")
par(mar=c(4.8,4.8,0.5,0.5))

plot(tps,x,type="l",pch=0,lwd=2,xlim=c(0,max(tps)),ylim=c(min(x),max(x)),xlab="Temps (h)",ylab="X (% bs)",cex.axis=1.5,cex.lab=1.5)
box()
dev.off()

png(file="C_t.png")
par(mar=c(4.8,4.8,0.5,0.5))

plot(tps,conw,type="l",pch=0,lwd=2,xlim=c(0,max(tps)),ylim=c(min(conw),max(conw)),xlab="Temps (h)",ylab=expression(paste("C (mol.m"^-3,")")),cex.axis=1.5,cex.lab=1.5)
box()
dev.off()

png(file="gradC_t.png")
par(mar=c(4.8,4.8,0.5,0.5))

plot(tps,grad_conw,type="l",pch=0,lwd=2,xlim=c(0,max(tps)),ylim=c(min(grad_conw),max(grad_conw)),xlab="Temps (h)",ylab=expression(paste("Taux de séchage (mol.s"^-1,".m"^-3,")")),cex.axis=1.5,cex.lab=1.5)
box()
dev.off()

png(file="T_t.png")
par(mar=c(4.8,4.8,0.5,0.5))

plot(tps,temp,type="l",pch=0,lwd=2,xlim=c(0,max(tps)),ylim=c(min(temp),max(temp)),xlab="Temps (h)",ylab="T (°C)",cex.axis=1.5,cex.lab=1.5)
box()
dev.off()
