
library(lemon)
library(patchwork)
library(scoring)

quantile.95<-function(data.tmp,row.names.tmp)
{
N.ROW.Q<-length(row.names.tmp)
mat.quant.tmp<-as.data.frame(matrix(0,N.ROW.Q,4))
names(mat.quant.tmp)<-c("Mean","Q2.5","Q50","Q97.5")

for (i.N.ROW.Q in 1:N.ROW.Q)
mat.quant.tmp[i.N.ROW.Q,]<-c(mean(data.tmp[,i.N.ROW.Q]),quantile(data.tmp[,i.N.ROW.Q],probs=c(0.025,0.5,0.975)))

row.names(mat.quant.tmp)<-row.names.tmp
mat.quant.tmp

}

file.to.numeric<-function(file.t)
{
N.c<-dim(file.t)[2]
for (j in 1:N.c)
{
file.t[,j]<-as.numeric(paste(file.t[,j]))
}
file.t
}



gen.df.freq<-function(df.t=fit.tmp)
{
fit.tmp.txt<-rep(0,length(df.t[,1]))
N.r<-dim(df.t)[1]
N.c<-dim(df.t)[2]
print(N.r)
print(N.c)
fit.tmp.tt<-df.t
for (j in 1:N.c)
fit.tmp.tt[,j]<-as.numeric(paste(df.t[,j]))


for (i in 1:N.r)
{
fit.tmp.txt[i]<-paste(as.numeric(fit.tmp.tt[i,]),collapse="")
}

df.tmp<-as.data.frame(table(fit.tmp.txt))
df.tmp<-df.tmp[sort(df.tmp$Freq,index.return=T,decreasing=T)$ix,]
df.tmp$p<-df.tmp$Freq/sum(df.tmp$Freq)
df.tmp
}



GMSIM.to.DF<-function(data.sim.tmp,name.vars.ord)
{
data.sim.tmp<-as.data.frame(data.sim.tmp)
data.sim.tmp<-data.sim.tmp[,name.vars.ord]
N.cols.sim<-dim(data.sim.tmp)[2]
for (i.ncols.sim in 1:N.cols.sim)
{
data.sim.tmp[,i.ncols.sim]<-as.numeric(data.sim.tmp[,i.ncols.sim])-1
}
data.sim.tmp
}

gen.sim.GM<-function(data.tmp.base,n.sim.tmp=10000,names.vars.t=c("hpv6","hpv18","hpv56","hpv11","hpv16","hpv51","infection"))
{
var.data.f<-names(data.tmp.base)

hpv.tab.f <- xtabs(~ ., data = data.tmp.base)
hpv.mod.f <- dmod(~ .^., data = hpv.tab.f, marginal = var.data.f)
# select the model (AIC)
hpv.sel.f <- stepwise(hpv.mod.f,details=1)

par(mfrow=c(1,2))
plot(hpv.mod.f)
plot(hpv.sel.f)
title("SIMULATED")

hpv.gc.f <- hpv.sel.f$glist
str(hpv.gc.f)
hpv.ug.f <- ugList(hpv.gc.f)
potlist.f <- extractPOT(data.tmp.base, hpv.ug.f, smooth=0.01)
# Compile Net
Comp.hpv.f<-grain(compilePOT(potlist.f))
hpv.sim.f<-simulate(Comp.hpv.f,nsim=n.sim.tmp)
N.col.sim<-dim(hpv.sim.f)[2]
if (N.col.sim<length(var.data.f))
{
print("ENTRA SAT")
hpv.sel.f <- stepwise(hpv.mod.f,details=1,k=0)
hpv.gc.f <- hpv.sel.f$glist
str(hpv.gc.f)
hpv.ug.f <- ugList(hpv.gc.f)
potlist.f <- extractPOT(data.tmp.base, hpv.ug.f, smooth=0.01)
# Compile Net
Comp.hpv.f<-grain(compilePOT(potlist.f))
hpv.sim.f<-simulate(Comp.hpv.f,nsim=n.sim.tmp)
}
hpv.sim.f<-GMSIM.to.DF(hpv.sim.f,names.vars.t)
hpv.sim.f
}

brierSc<-function(observed.tmp,predicted.tmp)
{
data.o.p<-as.data.frame(cbind(observed.tmp,predicted.tmp))
names(data.o.p)<-c("obs","pred")
data.o.p$dif2<-sqrt((data.o.p$obs-data.o.p$pred)^2)
#mean(brierscore(obs ~ pred, data=data.o.p))
mean(data.o.p$dif2)
}

Dist.vec<-function(observed.tmp,predicted.tmp)
{
#Euclidean Distance
data.o.p<-as.data.frame(cbind(observed.tmp*100000,predicted.tmp*100000))
names(data.o.p)<-c("obs","pred")
data.o.p$dif2<-sqrt((data.o.p$obs-data.o.p$pred)^2)
mean(data.o.p$dif2)
}

Chev.vec<-function(observed.tmp,predicted.tmp)
{
#Chebyshev distance
data.o.p<-round(as.data.frame(cbind(observed.tmp*100000,predicted.tmp*100000)))
names(data.o.p)<-c("obs","pred")
data.o.p$difchev<-abs(data.o.p$obs-data.o.p$pred)
max(data.o.p$difchev)
}

#############################################
#############################################
#############################################
#############################################
#############################################
#install.packages("remotes")
#remotes::install_github("RomeroBarata/bimba")

#install.packages("installr")
#install.packages("gRim")
#install.packages("Rgraphviz")
#install.packages("gRain")
#install.packages("scoring")
#install.packages("bimba")

library(installr)
library(gRim)
library(Rgraphviz)
library(gRain)
library(scoring)
library(bimba)


dir.in<-"/PATH/"
file.in<-"hpvsimTP.txt"
fit.in.p<-paste(dir.in,file.in,sep="")
hpv.sim<-as.data.frame(read.table(fit.in.p,header=T))
summary(hpv.sim)
names.data.sel<-names(hpv.sim)
str(hpv.sim)
########################################################################################
########################################################################################
##### STARTS SIMULATION
########################################################################################
########################################################################################

}

N.sim.MAT<-500
Cond.N<-c(50,100,250)


mat.Dist.vec<-as.data.frame(matrix(0,length(Cond.N),8))
names(mat.Dist.vec)<-c("GM","BDS1","BDS2","SM","ADA","MW","ROS","SLS")
mat.Dist.vec.R<-mat.Dist.vec
mat.Chev.vec<-as.data.frame(matrix(0,length(Cond.N),8))
names(mat.Chev.vec)<-c("GM","BDS1","BDS2","SM","ADA","MW","ROS","SLS")
mat.Chev.vec.R<-mat.Chev.vec

for (i.cond in 1:length(Cond.N))

{


N.names<-length(names.data.sel)


sample.size<-Cond.N[i.cond]

for (i.sim.DB in 1:N.sim.MAT)
{
hpv.sim.sam<-hpv.sim[sample(nrow(hpv.sim),1000),]
runif.N<-runif(0.1,0.4,100)
sample.N<-sample(1:sample.size,size=sample.size,replace=T)

rand.sam<-runif(1,45,55)

hpv.sim.sam.GM<-try(gen.sim.GM(hpv.sim.sam[sample.N,],n.sim.tmp=sample.size,names.vars.t=names.data.sel),silent=TRUE)
if (!is.data.frame(hpv.sim.sam.GM)) 
{
print("ERROR GM")
hpv.sim.sam.GM<-hpv.sim.sam[1,]
}

hpv.sim.sam.BDS1<-try(as.data.frame(BDLSMOTE(hpv.sim.sam[sample.N,],perc_min=rand.sam,k=5,borderline=1),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.BDS1)) 
{
print("ERROR BDS1")
hpv.sim.sam.BDS1<-hpv.sim.sam[1,]
}


hpv.sim.sam.BDS2<-try(as.data.frame(BDLSMOTE(hpv.sim.sam[sample.N,],perc_min=rand.sam,k=5,borderline=2),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.BDS2)) 
{
print("ERROR BDS2")
hpv.sim.sam.BDS2<-hpv.sim.sam[1,]
}


hpv.sim.sam.SM<-try(as.data.frame(SMOTE(hpv.sim.sam[sample.N,],perc_min=rand.sam,k=5),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.SM)) 
{
print("ERROR SM")
hpv.sim.sam.SM<-hpv.sim.sam[1,]
}

hpv.sim.sam.ADA<-try(as.data.frame(ADASYN(hpv.sim.sam[sample.N,],perc_min=rand.sam,k=5),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.ADA)) 
{
print("ERROR ADA")
hpv.sim.sam.ADA<-hpv.sim.sam[1,]
}

hpv.sim.sam.MW<-try(as.data.frame(MWMOTE(hpv.sim.sam[sample.N,],perc_min=rand.sam),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.MW)) 
{
print("ERROR MW")
hpv.sim.sam.MW<-hpv.sim.sam[1,]
}


hpv.sim.sam.ROS<-try(as.data.frame(ROS(hpv.sim.sam[sample.N,],perc_min=rand.sam),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.ROS))
{
print("ERROR ROS")
hpv.sim.sam.ROS<-hpv.sim.sam[1,]
}

hpv.sim.sam.SLS<-try(as.data.frame(SLSMOTE(hpv.sim.sam[sample.N,],perc_min=rand.sam,k=5),silent=TRUE))
if (!is.data.frame(hpv.sim.sam.SLS)) 
{
print("ERROR SLS")
hpv.sim.sam.SLS<-hpv.sim.sam[1,]
}





if (i.sim.DB==1)
{
hpv.GM<-hpv.sim.sam.GM
hpv.BDS1<-hpv.sim.sam.BDS1
hpv.BDS2<-hpv.sim.sam.BDS2
hpv.SM<-hpv.sim.sam.SM
hpv.ADA<-hpv.sim.sam.ADA
hpv.MW<-hpv.sim.sam.MW
hpv.ROS<-hpv.sim.sam.ROS
hpv.SLS<-hpv.sim.sam.SLS
}

if (i.sim.DB>1)
{
hpv.GM<-rbind(hpv.GM,hpv.sim.sam.GM)
hpv.BDS1<-rbind(hpv.BDS1,hpv.sim.sam.BDS1)
hpv.BDS2<-rbind(hpv.BDS2,hpv.sim.sam.BDS2)
hpv.SM<-rbind(hpv.SM,hpv.sim.sam.SM)
hpv.ADA<-rbind(hpv.ADA,hpv.sim.sam.ADA)
hpv.MW<-rbind(hpv.MW,hpv.sim.sam.MW)
hpv.ROS<-rbind(hpv.ROS,hpv.sim.sam.ROS)
hpv.SLS<-rbind(hpv.SLS,hpv.sim.sam.SLS)
}

print("CONDITION")
print(i.cond)
print("SIMULATION")
print(i.sim.DB)
rm(hpv.sim.sam)

}

#Norm.vec
mat.Dist.vec[i.cond,1]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.GM,2,mean)))
mat.Dist.vec[i.cond,2]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.BDS1,2,mean)))
mat.Dist.vec[i.cond,3]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.BDS2,2,mean)))
mat.Dist.vec[i.cond,4]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.SM,2,mean)))
mat.Dist.vec[i.cond,5]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.ADA,2,mean)))
mat.Dist.vec[i.cond,6]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.MW,2,mean)))
mat.Dist.vec[i.cond,7]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.ROS,2,mean)))
mat.Dist.vec[i.cond,8]<-mean(Dist.vec(apply(hpv.sim,2,mean),apply(hpv.SLS,2,mean)))
mat.Dist.vec.R[i.cond,]<-rank(mat.Dist.vec[i.cond,])


mat.Chev.vec[i.cond,1]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.GM,2,mean)))
mat.Chev.vec[i.cond,2]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.BDS1,2,mean)))
mat.Chev.vec[i.cond,3]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.BDS2,2,mean)))
mat.Chev.vec[i.cond,4]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.SM,2,mean)))
mat.Chev.vec[i.cond,5]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.ADA,2,mean)))
mat.Chev.vec[i.cond,6]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.MW,2,mean)))
mat.Chev.vec[i.cond,7]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.ROS,2,mean)))
mat.Chev.vec[i.cond,8]<-mean(Chev.vec(apply(hpv.sim,2,mean),apply(hpv.SLS,2,mean)))
mat.Chev.vec.R[i.cond,]<-rank(mat.Chev.vec[i.cond,])

#rm(hpv.GM)
#rm(hpv.BDS1)
#rm(hpv.BDS2)
#rm(hpv.SM)
#rm(hpv.ADA)
#rm(hpv.MW)
#rm(hpv.ROS)
#rm(hpv.SLS)

}

