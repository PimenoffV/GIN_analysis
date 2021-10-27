# R code for GIN simulation and validation with knonw oversampling methods.

rm(list = ls())
#install.packages(c("lemon", "patchwork", "scoring"))
#install.packages("lemon")
#install.packages("patchwork")
#install.packages("scoring")
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


dir.in<-"PATH to the FOLDER"
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


##########################################################################################
################## FIGURE
#########################################################################################

table10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
table20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5","#9c9ede", "#7375b5", "#4a5584", "#cedb9c", "#b5cf6b", "#8ca252","#7b4173")
manyEyes <- c("#9c9ede", "#7375b5", "#4a5584", "#cedb9c", "#b5cf6b", "#8ca252", "#637939", "#e7cb94", "#e7ba52", "#bd9e39", "#8c6d31", "#e7969c", "#d6616b", "#ad494a", "#843c39", "#de9ed6", "#ce6dbd", "#a55194", "#7b4173")
table4 <- c("wheat1", "tan1","salmon1","orangered1")

dir.in<-"PATH to the FOLDER"
#install.packages("ggplot2")
library(ggplot2)
 
# create a dataset
# c("GIN","BDS1","BDS2","SM","ADA","MW","ROS","SLS")
method <- c(rep("1.GIN" , 3) , rep("3.BDS1" , 3) , rep("4.BDS2" , 3) , rep("6.SM" , 3), rep("2.ADASYN" , 3), rep("8.MW" , 3), rep("5.ROS" , 3), rep("7.SLS" , 3))
Scheme <- rep(c("SC1(N=050)" , "SC2(N=100)" , "SC3(N=250)") , 8)
distance <- c(mat.Dist.vec.R$GM, mat.Dist.vec.R$BDS1, mat.Dist.vec.R$BDS2, mat.Dist.vec.R$SM, mat.Dist.vec.R$ADA, mat.Dist.vec.R$MW, mat.Dist.vec.R$ROS, mat.Dist.vec.R$SLS)
Euclidean <- c(mat.Dist.vec.R$GM, mat.Dist.vec.R$BDS1, mat.Dist.vec.R$BDS2, mat.Dist.vec.R$SM, mat.Dist.vec.R$ADA, mat.Dist.vec.R$MW, mat.Dist.vec.R$ROS, mat.Dist.vec.R$SLS)
measure <- c("Euclidean")
data1 <- data.frame(method,Scheme,distance,measure)
 
# Grouped

dist.euc<-ggplot(data1, aes(fill=Scheme, y=distance, x=method)) + 
    geom_bar(position="dodge", stat="identity")+ylab("Ranking Euclidean") + coord_flip() +
  scale_fill_manual(values = manyEyes) +
  facet_grid(.~Scheme)
dist.euc
library(ggplot2)
 
# create a dataset
# c("GM","BDS1","BDS2","SM","ADA","MW","ROS","SLS")
method <- c(rep("1.GIN" , 3) , rep("3.BDS1" , 3) , rep("4.BDS2" , 3) , rep("5.SM" , 3), rep("2.ADASYN" , 3), rep("6.MW" , 3), rep("6.ROS" , 3), rep("8.SLS" , 3))
Scheme <- rep(c("SC1(N=050)" , "SC2(N=100)" , "SC3(N=250)") , 8)
ranking <- c(mat.Chev.vec.R$GM, mat.Chev.vec.R$BDS1, mat.Chev.vec.R$BDS2, mat.Chev.vec.R$SM, mat.Chev.vec.R$ADA, mat.Chev.vec.R$MW, mat.Chev.vec.R$ROS, mat.Chev.vec.R$SLS)
Chebyshev <- c(mat.Chev.vec.R$GM, mat.Chev.vec.R$BDS1, mat.Chev.vec.R$BDS2, mat.Chev.vec.R$SM, mat.Chev.vec.R$ADA, mat.Chev.vec.R$MW, mat.Chev.vec.R$ROS, mat.Chev.vec.R$SLS)
measure <- c("Chebyshev")
data2 <- data.frame(method,Scheme,distance, measure)
 
# Grouped
dist.cheb<-ggplot(data2, aes(fill=Scheme, y=ranking, x=method)) + 
    geom_bar(position="dodge", stat="identity")+ylab("Ranking Chebyshev")  + coord_flip() +
  scale_fill_manual(values = manyEyes) +
    facet_wrap(vars(Scheme))
dist.cheb

# Both separate together
#############################
data <- data.frame(method,Scheme,Euclidean, Chebyshev)
str(data)
distance <- c("Euclidean","Chebyshev")

pdf(paste(dir.in,"Figure_4AB_individual.pdf", sep = ""), width = 6, height = 6)
for(v in distance){
  d=data[,c("method","Scheme",v)]
  p <- ggplot(d, aes(fill=Scheme, y=get(v), x=method)) + 
    geom_bar(position="dodge", width = 0.7, stat="identity") +
    ylab(v) + scale_fill_manual(values = manyEyes) +
  facet_wrap(vars(Scheme)) + coord_flip() 
  print(p)
  }
dev.off()

# Both plotted parallel
#############################
#install.packages("gridExtra")
#install.packages("patchwork")
#install.packages("arrangeGrob")
library("ggplot2")
library("gridExtra")
library(grid)
library(patchwork)
p1 <- ggplot(data1, aes(fill=Scheme, y=distance, x=method)) + 
  geom_bar(position="dodge", width = 0.7, stat="identity", colour="black")+ylab("Ranking Euclidean") + coord_flip() +
  scale_fill_manual(values = manyEyes) + labs(tag = "A") +
  theme(text = element_text(size=14))
  #facet_grid(.~Scheme)

p2 <- ggplot(data2, aes(fill=Scheme, y=ranking, x=method)) + 
  geom_bar(position="dodge", width = 0.7, stat="identity", colour="black")+ylab("Ranking Chebyshev")  + coord_flip() +
  scale_fill_manual(values = manyEyes) + labs(tag = "B") +
  theme(text = element_text(size=14))
  #facet_wrap(vars(Scheme))

pdf(paste(dir.in,"Figure_4AB_patch.parallel.pdf", sep = ""), width = 12, height = 6)
comb <- p1 + p2 & theme(legend.position = "right")
comb + plot_layout(guides = "collect") 
dev.off()

postscript(paste(dir.in,"Figure_4AB_patch.parallel.ps", sep = ""), width = 12, height = 6)
comb <- p1 + p2 & theme(legend.position = "right")
comb + plot_layout(guides = "collect") 
dev.off()

# face_grip version
#install.packages("lemon")

library(lemon)
pdf(paste(dir.in,"Figure_4AB_grid.parallel.pdf", sep = ""), width = 15, height = 8)
#p1A <- arrangeGrob(p1, top = textGrob("A", x = unit(0, "npc")
#                                               , y   = unit(1, "npc"), just=c("left","top"),
#                                               gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))

#p2B <- arrangeGrob(p2, top = textGrob("B", x = unit(0, "npc")
#                                               , y = unit(1, "npc"), just=c("left","top"),
#                                               gp=gpar(col="black", fontsize=18, fontfamily="Helvetica")))
grid.arrange(p1, p2, nrow = 1) +
grid_arrange_shared_legend(p1, p2, position = "bottom")
dev.off()

# face_wrap version
#############################
data0 <- rbind(data1, data2)
str(data0)
data0 <- data.frame(data0,Df = rep(c("data1","data2"),
                                     times=c(nrow(data1),nrow(data2))))

p <-  ggplot(data0, aes(fill=Scheme, y=distance, x=method)) + 
  geom_bar(position="dodge", stat="identity")+ylab(measure)  + coord_flip() +
  scale_fill_manual(values = manyEyes) +
  facet_wrap(~Df)
  #facet_grid(measure ~.)
p

install.packages("tidyverse")
install.packages("dplyr")
library(ggplot2)
library(dplyr)
library(tidyverse)

p <- ggplot(mpg, aes(displ, hwy)) + geom_point()
# Use vars() to supply faceting variables:
p + facet_wrap(vars(class))
# Control the number of rows and columns with nrow and ncol
p + facet_wrap(vars(class), nrow = 4)
# You can facet by multiple variables
ggplot(mpg, aes(displ, hwy)) +
  geom_point() +
  facet_wrap(vars(cyl, drv))



save.image("Sim.500.Samples.RData")
