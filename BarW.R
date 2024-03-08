# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DNABarcodes")
# 
# library(DNABarcodes)
# barcodes <- c("AGGT", "TTCC", "CTGA", "GCAA")
# barcode.set.distances(barcodes)
# barcode.set.distances(barcodes,metric="seqlev")
##########distancia robusta
##base---datos
##K--numero de grupos
library(spatstat.random)
library(TDA)
library(reshape2)
library(ggplot2)
rm(list=ls())
distanciaR<-function(base, K){
  base<-as.data.frame(base)
  ##numero de datos
  NN<-nrow(base)
  x1<-rep(1:K,1000)[1:NN]  
  base$etiquetas<-sample(x1,length(x1))  
  persitencias<-list()
  for(ii in 1:K){
    baseK<-subset(base, etiquetas==ii)[,1:2]
    
    persitencias[[ii]]<-ripsDiag(baseK, maxdimension, DiagLim, printProgress = FALSE) }
  
  MM<-matrix(NA, ncol=K,nrow=K)
  
  for(i in 1:K){for(j in 1:K){
    MM[i,j]<-wasserstein(persitencias[[i]][["diagram"]], persitencias[[j]][["diagram"]], p = 1,
                         dimension = 1)
  }}
  
  distanciass<-c()
  for(it in 1:K){
    distanciass[it]<-sum(sort(MM[it,])[2:(1+round(K/2))]) / round(K/2)
  }
  persitencias[[which.min(distanciass)]]
}
################################################################
################################################################
wassersteinDist12<-c()
wassersteinDist13<-c()
wassersteinDist1R2<-c()
wassersteinDist1R3<-c()


for(jj in 1:100){ print(jj)
N<-600
epsil<-0.1
n<-round(N*epsil)
N1<-N-n

######Datos sin ruido
XX1 <- circleUnif(N)
#######datos con ruido espilon
XX2<-XX1+matrix(rnorm(2*N,0,0.05), ncol=2)
######


grupos<-3
X <- rMatClust(grupos, scale=0.25, mu=round(n/grupos), win =owin(c(-0.5,0.5), c(-0.5,0.5)))

ruido<-cbind(X$x,X$y)
######datos con ruido Matern
XX3<-rbind(XX2[1:N1,],ruido)


#plot(XX3)
#lines(XX2, type="p", col="red")

DiagLim <- 5
maxdimension <- 1

#######diagrama de persisitencia para los originales 
Diag1 <- ripsDiag(XX1, maxdimension, DiagLim, printProgress = FALSE)

#######diagrama de persisitencia para los originales con ruido epsilon
Diag2 <- ripsDiag(XX2, maxdimension, DiagLim, printProgress = FALSE)

#######diagrama de persisitencia para los originales con ruido Matern
Diag3 <- ripsDiag(XX3, maxdimension, DiagLim, printProgress = FALSE) 
K=6
DiagR2<-distanciaR(XX2,K)
DiagR3<-distanciaR(XX3,K)





wassersteinDist12[jj] <- wasserstein(Diag1[["diagram"]], Diag2[["diagram"]], p = 1,
                               dimension = 1)

wassersteinDist13[jj]  <- wasserstein(Diag1[["diagram"]], Diag3[["diagram"]], p = 1,
                               dimension = 1)
wassersteinDist1R2[jj]<-wasserstein(Diag1[["diagram"]], DiagR2[["diagram"]], p = 1,
                                    dimension = 1)

wassersteinDist1R3[jj]<-wasserstein(Diag1[["diagram"]], DiagR3[["diagram"]], p = 1,
                                    dimension = 1)

  
  

}

setwd("~/Dropbox/2023_Ale_Emilien_Leo/Tema1-MedianofMeans/Simulaciones en R")
pdf("diag1.pdf",height = 4,width = 4 )
plot(Diag1[["diagram"]],add=F) 
dev.off()
pdf("diag2.pdf",height = 4,width = 4 )
plot(Diag2[["diagram"]],add=F) 
dev.off()
pdf("diag3.pdf" ,height = 4,width = 4 )
plot(Diag3[["diagram"]],add=F) 
dev.off()

datt<-data.frame(XX1)
names(datt)<-c("x","y")
pdf("muestra1.pdf",height = 4,width = 4)
ggplot(datt, aes(x=x, y=y ))+geom_point(alpha=0.3, color="#31a354") +theme_bw()
dev.off()


datt<-data.frame(XX2)
names(datt)<-c("x","y")
pdf("muestra2.pdf",height = 4,width = 4)
ggplot(datt, aes(x=x, y=y ))+geom_point(alpha=0.3, color="#fdae6b") +theme_bw()
dev.off()

datt<-data.frame(XX3)
datt$out<- as.factor(c(rep(1,N1),rep(0,nrow(ruido)) ))
names(datt)<-c("x","y", "out")
pdf("muestra3.pdf",height = 4,width = 4)
ggplot(datt, aes(x=x, y=y, color= out ))+geom_point(alpha=0.3) +theme_bw()+ theme(legend.position="none")+scale_color_manual(values=c("#de2d26","#fdae6b"))
dev.off()




ASAS1<-as.data.frame(cbind(`Scenario 1`= wassersteinDist12, `Scenario 2`=wassersteinDist13))
ASAS1$Method<-"Wasserstein"
ASASR<-as.data.frame(cbind(`Scenario 1`= wassersteinDist1R2, `Scenario 2`=wassersteinDist1R3))
ASASR$Method<-"Robust_Wasserstein"

ASAS<-rbind(ASAS1,ASASR)

library(reshape2)
library(ggplot2)
basef<-melt(ASAS)
names(basef)<-c( "Method", "Persistence_Diagrams", "Distance to Dgm")
basef$Persistence_Diagrams<-as.factor(basef$Persistence_Diagrams)
basef$Method<-as.factor(basef$Method)

library(RColorBrewer)
pdf("boxdiag.pdf",height = 4,width = 6)
ggplot(basef, aes(x=Persistence_Diagrams, y=`Distance to Dgm`, color=Method ))+geom_boxplot()+theme_bw() +scale_color_manual(values=c("#a8ddb5","#43a2ca"))
dev.off()

save.image("persistencia.Rdata")
