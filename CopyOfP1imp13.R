library(Matrix)
library(tidyverse)
library(MASS)
library(MVA)
library(maps)
library(MCMCpack)
library(doBy)
library(utilities)
library(irlba)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(vegan3d)
library(car)
library(igraph)
library(Rdimtools)
library(combinat)
library(gridExtra)

#loading the dataset
load("/cloud/project/df-forAranyak.RData")

#number of obs
n<-nrow(df)

#naming the columns
X1<-df$X.1
X2<-df$X.2
X3<-df$X.3
X4<-df$X.4
X5<-df$X.5
X6<-df$X.6
y<-df$dist
w<-df$claw


#matrix of ASE estimates
X_hat<-cbind(df$X.1,df$X.2,df$X.3,df$X.4,df$X.5,df$X.6)


claw<-as.character(df$claw)


#number of dimensions
d<-6

#neighbourhood parameter of localization graph
lambda<-0.50


#matrix of interpoint distances of 6d obs
B0<-as.matrix(dist(X_hat,method = "euclidean",
                   diag=TRUE,upper=TRUE))

#matrix of localization graph
BB<-ifelse(B0<lambda,B0,0)


colnames(BB)<-as.character(seq(1,nrow(B0),1))


#creating localization graph from adjacency matrix
g<-graph_from_adjacency_matrix(BB,
                               mode="undirected",
                               weighted=TRUE,
                               diag=FALSE,
                               add.colnames = NULL,
                               add.rownames = NA)





#matrix of shortest path distances
D<-shortest.paths(g,v=V(g),to=V(g))
Ds<-as.matrix(D)



#raw-stress minimization
MM<-mds(Ds,ndim = 1,type = "ratio",
        weightmat = NULL,
        init = "torgerson")


#raw-stress embeddings
z_hat<-as.vector(MM$conf)
print(z_hat)


#regression model of response and embeddings
linmod1<-lm(y~z_hat)
print(summary(linmod1))


#plot of response and embeddings with 
# fitted regression line
p1<-ggplot(df, aes(x=z_hat, y=y)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$")) 
p1

ggsave(file="P1plot13.eps", plot = p1,
       width = 4, height = 2,
       units = "in", dpi = 300)


#X1vsX2plot, sized by response y, for INTRO
YX12_hat<-cbind(X1,X2,y)
dfyx<-as.data.frame(YX12_hat)
p2<-ggplot(data = dfyx, aes(x=X1,y=X2,size=y)) +
  geom_point(alpha=0.5) +
  ylab(TeX("embedding dimension $2$"))+
  xlab(TeX("embedding dimension $1$"))+
  theme(legend.position = "none")
p2

g<-grid.arrange(p2,p1,ncol=2)


ggsave(file="P1plot13A.eps", plot = g,
       width = 5, height = 2.1,
       units = "in", dpi = 500,
       device = cairo_ps)









#plot as a whole, with regression line
# without separation by no. of claws
p2<-ggplot(df, aes(x=z_hat, y=y)) +
     geom_point() +
     geom_smooth(method = "lm", se=FALSE) +
     ylab(TeX("responses  $y_i$")) +
     xlab(TeX("isomap embeddings $\\hat{z}_i$")) 
p2


#plot as a whole, without regression line
# without separation by no. of claws
p3<-ggplot(df1, aes(x=z_hat, y=w)) +
  geom_point() +
  ylab(TeX("responses $y_i$ $\\rightarrow$")) +
  xlab(TeX("isomap embeddings $\\hat{z}_i$ $\\rightarrow$")) 
p3



#X1vsX2plot, sized by response y, for INTRO
YX12_hat<-cbind(X1,X2,y)
dfyx<-as.data.frame(YX12_hat)
p4<-ggplot(data = dfyx, aes(x=X1,y=X2,size=y)) +
  geom_point(alpha=0.5) +
  ylab(TeX("spectral embedding dimension $2$"))+
  xlab(TeX("spectral embedding dimension $1$"))+
  theme(legend.position = "none")
p4

grid.arrange(p4,p2,ncol=2)


#regression of y (distance)
linmody<-lm(y~z_hat,data = df1)
print(summary(linmody))


#regression of w (no. of claws)
linmodw<-lm(w~z_hat,data = df1)
print(summary(linmodw))