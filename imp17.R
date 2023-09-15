library("RColorBrewer")
library(tidyverse)
library(ggplot2)
library(GGally)
library(grid)
library(igraph)
library(gridExtra)

p_<-GGally::print_if_interactive
points_legend<-gglegend(ggally_points)

#loading the dataset
load("/cloud/project/df-forAranyak.RData")




n<-nrow(df)

X1<-df$X.1
X2<-df$X.2
X3<-df$X.3
X4<-df$X.4
X5<-df$X.5
X6<-df$X.6
y<-df$dist

X_hat<-cbind(X1,X2,X3,X4,X5,X6)
y<-df$dist
w<-df$claw
print(w)
claw<-as.character(df$claw)
print(X_hat)
dfx<-as.data.frame(X_hat)


#X1vsX2plot, sized by response y, for INTRO
YX12_hat<-cbind(X1,X2,y)
dfyx<-as.data.frame(YX12_hat)
p1<-ggplot(data = dfyx, aes(x=X1,y=X2,size=y)) +
     geom_point() +
     theme(legend.position = "none")


WX12_hat<-cbind(X1,X2,w)
dfwx<-as.data.frame(YX12_hat)
p1<-ggplot(data = dfwx, aes(x=X1,y=X2,size=y)) +
  geom_point() +
  theme(legend.position = "none")



  




eta<-0.6

#pairwise plot without colours
ggpairs(dfx,
        lower = list(continuous="points"),
        upper = list(continuous="points")
        ) +
  scale_x_continuous(guide = 
                       guide_axis(n.dodge = 1)) +
  scale_y_continuous(guide = 
                       guide_axis(check.overlap = 
                                    TRUE)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 


ggsave("P1plot17B.eps",
       height=4,width=4,
       units="in",dpi=300)


#pairwise plot with coloured by claws
ggpairs(df[,6:11],
        lower = list(continuous="points"),
        upper = list(continuous="points"),
        legend = 1,
        mapping=ggplot2::
          aes(colour = claw)) +
  scale_x_continuous(guide = 
                       guide_axis(n.dodge = 1)) +
  scale_y_continuous(guide = 
                       guide_axis(check.overlap = 
                                    TRUE)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 



X_hat<-cbind(df$X.1,df$X.2,df$X.3,df$X.4,df$X.5,df$X.6)
y<-df$dist



d<-6
#neighbourhood parameter
eta<-0.6

#construction of adjacency matrix of 
 #localozation graph
A<-matrix(nrow=n,ncol=n)
for(i in 1:n)
{
  for(j in 1:n)
  {
    if(i==j)
    {
      A[i,j]<-0
    }
    if(i!=j)
    {
      if(norm(X_hat[i,]-X_hat[j,],type="2")<eta)
      {
        A[i,j]<-1
      }
      if(norm(X_hat[i,]-X_hat[j,],type="2")>=eta)
      {
        A[i,j]<-0
      }
    }
  }
}
print(A)

#naming the columns of the adjacency matrix
nm<-vector()
for(j in 1:n)
{
  nm[j]<-as.character(j)
}
print(nm)
colnames(A)<-nm
print(colnames(A))

#forming the localization graph
 #from its adjacency matrix
g<-graph_from_adjacency_matrix(A,mode = "undirected",weighted = NULL,
                               diag=FALSE,add.colnames = NULL,
                               add.rownames = NA)
#is.connected(g,connected="weak")
#is.connected(g, connected = "weak", 
             #comp.dist.precomp = NULL)
print(g)
print(V(g))
tmp2 = shortest.paths(g, v='2')
print(tmp2)
tmp2[1,which(V(g)$name=='4')]

#forming the dissimilarity matrix of
 #shortest path graph distances
DM<-matrix(nrow=n,ncol=n)
for(i in 1:n)
{
  a<-as.character(i)
  tmp2 = shortest.paths(g, v=a)
  for(j in 1:n)
  {
   b<-as.character(j)
   DM[i,j]<-tmp2[1,which(V(g)$name==b)]
  }
}
print(DM)

#squared dissimilarity matrix
DMSq<-DM^2
print(DMSq)

d<-6
#extracting its eigenvalues
eigen(DMSq)
eval<-sort(abs(eigen(DMSq)$values),
           decreasing = TRUE)
yvec<-eval[1:d]
print(length(yvec))
xvec<-seq(1,d,1)
print(length(xvec))

dfnew<-data.frame(xvec,yvec)

ggplot(data = dfnew, aes(x=xvec,y=yvec)) +
  geom_point()+
  geom_line(colour = "red") +
  xlab(TeX("dimension \\rightarrow"))+
  ylab(TeX("eigenvalues \\rightarrow"))