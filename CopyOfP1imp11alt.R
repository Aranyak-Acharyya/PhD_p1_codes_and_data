library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()

library(Matrix)
library(irlba)
library(vegan3d)
library(MASS)
library(igraph)
library(smacof)


e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d","smacof",
            "MASS","igraph",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))




d<-4

s<-20
ss<-s+1

n_vec<-seq(500,3500,250)

lambda_vec<-vector()
for(i in 1:length(n_vec))
{
  lambda_vec[i]<-0.85*0.99^(i-1)
}
print(lambda_vec)







alpha<-2.0
beta<-5.0
sig_ep<-0.01

lsg<-0.01
thres<-qf(lsg,1,s-2,lower.tail = FALSE)

u_vec<-seq(0,100,0.00001)

clusterExport(clust,list("d","s","ss",
                         "alpha","beta","sig_ep"))

RP<-matrix(,ncol=3)

power_diff_vec<-vector()

for(n in n_vec)
{
  
  lambda<-lambda_vec[which(n_vec==n)]
  
  
  l<-as.integer(n/10)
  
  clusterExport(clust,list("n","lambda","l"))
  
  
  
  
  B<-foreach(trial=1:100,.combine='c') %dopar%
    {
      #generating regressors
      ts<-runif(s,min=0,max=1)
      
      #generating pre-images for auxiliary latent positions
      tm<-runif(n-s,min=0,max=1)
      
      #combining to form set of all pre-images
      t<-c(ts,tm)
      
      
      #forming the matrix whose rows are latent positions
      X<-cbind(t/2,t/2,t/2,t/2)
      
      #forming the probability matrix
      P<-X%*%t(X)
      
      #generating the adjacency matrix
      pvec<-c(P[lower.tri(P,diag=FALSE)])
      avec<-rbinom(length(pvec),1,pvec)
      A<-matrix(nrow=nrow(P),ncol=ncol(P))
      A[lower.tri(A,diag=FALSE)]<-avec
      A[upper.tri(A)]<-t(A)[upper.tri(A,diag=FALSE)]
      diag(A)<-rep(0,n)
      #diag(A)<-apply(A,1,sum)/(n-1)
      
      
      #finding ASE estimates
      A_irlba<-irlba(A,d)
      X_hat<-A_irlba$u%*%diag(A_irlba$d)^0.5
      
      
      
      
      
      
      B0<-as.matrix(dist(X_hat[1:l,],method = "euclidean",
              diag=TRUE,upper=TRUE))
      
      BB<-ifelse(B0<lambda,B0,0)
      
      
      colnames(BB)<-as.character(seq(1,nrow(B0),1))
      
      
      
      g<-graph_from_adjacency_matrix(BB,
                                      mode="undirected",
                                      weighted=TRUE,
                                      diag=FALSE,
                                      add.colnames = NULL,
                                      add.rownames = NA)
      
     
      
      
      
      #matrix of shortest path distances
      D<-shortest.paths(g,v=V(g)[1:ss],to=V(g)[1:ss])
      Ds<-as.matrix(D)
      
      
      
      #raw-stress minimization
      MM<-mds(Ds,ndim = 1,type = "ratio",
          weightmat = NULL,
          init = "torgerson")
      
      
      #raw-stress embeddings
      z<-as.vector(MM$conf)
      zs<-z[1:s]
      z0<-z[ss]
      
      #generating responses
      ys<-alpha+beta*ts+rnorm(s,mean=0,sd=sig_ep)
      
      
      #true model
      beta_true<-cov(ys,ts)/var(ts)
      alpha_true<-mean(ys)-beta_true*mean(ts)
      y_true<-alpha_true+beta_true*t[s+1]
      
      #substitute model
      beta_sub<-cov(ys,zs)/var(zs)
      alpha_sub<-mean(ys)-beta_sub*mean(zs)
      y_sub<-alpha_sub+beta_sub*z0
        
      
      
      
      dec<-(y_true-y_sub)^2
      dec
      
      
    }
  

  
  new_row<-c(n,length(B),mean(B))
  print(new_row)
  
  

  RP<-rbind(RP,new_row)
}

stopCluster(clust)

RP<-RP[-1,]
df<-data.frame(RP)

save(df,file="P1new11alt.RData")
power_diff_vec<-RP[,3]

load("P1new11alt.RData")
RP<-as.matrix(df)


library(ggplot2)
library(reshape2)
library(latex2exp)

ggplot(data = df, aes(x=RP[,1],y=RP[,3])) +
  geom_point()+
  geom_line() +
  xlab(TeX("number of nodes (n)"))+
  ylab(TeX("sample E$(\\hat{y}_{sub}-\\hat{y}_{true})^2$")) +
  scale_x_continuous(
    breaks=round(seq(min(RP[,1]),max(RP[,1]),by=500),1)) 

ggsave(file="P1plot11alt.eps", 
       width = 5, height = 2,
       units = "in", dpi = 300)