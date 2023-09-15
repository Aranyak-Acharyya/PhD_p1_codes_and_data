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

n_vec<-seq(100,1000,50)
lambda_vec<-vector()
for(i in 1:length(n_vec))
{
  lambda_vec[i]<-0.9*0.99^(i-1)
}


alpha<-2.0
beta<-5.0
sig_ep<-0.01

lsg<-0.01
thres<-qf(lsg,1,s-2,lower.tail = FALSE)

u_vec<-seq(0,100,0.00001)

clusterExport(clust,list("d","s",
                         "alpha","beta","sig_ep"))

RP<-matrix(,ncol=3)

power_diff_vec<-vector()

for(n in n_vec)
{
  
  lambda<-lambda_vec[which(n_vec==n)]
  
  l<-n
  
  
  
  clusterExport(clust,list("n","lambda","l"))
  
  
  #opts<-list(preschedule=FALSE)
  
  B<-foreach(trial=1:100,.combine='rbind') %dopar%
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
      
      #finding ASE estimates
      A_irlba<-irlba(A,d)
      X_hat<-A_irlba$u%*%diag(A_irlba$d)^0.5
      
      
      B0<-as.matrix(dist(X_hat[1:l,],
                         method = "euclidean",
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
      D<-shortest.paths(g, v=V(g)[1:s],to=V(g)[1:s])
      Ds<-as.matrix(D)
      
      
      
      #raw-stress minimization
      MM<-mds(Ds,ndim = 1,type = "interval",
              weightmat = NULL,
          init = "torgerson")
      
    
      
      #raw-stress embeddings
      zs<-as.vector(MM$conf)
      
      
      #generating responses from regression model
      ys<-alpha+beta*ts+rnorm(s,mean=0,sd=sig_ep)
      
      #true regression parameter estimates
      beta_true<-cov(ys,ts)/var(ts)
      alpha_true<-mean(ys)-beta_true*mean(ts)
      
      #naive regression parameter estimates
      beta_naive<-cov(zs,ys)/var(zs)
      alpha_naive<-mean(ys)-beta_naive*mean(zs)
      
      #true predicted response
      ys_pred_true<-alpha_true+beta_true*ts
      
      #naive predicted responses
      ys_pred_naive<-alpha_naive+beta_naive*zs
      
      #true F-statistic
      num_true<-sum((ys_pred_true-mean(ys))^2)
      den_true<-sum((ys-ys_pred_true)^2)
      F_true<-(s-2)*(num_true/den_true)
      
      #naive F-statistic
      num_naive<-sum((ys_pred_naive-mean(ys))^2)
      den_naive<-sum((ys-ys_pred_naive)^2)
      F_naive<-(s-2)*(num_naive/den_naive)
      
      
      
      dec<-c(F_true,F_naive)
      dec
      
    
      
    }
  
  
  
  
  
  power_true<-mean(ifelse(B[,1]<thres,1,0))
  power_sub<-mean(ifelse(B[,2]<thres,1,0))
  
  new_row<-c(n,nrow(B),abs(power_true-power_sub))
  print(new_row)
  
  

  RP<-rbind(RP,new_row)
}

stopCluster(clust)




RP<-RP[-1,]
df<-data.frame(RP)
save(df,file="P1new18alt.RData")

load("P1new18alt.RData")
RP<-as.matrix(df)
power_diff_vec<-RP[,3]

library(ggplot2)
library(reshape2)
library(latex2exp)

ggplot(data = df, aes(x=RP[,1],y=RP[,3])) +
  geom_point()+
  geom_line() +
  theme(legend.title = element_blank()) +
  xlab(TeX("number of nodes (n)"))+
  ylab(TeX(" $|\\hat{\\pi}- \\pi^*|$")) 

ggsave(file="P1plot18alt.eps", 
       width = 5, height = 2,
       units = "in", dpi = 300)
