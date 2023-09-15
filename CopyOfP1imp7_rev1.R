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
library(np)
library(caret)
library(neuralnet)


e <- new.env()
e$libs <- c("irlba","Matrix","vegan3d","smacof",
            "np","caret","neuralnet",
            "MASS","igraph",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))


qq<-1/(2^0.5)

d<-4

s<-5
ss<-s+1

n_vec<-seq(500,4000,250)

lambda_vec<-vector()
for(i in 1:length(n_vec))
{
  lambda_vec[i]<-0.8*0.99^(i-1)
}
print(lambda_vec)

grad_desc<-function(X,y,learn_rate,threshold,
                    max_iter)
{
  d<-ncol(X)
  Dat<-X
  beta<-runif(d,min=0,max=1)
  loss<-t(y-Dat%*%beta)%*%(y-Dat%*%beta)
  converged=F
  iter=0
  while(converged==F)
  {
    beta_new<-beta-
      learn_rate*(2*t(Dat)%*%Dat%*%beta-2*t(Dat)%*%y)
    beta<-beta_new
    loss_new<-t(y-Dat%*%beta)%*%(y-Dat%*%beta)
    if(loss-loss_new<=threshold)
    {
      converged=T
    }
    iter=iter+1
    if(iter>max_iter)
    {
      converged=T
    }
  }
  return(beta)
}








alpha<-2.0
beta<-5.0
sig_ep<-0.01





clusterExport(clust,list("d","s","ss","qq",
                         "alpha","beta","sig_ep",
                         "grad_desc"))

RP<-matrix(,ncol=3)

power_diff_vec<-vector()


  
  n<-500
  lambda<-lambda_vec[which(n_vec==n)]
  
  l<-n
  #l<-as.integer(n/2)
  
  
  clusterExport(clust,list("n","lambda","l"))
  
  
  
  
  B<-foreach(trial=1:100,.combine='rbind') %dopar%
    {
      #generating regressors
      ts<-runif(s,min=0,max=1)
      
      #generating pre-images for auxiliary latent positions
      tm<-runif(n-s,min=0,max=1)
      
      #combining to form set of all pre-images
      t<-c(ts,tm)
      
      
      #forming the matrix whose rows are latent positions
      X<-cbind(qq*cos(t),
               qq*sin(t),
               qq*cos(t),
               qq*sin(t))
      
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
      zss<-z[1:ss]
      
      #generating responses
      tss<-t[1:ss]
      yss<-alpha+beta*tss+rnorm(ss,mean=0,sd=sig_ep)
      ys<-yss[1:s]
      
      y0<-yss[ss]
      
      
      #prediction by true regressors
      df<-as.data.frame(cbind(tss,yss))
      df_train<-df[1:s,]
      df_test<-df[ss,]
      
      
      
      modt<-lm(yss~tss,data = df_train)
      y_true<-predict(modt,df_test)
    
      
      #prediction by embeddings
      df<-as.data.frame(cbind(zss,yss))
      df_train<-df[1:s,]
      df_test<-df[ss,]
      
      modz<-lm(yss~zss,data=df_train)
      y_sub<-predict(modz,df_test)
      
      
      
        
      
      #prediction by one-stage 
      beta_hat<-grad_desc(A[1:s,],yss[1:s],0.01,
                          10^(-3),
                          10^2)
      
      
      y_os<-A[ss,]%*%beta_hat
      
      
      
    
      
      
      
      
      
      
      
      
      
      
      
      
      
      dec<-c((y0-y_true)^2,
             (y0-y_sub)^2,
             (y0-y_os)^2)
      
      
      
      
    }
  
  
  main<-apply(B,2,mean)
  
  new_row<-c(n,nrow(B),main)
  print(new_row)
  
  

  

stopCluster(clust)



dff<-data.frame(B)
save(dff,file="P1new7_rev1.RData")

load("P1new7_rev1.RData")


colnames(dff)<-NULL
rownames(dff)<-NULL






B<-as.matrix(dff)
B<-log(B)





library(ggplot2)
library(reshape2)
library(latex2exp)

tr_vec<-(1:nrow(B))
loss1_vec<-B[,1]
loss2_vec<-B[,2]
loss3_vec<-B[,3]


df<-data.frame(tr_vec,loss1_vec,loss2_vec,
               loss3_vec)
dfm<-reshape2::melt(df,id.var='tr_vec')
ggplot(data = dfm, aes(x=variable, y=value, 
                       fill = variable)) +
  geom_boxplot() +
  theme(legend.title = element_blank())  +
  scale_fill_manual(values = c("red", "orange", "green"),
                    labels=unname(TeX(c("log$(\\hat{y}_{true}-y)^2$",
                                        "log$(\\hat{y}_{sub}-y)^2$",
                                        "log$(\\hat{y}_{os}-y)^2$"
                    )))) +
  labs(x="predicted responses",
       y="log of squared errors") +
  theme(axis.text.x = element_blank()) +
  theme(legend.box.background = element_rect(color="red")) 




ggsave(file="P1plot7_rev1.eps",height = 3,
       width = 5, units = "in", dpi = 300)





