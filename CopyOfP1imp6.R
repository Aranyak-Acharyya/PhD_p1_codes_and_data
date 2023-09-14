library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()



library(Matrix)
library(MASS)
library(irlba)
library(princurve)


e <- new.env()
e$libs <- c("Matrix","MASS","irlba","princurve",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))






p<-seq(0,1,length.out=10^5)
CC<-cbind(p^2,2*p*(1-p),(1-p)^2)
d<-3
n_vec<-seq(600,2500,100)


beta<-5.0
alpha<-2.0
sig_ep<-0.1

clusterExport(
  clust,list("d","p","CC","alpha","beta","sig_ep"))

RS<-matrix(,ncol=4)

for(n in n_vec)
{
  
  clusterExport(clust,"n")
  
  L<-foreach(i=1:100,.combine = 'rbind') %dopar%
    {
      #generating regressors
      t<-runif(n,min=0,max=1)
      
      #generating response variable values
      e<-rnorm(n,mean=0,sd=sig_ep)
      y<-alpha+beta*t+e
      
      
      #forming the matrix of latent positions
      X<-cbind(t^2,2*t*(1-t),(1-t)^2)
      
      #forming probability matrix
      P<-X%*%t(X)
      
      #generating adjacency matrix
      pvec<-c(P[lower.tri(P,diag=FALSE)])
      avec<-rbinom(length(pvec),1,pvec)
      A<-matrix(nrow=nrow(P),ncol=ncol(P))
      A[lower.tri(A,diag=FALSE)]<-avec
      A[upper.tri(A)]<-t(A)[upper.tri(A,diag=FALSE)]
      diag(A)<-rep(0,n)
      
      
      
      #finding adjacency spectral embedding
      A_irlba<-irlba(A,d)
      X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
      
      #rotating ASE to find consistent estimates of latent positions
      vals<-svd(t(X_hat_raw)%*%X)
      W<-vals$u%*%t(vals$v)
      X_hat<-X_hat_raw%*%W
      
      
      #projecting ASE estimates on manifold 
      proj<-project_to_curve(X_hat,CC)
      X_tilde<-proj$s
      
      
      #estimating regressors from projections
      t_hat<-0.5*X_tilde[,2]+X_tilde[,1]
      
      
      t_ct<-t-mean(t)
      t_hat_ct<-t_hat-mean(t_hat)
      
      reg_par<-c(alpha,beta)
      
      
      #true estimator
      beta_true<-(t(t_ct)%*%y)/(t(t_ct)%*%t_ct)
      alpha_true<-mean(y)-beta_true*mean(t)
      reg_par_true<-c(alpha_true,beta_true)
      loss1<-(norm(reg_par_true-reg_par,type = "2"))^2
      loss1a<-(alpha_true-alpha)^2
      loss1b<-(beta_true-beta)^2
      
      
      
      #naive estimator
      beta_naive<-(t(t_hat_ct)%*%y)/(t(t_hat_ct)%*%t_hat_ct)
      alpha_naive<-mean(y)-beta_naive*mean(t_hat)
      reg_par_naive<-c(alpha_naive,beta_naive)
      loss2<-(norm(reg_par_naive-reg_par,type = "2"))^2
      loss2a<-(alpha_naive-alpha)^2
      loss2b<-(beta_naive-beta)^2
      
      
      
      
      
      dec<-c(loss1,loss2)
      dec
      
    
     
      
    }
  
  risk_all<-apply(L,2,mean)
  new_row<-c(n,nrow(L),risk_all)
  RS<-rbind(RS,new_row)
  
  
  print(new_row)
  
  
}

stopCluster(clust)


RS<-RS[-1,]
df<-data.frame(RS)
save(df,file = "P1new6.RData")
load("P1new6.RData")
RS<-as.matrix(df)



risk1_vec<-RS[,3]
risk2_vec<-RS[,4]


library(ggplot2)
library(reshape2)
library(latex2exp)

df<-data.frame(n_vec,risk1_vec,risk2_vec)
dfm<-melt(df, id.vars = 'n_vec')
print(dfm)
ggplot(dfm, aes(x=n_vec, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("sample MSEs")) +
  xlab(TeX("number of nodes(n)")) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = c("red","orange"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\theta}_{true}$",
                        "sample MSE of $\\hat{\\theta}_{sub}$"))))



ggsave(file="P1plot6.eps", width = 5, height = 2,
       units = "in", dpi = 300)