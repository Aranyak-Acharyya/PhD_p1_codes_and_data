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
e$libs <- c("irlba","Matrix","princurve",
            "MASS",
            .libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))


n<-800



p<-seq(0,1,length.out=10^5)
CC<-cbind(p^2,2*p*(1-p),(1-p)^2)
d<-3
sig_ep<-0.1


beta<-5.0


clusterExport(clust,list("d","n","p","CC","beta","sig_ep"))




U<-foreach(trial=1:100,.combine='cbind') %dopar%
  {
    #pre-images of latent positions
    t<-runif(n,min=0,max=1)
    
    
    
    #forming matrix of latent positions
    X<-cbind(t^2,2*t*(1-t),(1-t)^2)
    
    #probability matrix
    P<-X%*%t(X)
    
    
    
    #generating adjacency matrix
    pvec<-c(P[lower.tri(P,diag=TRUE)])
    avec<-rbinom(length(pvec),1,pvec)
    A<-matrix(nrow=nrow(P),ncol=ncol(P))
    A[lower.tri(A,diag=TRUE)]<-avec
    A[upper.tri(A)]<-t(A)[upper.tri(A)]
    
    
    #finding adjacency spectral embedding
    A_irlba<-irlba(A,d)
    X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
    
    
    #rotating to find consistent ASE estimates
    vals<-svd(t(X_hat_raw)%*%X)
    W<-vals$u%*%t(vals$v)
    X_hat<-X_hat_raw%*%W
    
    
    #projecting ASE estimates on manifold
    proj<-project_to_curve(X_hat,CC)
    X_tilde<-proj$s
    
    
    
    
    #estimating regressors from projections
    t_hat<-0.5*X_tilde[,2]+X_tilde[,1]
    
    u<-(t_hat-t)
    u
    
  }


#stopCluster(clust)
sum_Gamma<-sum(apply(U,1,var))

clusterExport(clust,"sum_Gamma")

R<-foreach(i=1:100,.combine = 'rbind') %dopar%
  {
    #generating regressors
    t<-runif(n,min=0,max=1)
    
    #generating response variable values
    e<-rnorm(n,mean=0,sd=sig_ep)
    y<-beta*t+e
    
    
    #forming the matrix of latent positions
    X<-cbind(t^2,2*t*(1-t),(1-t)^2)
    
    #forming probability matrix
    P<-X%*%t(X)
    
    #generating adjacency matrix
    pvec<-c(P[lower.tri(P,diag=TRUE)])
    avec<-rbinom(length(pvec),1,pvec)
    A<-matrix(nrow=nrow(P),ncol=ncol(P))
    A[lower.tri(A,diag=TRUE)]<-avec
    A[upper.tri(A)]<-t(A)[upper.tri(A)]
    
    
    
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
    
    
    
    delta<-matrix(0,d,d)
    
    
    for(j in 1:n)
    {
      delta<-delta+(X_hat[j,]%*%t(X_hat[j,]))
    }
    Delta<-(1.0/n)*delta
    
    Delinv<-ginv(Delta,tol = sqrt(.Machine$double.eps))
    
    incr_vec<-vector()
    for(i in 1:n)
    {
      sum_mat<-matrix(0,d,d)
      
      grad_vec<-c(1,0.5,0)
      
      for(j in 1:n)
      {
        
        F1<-(t(X_hat[i,])%*%X_hat[j,])*(1-t(X_hat[i,])%*%X_hat[j,])
        f1<-as.double(F1)
        num<-f1*(X_hat[j,]%*%t(X_hat[j,]))
        sum_mat<-sum_mat+num
      }
      
      Sum_mat<-(1.0/n)*sum_mat
      
      incr<-t(grad_vec)%*%Delinv%*%Sum_mat%*%Delinv%*%grad_vec
      
      incr_vec<-c(incr_vec,incr)
      
    }
    
    sum_Gamma_hat<-mean(incr_vec)
    
    
    
    #true estimator
    beta_true<-(t(t)%*%y)/(t(t)%*%t)
    loss1<-(beta_true-beta)^2
    
    
    
    #naive estimator
    beta_naive<-(t(t_hat)%*%y)/(t(t_hat)%*%t_hat)
    loss2<-(beta_naive-beta)^2
    
    
    #ME-adjusted estimator with known error variance
    beta_adj_par<-(t(t_hat)%*%y)/(t(t_hat)%*%t_hat-sum_Gamma)
    loss3<-(beta_adj_par-beta)^2

    #ME-adjusted estimator with estimated unknown variance
    beta_adj_par_hat<-(t(t_hat)%*%y)/(t(t_hat)%*%t_hat-sum_Gamma_hat)
    loss4<-(beta_adj_par_hat-beta)^2
    
    
    dec<-c(beta_true,beta_naive,beta_adj_par,beta_adj_par_hat)
    dec
   
  }

print(beta)

stopCluster(clust)

dff<-data.frame(R)
save(dff,file="P1new3.RData")
load("P1new3.RData")
R<-as.matrix(dff)


tr_vec<-1:nrow(R)
beta1_vec<-R[,1]
beta2_vec<-R[,2]
beta3_vec<-R[,3]
beta4_vec<-R[,4]

risk1_vec<-(beta1_vec-beta)^2
risk2_vec<-(beta2_vec-beta)^2
risk3_vec<-(beta3_vec-beta)^2
risk4_vec<-(beta4_vec-beta)^2

library(ggplot2)
library(reshape2)
library(latex2exp)

df<-data.frame(tr_vec,risk1_vec,risk2_vec,risk3_vec,
               risk4_vec)
print(df)
dfm<-melt(df, id.var='tr_vec')
print(dfm)
ggplot(data = dfm, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot() +
  theme(legend.title = element_blank())  +
  scale_fill_manual(values = c("red", "orange", "green","blue"),
                    labels=unname(TeX(c("$(\\hat{\\beta}_{true}-\\beta)^2$",
                                        "$(\\hat{\\beta}_{naive}-\\beta)^2$",
                                        "$(\\hat{\\beta}_{adj,\\sigma}-\\beta)^2$",
                                        "$(\\hat{\\beta}_{adj,\\hat{\\sigma}}-\\beta)^2$"
                    )))) +
  labs(x="regression parameter estimators",
       y="squared errors") +
  theme(axis.text.x = element_blank()) +
  theme(legend.box.background = element_rect(color="red")) 



ggsave(file="P1plot3.eps",height = 2,
       width = 5, units = "in", dpi = 300)
















