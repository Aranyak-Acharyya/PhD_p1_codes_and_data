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

s<-15
ss<-s+1

n_vec<-seq(500,4000,250)

lambda_vec<-vector()
for(i in 1:length(n_vec))
{
  lambda_vec[i]<-1.5*0.99^(i-1)
}
print(lambda_vec)







alpha<-2.0
beta<-5.0
sig_ep<-0.01





clusterExport(clust,list("d","s","ss","qq",
                         "alpha","beta","sig_ep"))

RP<-matrix(,ncol=3)

power_diff_vec<-vector()


  
  n<-1500
  #lambda<-25.0
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
      
      
      
        
      
      
      
      
    
      
  
      
      #prediction by local linear regression
      df<-as.data.frame(cbind(X_hat[1:ss,],yss))
      df_train<-df[1:s,]
      df_test<-df[ss,]
      bw0<-npregbw(formula=yss~V1+V2+V3+V4,
                   regtype="ll",
                   bwmethod="cv.aic",
                   data=df_train)
      kre0<-npreg(bws=bw0)
      
      y_ase_kreg<-predict(kre0,
                          data=df_train,
                          newdata=df_test,
                          type="response")
      
      
      #prediction by KNN regression
      knn_mod<-knnreg(df_train[,-5],yss[1:s],k=3)
      y_ase_knn<-predict(knn_mod,df_test[,-5])
      
      
      
      
      
      dec<-c((y0-y_true)^2,
             (y0-y_sub)^2,
             (y0-y_ase_kreg)^2,
             (y0-y_ase_knn)^2)
      
      
      print(dec)
      
    }
  
  
  main<-apply(B,2,mean)
  
  new_row<-c(n,nrow(B),main)
  print(new_row)
  

stopCluster(clust)


data_row<-t(c(d,n,s,main))
print(data_row)


df<-data.frame(data_row)
colnames(df)<-c("d","n","s","true","sub",
                "llr","kNN")
print(df)
save(df,file="P1new10_rev1.RData")



load("P1new10_rev1.RData")
print(df)
df_new<-data.frame(data_row)
colnames(df_new)<-c("d","n","s","true","sub",
                    "llr","kNN")

df<-rbind(df,df_new)
save(df,file="P1new10_rev1.RData")


load("P1new10_rev1.RData")
print(df)




#B<-B[B[,3]<=50,]
B<-as.matrix(B)
B<-log(B)
print(B)
print(nrow(B))



library(ggplot2)
library(reshape2)
library(latex2exp)

tr_vec<-(1:nrow(B))
loss1_vec<-B[,1]
loss2_vec<-B[,2]
loss3_vec<-B[,3]
loss4_vec<-B[,4]


df<-data.frame(tr_vec,loss1_vec,loss2_vec,
               loss3_vec,loss4_vec)
dfm<-melt(df,id.var='tr_vec')
ggplot(data = dfm, aes(x=variable, y=value, 
                       fill = variable)) +
  geom_boxplot() +
  theme(legend.title = element_blank())  +
  scale_fill_manual(values = c("red", "orange", "green","blue"),
                    labels=unname(TeX(c("log$(\\hat{y}_{true}-y)^2$",
                                        "log$(\\hat{y}_{sub}-y)^2$",
                                        "log$(\\hat{y}_{llr}-y)^2$",
                                        "log$(\\hat{y}_{knn}-y)^2$"
                    )))) +
  labs(x="predictors",
       y="log of squared errors") +
  theme(axis.text.x = element_blank()) +
  theme(legend.box.background = element_rect(color="red")) 




ggsave(file="P1plot6_rev1.eps",height = 3,
       width = 5, units = "in", dpi = 300)





