# This script is adapted from Yuan Yuan 2012-12-13
# This script calculates screen for the significant features against the clinical outcome, using univariate cox model
library(foreach)
source("models.R")
cox.screen <- function(y.train,x.train, pvalue.cutoff=0.05, qvalue.cutoff=0.2, top=100) 
  {
    # calculate the p-value from likelihood ratio test of univariate cox model
    # for categorical variable, such as mutation, should use log-rank test
    x.train<-as.matrix(x.train)
    feature.pvalue <- c()
    feature.name <- c()
    feature.col <- c()
    
      
    result<-c()
    for (j in 1:ncol(x.train))
    {
      feature <- colnames(x.train)[j]
      x <- x.train[,j]
      cox <- try(coxph(y.train~x))
      if (class(cox)=="try-error")
        {
          stop() #next# note next does not work with foreach %dopar%    
        }
      p.value <-summary(cox)$logtest["pvalue"] # likelihood ratio test
      #print(p.value)
      result<-rbind(result,list(p.value, feature, j))
    }

    feature.pvalue <- unlist(result[,1])
    print(paste("Total number of valid features (after removal of potential flat records):", length(feature.pvalue)))
    feature.name <- unlist(result[,2])
    feature.col <- unlist(result[,3])
          
    names(feature.pvalue)=c()
    names(feature.name)=c()
    names(feature.col)=c()
    
    q.value <- p.adjust(feature.pvalue, method="fdr")
    for (i in 1:10/10)
      {
        print(paste("qvalue <=", i, ":", length(which(q.value<=i))))
      }
    col.sig <- which( feature.pvalue < pvalue.cutoff) # & q.value < qvalue.cutoff  )
    print(paste("Significant records: p-value < 0.05:", length(col.sig)))  #and FDR <",qvalue.cutoff,":", length(col.sig)))
   
    name.sig <- feature.name[col.sig]
    pvalue.sig <- feature.pvalue[col.sig]
    
    col.retain <- feature.col[col.sig]
    if (length(col.retain)> top) # only keep the top significant ones if there are too many
      {
        col.retain <- col.retain[head(sort(pvalue.sig, index.return=TRUE)$ix, n=top)]
      }
    
    return(col.retain)
}


screening.cox<-function(pocket,embs,pofFlag){
 if(pofFlag==1){
   dataX <- embs
   f.names<-colnames(dataX)
   mySurv<-pocket$mySurv
   
   tr.inds<-pocket$BOO$tr.inds
   N<-dim(tr.inds)[1]
   
   fset<-matrix(0,N,dim(dataX)[2])
   for(i in c(1:N)){
     y.train<-mySurv[tr.inds[i,],]
     x.train<-dataX[tr.inds[i,],]
     fset[i,cox.screen(y.train,x.train,top=8515)]<-1
   } 
   colnames(fset)<-f.names
   return(fset)
 }
  else
  {
    dataX<-pocket$data

    f.names<-colnames(dataX)
    mySurv<-pocket$mySurv
    
    tr.inds<-pocket$BOO$tr.inds
    N<-dim(tr.inds)[1]
    
    fset<-matrix(0,N,dim(dataX)[2])
    for(i in c(1:N)){
      y.train<-mySurv[tr.inds[i,],]
      x.train<-dataX[tr.inds[i,],]
      fset[i,cox.screen(y.train,x.train,top=8515)]<-1
    } 
    colnames(fset)<-f.names
    return(fset)
  }
}




