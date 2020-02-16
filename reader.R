
reader<-function(nf,wd,folder,saved.filename){
  # nf=131
  file=paste(wd,folder,"/genomicMatrix",sep = "");
  file.cil=paste(wd,folder,"/clinical_data",sep = "");
  
  # read  file 
  datx<- read.table(file,nrows=nf, header=T,row.names=1,sep="\t", quote="")
  datx<-t(datx)
  
  clinical <- read.table(file.cil,header=T, row.names=1, sep="\t",quote="")
  
  # # match the patient ID in clinical data with the colnames of z_rna
  clinical<-t(data.frame(t(clinical)))
  
  
  sum(rownames(clinical) %in% rownames(datx)) 
  ind_sample=which(rownames(clinical) %in% rownames(datx)==TRUE)
  
  # get the columns that contain data we can use: days to death, event
  ind_keep <-grep("X_OS",colnames(clinical))
  
  clinical=clinical[ind_sample,ind_keep]
  name.rm=rownames(clinical)[which(is.na(clinical[,1]))]
  
  datx=datx[which(is.na(match(rownames(datx),name.rm))),]
  clinical=clinical[which(is.na(match(rownames(clinical),name.rm))),]
  
  clinical=clinical[order(rownames(clinical)),]
  datx=datx[order(rownames(datx)),]
  
  clinical=matrix(as.numeric(unlist(clinical)),nrow=nrow(clinical))
  
  ind0=which(clinical[,1]>0)
  
  clinical=clinical[ind0,]
  datx=datx[ind0,]
  
  
  library(survival)
  mySurv <- Surv(clinical[,1], clinical[,2],type='right')
  data=datx
  
  library(randomForest)
  data <- na.roughfix(data)
  
  ind<-rem(t(data))
  
  if(length(ind)!=0){data<-data[,-ind]}

  
  re<-list()
  re$data<-data
  re$mySurv<-mySurv
  re$folder<-folder
  re$dims<-dim(data)
  re$c.e<-c(sum(1-mySurv[,2]),sum(mySurv[,2]))
  
  k<-floor(re$c.e/5)
  bo.size<-c(4*k,k)
  re$BOO<-boo(re$mySurv, bo.size)
  
  #save(file=saved.filename,re)
  
  return(re)
}

#remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x[is.na(x)]<-0
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

boo<-function(mySurv,bo.size){
  #bo.size=[tr.c,tr.e,ts.c,ts.e]
  #va.size=[va.c,va.e]
  
  c.inds<-which(mySurv[,2]==0)
  e.inds<-which(mySurv[,2]==1)
  print(c(length(c.inds),length(e.inds)))
  
  tr.inds<-matrix(0,100,sum(bo.size[c(1,2)]))
  ts.inds<-matrix(0,100,sum(bo.size[c(3,4)]))
  for(i in c(1:100)){
    set.seed(i)
    cts<-sample(c.inds,bo.size[3])
    ctr<-sample(setdiff(c.inds,cts), bo.size[1])
    set.seed(-1*i)
    ets<-sample(e.inds,bo.size[4])
    etr<-sample(setdiff(e.inds,ets),bo.size[2])
    tr.inds[i,]<-c(ctr,etr)
    ts.inds[i,]<-c(cts,ets)
  }
  poo<-list()
  poo$tr.inds<-tr.inds
  poo$ts.inds<-ts.inds
  poo$split.size<-c(sum(bo.size[c(1,2)]),sum(bo.size[c(3,4)]))
  return(poo)
}

booer<-function(zubi){
  #bo.size=[tr.c,tr.e,ts.c,ts.e]
  k<-floor(zubi$c.e/5)
  bo.size<-c(4*k,k)
  return(boo(zubi$mySurv, bo.size))
}

l1co<-function(data,mySurv,seed,tri,tsi){
  
  x.train=data[tri,]
  y.train=mySurv[tri,]
  x.test=data[tsi,]
  y.test=mySurv[tsi,]
  
  coefs<-0*(1:dim(data)[2])
  clnm<-colnames(data)
  
  cols.include<-pres(x.train,sum(y.train[,2]),0.8)
  
  x.train <- x.train[,cols.include]
  x.test <- x.test[,cols.include]
  print(paste("After univariate cox screen, features remain:", length(cols.include)))
  
  if(length(cols.include) > 5) #useLASSO &  further shrink by LASSO, if only a few features, no need to use LASSO, note 5 is a quite arbitrary setting
  {
    # do LASSO without cross validation to get the features to include in the model
    # change x.train and x.test to only include those features
    cols.include <- c()
    iter <- 0
    while(length(cols.include)<1)
    {
      set.seed(iter+1)
      cols.include <- try(lasso(x=x.train,y=y.train,above=0))
      if (class(cols.include)=="try-error")
      {
        print("Errors occur while calculating by LASSO, recalculating...")
        cols.include <- c()
      }
      iter <- iter+1
      if(iter> 1)
      {
        print(paste(length(cols.include)," features selected. Recalculated by LASSO:", iter))
      }
      if(iter>100) # maximum number of iterations allowed
      {
        stop("No significant features can be selected by LASSO: exit.")
      }
    }    
    print(paste("After LASSO, features remain:", length(cols.include)))
    x.train <- x.train[,cols.include]
    x.test <- x.test[,cols.include]
  }
  print(paste("seed =", seed, ": final features:"))
  print("------------------------------------------------")
  print(colnames(x.train))
  print("------------------------------------------------")
  
  # cox model for prediction
  if (length(cols.include)==1)
  {
    cox <- coxph(y.train~x.train)
  }else
  {
    cox <- coxph(y.train~., data= data.frame(x.train))
  }
  
  cox$coefficients[is.na(cox$coefficients)]=0 # convert NAs to zero if any       
  cox.predict <- as.matrix(x.test)%*%cox$coefficients
  
  coefs[match(colnames(x.train),clnm)]<-cox$coefficients #here3
  #fxx[which(abs(coefs[seed,])>0)]<-3#here4
  
  c1<-concordance.index(cox.predict,y.test[,1], y.test[,2])    
  cx=c1$c.index
  print(paste("seed ",seed,": ",c1$c.index))
  ree<-list()
  ree$C
  ree$coefs
  return(ree)
}



