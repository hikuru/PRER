library(foreach)
library(survcomp);
require(survival);


pres<-function(data,top,co.th){

  md<-apply(data,2,mad)
  col.retain<-head(sort(md,index.return=T,decreasing = T)$ix,n=top)
  
  dati<-data[,col.retain]
  cr<-cor(dati)
  k<-which(cr>co.th,arr.ind = T)
  k<-k[which(k[,1]>k[,2]),]
  
  
  if(is.null(dim(k))){
    col.retain<-col.retain[-unique(k)]
    return(col.retain)
  }else{  
    col.retain<-col.retain[-unique(k[,1])]
    return(col.retain)
  }
}

rsf<- function(re,embs,fset=0,coxscr,pofFlag){
  require(randomForestSRC)
  
  ntree=1000
  dataa<-embs
  mySurv<-re$mySurv
  tri<-re$BOO$tr.inds
  tsi<-re$BOO$ts.inds

  l<-dim(tri)[1]
 
  
  if(pofFlag == 1){
    clnm <- colnames(embs)
    cx <-NA*c(1:l)
    
    coefs<-matrix(0, l, ncol(embs))

    for(seed in 1:l)
    { 
      x.train<-dataa[tri[seed,],]
      colnames(x.train)<-clnm
      
      x.train<-x.train[,which(as.numeric(coxscr[seed,])==1)]
      y.train<-mySurv[tri[seed,],]
      fs <- dim(x.train)[2]
      x.test<-dataa[tsi[seed,],]
      
      colnames(x.test)<-clnm
      x.test<-x.test[,which(as.numeric(coxscr[seed,])==1)]
      y.test=mySurv[tsi[seed,],]

      data <- data.frame(time=y.train[,1],status=y.train[,2])
      data <- cbind(data,x.train)
      
     
      print('here')
      rf <- rfsrc(Surv(time, status)~., data=data, ntree=ntree, seed=-1,importance = T)
      
      feature.imp <- which(rf$importance>0)
      
      print(paste("After random forest, features remain:", length(feature.imp)))
      print("------------------------------------------------")
      
      
      rf$importance[which(is.na(rf$importance))]<-0 #here final
      
      coefs[seed,match(names(feature.imp),clnm)]<-rf$importance[which(rf$importance>0)] #here3
      colnames(coefs)<-clnm
      data.test <- data.frame(time=y.test[,1],status=y.test[,2])
      data.test <- cbind(data.test,x.test)
      rsf.predict <- predict(rf, data.test[,3:(fs+2)], seed=-1,na.action="na.impute")$predicted
      
      c1<-concordance.index(rsf.predict,y.test[,1], y.test[,2])    
      cx[seed]=c1$c.index
      
      print(paste("seed ",seed,": ",c1$c.index))
    }
  }
  else{
    com <- colnames(dataa)
    cx <-NA*c(1:l)
    
    clnm<-colnames(dataa)
    coefs<-matrix(0, l, 131)
    for(seed in 1:l)
    { 
      x.train<-dataa[tri[seed,],]
      colnames(x.train)<-clnm
      x.train<-x.train[,which(as.numeric(coxscr[seed,])==1)]
      y.train<-mySurv[tri[seed,],]
      fs <- dim(x.train)[2]

      x.test<-dataa[tsi[seed,],]
      
      colnames(x.test)<-clnm
      x.test<-x.test[,which(as.numeric(coxscr[seed,])==1)]
      y.test=mySurv[tsi[seed,],]

      
      if(is.null(fs)){
        data <- data.frame(time=y.train[,1],status=y.train[,2],x.train=x.train)
        
        print('here')
        
        rf <- rfsrc(Surv(time, status)~., data=data, ntree=ntree, seed=-1,importance = T)
        
        feature.imp <- which(rf$importance>0)
        
        print(paste("After random forest, features remain:", length(feature.imp)))
        print("------------------------------------------------")
        
        
        rf$importance[which(is.na(rf$importance))]<-0 #here final
        
        coefs[seed,match(names(feature.imp),clnm)]<-rf$importance[which(rf$importance>0)] #here3
        colnames(coefs)<-clnm
        data.test <- data.frame(time=y.test[,1],status=y.test[,2])
        data.test <- cbind(data.test,x.test)
        cx[seed]= c1$c.index
      }
      else{
        data <- data.frame(time=y.train[,1],status=y.train[,2])
        data <- cbind(data,x.train)
        
        print('here')
        
        rf <- rfsrc(Surv(time, status)~., data=data, ntree=ntree, seed=-1,importance = T)
        
        feature.imp <- which(rf$importance>0)
        
        print(paste("After random forest, features remain:", length(feature.imp)))
        print("------------------------------------------------")
        
        
        rf$importance[which(is.na(rf$importance))]<-0 #here final
        
        coefs[seed,match(names(feature.imp),clnm)]<-rf$importance[which(rf$importance>0)] #here3
        colnames(coefs)<-clnm
        data.test <- data.frame(time=y.test[,1],status=y.test[,2])
        data.test <- cbind(data.test,x.test)
        rsf.predict <- predict(rf, data.test[,3:(fs+2)], seed=-1,na.action="na.impute")$predicted
        c1<-concordance.index(rsf.predict,y.test[,1], y.test[,2])    
        cx[seed]=c1$c.index
      }
      
      print(paste("seed ",seed,": ",c1$c.index))
    }
  }
  ree<-list()
  ree$coefs<-coefs
  ree$cinds<-cx
  
  print("---END---")#print(paste(colnames(x.train[,1]),"---END---"))
  
  return(ree)
}

comparator3 <- function(rppa, dd,randomWalks,nnetwork) {
  rppa2 <- data.frame(uniprot=rppa$uniprot, num = numeric(length = nrow(rppa)),count = numeric(length = nrow(rppa)),neig = vector(length = nrow(rppa)),count2 = numeric(length = nrow(rppa)))
  for (i in 1:nrow(rppa2)) {
    gen <- rppa2$uniprot[i]
    ind <- as.character(unlist(nnetwork$uniprots)) == gen
    rppa2$num[i] <- nnetwork$nodeNumber[ind]
  }
  nn <- ncol(randomWalks)
  for (i in 1:nrow(rppa2)) {
    genNum <- rppa2$num[i]
    genWalks <- randomWalks[,1] == genNum
    walks <- as.numeric(unlist(randomWalks[genWalks,2:nn]))
    
    tt <- table(intersect(walks,rppa2$num))
    tt <- as.numeric(names(tt))
    neig <- c()
    cc = 0
    for (j in 1:length(tt)) {
      if (sum(walks == tt[j]) > 1) {
        neig <- append(neig, tt[j])
        cc <- cc + 1
      }
    }
    zz <- neig != genNum
    neig<-neig[zz]
    xx <- rppa2$num %in% neig
    rppa2$neig[i] <- toString(unique(uniprots[xx]))
    rppa2$count[i] <- cc
    rppa2$count2[i] <- sum(xx)
  }
  
  nc <- sum(rppa2$count2)
  nr <- nrow(dd)
  rs = matrix(0,nr,nc)
  cur.ind <- 1
  clnm <- c()
  
  for (i in 1:nrow(rppa2)) {
    genNum <- rppa2$num[i]
    genWalks <- randomWalks[,1] == genNum
    walks <- as.numeric(unlist(randomWalks[genWalks,2:nn]))
    
    tt <- table(intersect(walks,rppa2$num))
    tt <- as.numeric(names(tt))
    neig <- c()
    for (j in 1:length(tt)) {
      if (sum(walks == tt[j]) > 1) {
        neig <- append(neig, tt[j])
      }
    }
    zz <- neig != genNum
    neig<-neig[zz]
    inds <- which(rppa2$num %in% neig)
    for (ii in 1:length(inds)) {
      cname <- paste(rppa$RPPA[i], rppa$RPPA[inds[ii]], sep = "-")
      clnm <- append(clnm, cname)
      temp <- rep(1,nr)
      ind <- which(dd[,i] > dd[,inds[ii]])
      temp[ind] <- -1
      rs[,cur.ind] <- temp
      cur.ind <- cur.ind+1
    }
  }
  colnames(rs) <- clnm
  
  unNum <- unique(rppa2$num)
  coms <- 0
  clnm2 <- c()
  for(i in 1:length(unNum)){
    co <- sum(rppa2$num == unNum[i])
    if(co>1){
      coms <- coms + (co * (co-1) / 2)
    }
  }
  rs2 = matrix(0,nr,coms)
  cur.ind <- 1
  for(i in 1:length(unNum)){
    co <- sum(rppa2$num == unNum[i])
    if(co>1){
      ind.m <- which(rppa2$num == unNum[i])
      combina <- combn(ind.m,2)
      for (j in 1:ncol(combina)) {
        xxx <- as.numeric(combina[,j])
        cname <- paste(rppa$RPPA[xxx[1]], rppa$RPPA[xxx[2]], sep = "-")
        clnm2 <- append(clnm2, cname)
        temp <- rep(1,nr)
        ind <- which(dd[,xxx[1]] > dd[,xxx[2]])
        temp[ind] <- -1
        rs2[,cur.ind] <- temp
        cur.ind <- cur.ind+1
      }
    }
  }
  colnames(rs2) <- clnm2
  rs <- cbind(rs,rs2)
  return(rs)
}

comparator2 <- function(rppa, dd, randomWalks) {
  
  rppa2 <- data.frame(uniprot=rppa$uniprot, num = numeric(length = nrow(rppa)),count = numeric(length = nrow(rppa)),neig = vector(length = nrow(rppa)),count2 = numeric(length = nrow(rppa)))
  for (i in 1:nrow(rppa2)) {
    gen <- rppa2$uniprot[i]
    ind <- as.character(unlist(nnetwork$uniprots)) == gen
    rppa2$num[i] <- nnetwork$nodeNumber[ind]
  }
  
  for (i in 1:nrow(rppa2)) {
    genNum <- rppa2$num[i]
    genWalks <- randomWalks[,1] == genNum
    walks <- as.numeric(unlist(randomWalks[genWalks,2:20]))
    
    tt <- table(intersect(walks,rppa2$num))
    tt <- as.numeric(names(tt))
    neig <- c()
    cc = 0
    for (j in 1:length(tt)) {
      if (sum(walks == tt[j]) > 1) {
        neig <- append(neig, tt[j])
        cc <- cc + 1
      }
    }
    zz <- neig != genNum
    neig<-neig[zz]
    xx <- rppa2$num %in% neig
    rppa2$neig[i] <- toString(unique(uniprots[xx]))
    rppa2$count[i] <- cc
    rppa2$count2[i] <- sum(xx)
  }
  
  nc <- sum(rppa2$count2)
  nr <- nrow(dd)
  rs = matrix(0,nr,nc)
  cur.ind <- 1
  clnm <- c()
  nn <- ncol(randomWalks)
  for (i in 1:nrow(rppa2)) {
    genNum <- rppa2$num[i]
    genWalks <- randomWalks[,1] == genNum
    walks <- as.numeric(unlist(randomWalks[genWalks,2:nn]))
    
    tt <- table(intersect(walks,rppa2$num))
    tt <- as.numeric(names(tt))
    neig <- c()
    for (j in 1:length(tt)) {
      if (sum(walks == tt[j]) > 1) {
        neig <- append(neig, tt[j])
      }
    }
    zz <- neig != genNum
    neig<-neig[zz]
    inds <- which(rppa2$num %in% neig)
    for (ii in 1:length(inds)) {
      cname <- paste(rppa$RPPA[i], rppa$RPPA[inds[ii]], sep = "-")
      clnm <- append(clnm, cname)
      temp <- rep(1,nr)
      ind <- which(dd[,i] > dd[,inds[ii]])
      temp[ind] <- -1
      rs[,cur.ind] <- temp
      cur.ind <- cur.ind+1
    }
  }
  colnames(rs) <- clnm
  return(rs)
}

lasso <- function(x,y,above=0) {
  require(glmnet)
  # get rid of NAs
  # keep = indices in x that are not NA
  
  for (i in 1:ncol(x)) {
    keep <- c(1:nrow(x))[!is.na(x[,i])]
    x <- x[keep,]
    #convert the columns to numeric
    x[,i] <- as.numeric(x[,i])
    y <- y[keep]
  }
  
  x <- as.matrix(x)
  
  #cl=parallel::makeCluster(5,type="PSOCK")
  #doParallel::registerDoParallel(cl)
  fit.cv <- cv.glmnet(x=x,y=y,family="cox",alpha=1,standardize=FALSE,nfolds=10,parallel = TRUE) #data have been already standardized 
  #parallel::stopCluster(cl)
  
  lambda <- fit.cv$lambda.min
  final <- glmnet(x=x,y=y,family="cox",alpha=1,lambda=lambda, standardize=FALSE) # glmnet is standardized by default, so turn it off
  coef.fit <- coef(final,s=lambda)[2:(ncol(x)+1)] 	# exclude intercept
  cols.include <- which(abs(coef.fit) > above)
  return(cols.include)
}

#presmat<-function(data,mySurv,tri){
presmat<-function(re){
  data<-re$data
  mySurv<-re$mySurv
  tri<-re$BOO$tr.inds
  
  fmat<-matrix(0,100,dim(data)[2])
  
  for(seed in 1:dim(tri)[1])
  { 
    x.train=data[tri[seed,],]
    y.train=mySurv[tri[seed,],]
    
    fmat[seed,pres(x.train,sum(y.train[,2]),0.8)]<-1
  }
  return(fmat)
  
}


cox.apply <- function(re,coxscr, pofFlag){
  data<-re$data
  mySurv<-re$mySurv
  tri<-re$BOO$tr.inds
  tsi<-re$BOO$ts.inds
  l<-dim(tri)[1]
  if(pofFlag == 1){
    com<-combn(colnames(data),2)
    #com <- colnames(dataa)
    clnm <- paste(com[1,],com[2,],sep='-')
    coefs<-matrix(0, l, 8515)
    #c index results
    cx <-NA*c(1:l)
    for(seed in 1:l){
      xtrain=(data[tri[seed,],])
      xtrain= comparator(xtrain)
      colnames(xtrain) <- clnm
      xtrain<-xtrain[,which(as.numeric(coxscr[seed,])==1)]
      ytrain=mySurv[tri[seed,],]
      
      xtest=(data[tsi[seed,],])
      xtest=comparator(xtest)
      colnames(xtest)=clnm
      xtest<-xtest[,which(as.numeric(coxscr[seed,])==1)]
      ytest=mySurv[tsi[seed,],]
      if(ncol(xtrain) > 5) #useLASSO &  further shrink by LASSO, if only a few features, no need to use LASSO, note 5 is a quite arbitrary setting
      {
        # do LASSO without cross validation to get the features to include in the model
        # change x.train and x.test to only include those features
        cols.include <- c()
        iter <- 0
        while(length(cols.include)<1)
        {
          set.seed(iter+1)
          cols.include <- try(lasso(x=xtrain,y=ytrain,above=0))
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
        xtrain <- xtrain[,cols.include]
        xtest <- xtest[,cols.include]
      }
      
      cox <- coxph(ytrain~xtrain)
      print(paste("seed:",seed,sep = " "))
      cox$coefficients[is.na(cox$coefficients)]=0 # convert NAs to zero if any 
      
      cox.predict <- as.matrix(xtest)%*%cox$coefficients
      coefs[seed,match(colnames(xtrain),clnm)]<-cox$coefficients
      c1<-concordance.index(cox.predict,ytest[,1], ytest[,2]) 
      cx[seed]<-c1$c.index
    }
    
  }
  else{
    clnm<-colnames(data)
    coefs<-matrix(0, l, 131)
    cx <-NA*c(1:l)
    for(seed in 1:l){
      xtrain=(data[tri[seed,],])
      colnames(xtrain) <- clnm
      xtrain<-xtrain[,which(as.numeric(coxscr[seed,])==1)]
      ytrain=mySurv[tri[seed,],]
      
      xtest=(data[tsi[seed,],])
      colnames(xtest)=clnm
      xtest<-xtest[,which(as.numeric(coxscr[seed,])==1)]
      ytest=mySurv[tsi[seed,],]
      
      if(ncol(xtrain) > 5) #useLASSO &  further shrink by LASSO, if only a few features, no need to use LASSO, note 5 is a quite arbitrary setting
      {
        # do LASSO without cross validation to get the features to include in the model
        # change x.train and x.test to only include those features
        cols.include <- c()
        iter <- 0
        while(length(cols.include)<1)
        {
          set.seed(iter+1)
          cols.include <- try(lasso(x=xtrain,y=ytrain,above=0))
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
        xtrain <- xtrain[,cols.include]
        xtest <- xtest[,cols.include]
      }
      
      cox <- coxph(ytrain~xtrain)
      print(paste("seed:",seed,sep = " "))
      cox$coefficients[is.na(cox$coefficients)]=0 # convert NAs to zero if any 
      
      cox.predict <- as.matrix(xtest)%*%cox$coefficients
      coefs[seed,match(colnames(xtrain),clnm)]<-cox$coefficients
      c1<-concordance.index(cox.predict,ytest[,1], ytest[,2]) 
      cx[seed]<-c1$c.index
    }
    
  }
  ree<-list()
  ree$coefs<-coefs
  ree$cinds<-cx
  
  print("---END---")#print(paste(colnames(x.train[,1]),"---END---"))
  
  return(ree)
  
}


applyRSF <- function(data,embs, data.cox,pofFlag){
  data_COXresults<-rsf(data,embs,0,data.cox,pofFlag)
  col.sumModels <- colSums(data.cox)
  col.sumCoefs <- colSums(data_COXresults$coefs)
  imp <- col.sumCoefs/col.sumModels
  na.ind <- is.na(imp)
  imp[na.ind]<-0
  data_POF_COX <- list(cinds = data_COXresults$cinds,coefs=imp,numOfModels = col.sumModels)
  return(data_POF_COX)
}


individualRank <- function(pofRes,rawData){
  xx <- (pofRes$coefs * pofRes$numOfModels)
  data <- data.frame(POF_Feature=names(pofRes$coefs),featureImportance=xx)
  data$featureImportance <- data$featureImportance / max(data$featureImportance)
  ss <- strsplit(as.character(data$POF_Feature), split = '-')
  df <- data.frame(f1=vector(length = nrow(data)),f2=vector(length = nrow(data)),imp=numeric(length = nrow(data)))
  for(i in 1:nrow(data)){
    df$f1[i] <- ss[[i]][1]
    df$f2[i] <- ss[[i]][2]
    df$imp[i] <- data$featureImportance[i]
  }
  
  tot.feats <- unique(c(df$f1,df$f2))
  ind_features <- data.frame(feature=tot.feats, included_POF = vector(length = length(tot.feats)), 
                             importance = numeric(length = length(tot.feats)),
                             totImp = numeric(length = length(tot.feats)),
                             rankDifference = numeric(length = length(tot.feats)),
                             RankPOF = numeric(length = length(tot.feats)),
                             RankIndividual = numeric(length = length(tot.feats)))
  for(i in 1:length(tot.feats)){
    
    ff <- tot.feats[i]
    st <- ""
    im = 0
    counter = 0
    for(j in 1:nrow(data))
    {
      if(ff == ss[[j]][1]){
        im = im + data$featureImportance[j]
        st <- paste(ss[[j]][2],st, sep = ' , ')
        counter = counter + 1
      }
      if(ff == ss[[j]][2]){
        im = im + data$featureImportance[j]
        st <- paste(ss[[j]][1],st,sep = ' , ')
        counter = counter + 1
      }
    }
    ind_features$importance[i] <- im/counter
    ind_features$totImp[i] <- im
    ind_features$included_POF[i] <- st
  }
  ind_features <- ind_features[order((ind_features$importance)),]
  
  rawCoef <- rawData$coefs
  rawCoef <- rawCoef[order(-rawCoef)]
  ind <- which(names(rawCoef) == ind_features$feature[1])
  for (i in 1:nrow(ind_features)) {
    ind <- which(names(rawCoef) == ind_features$feature[i])
    ind_features$rankDifference[i] <- i - ind
    ind_features$RankPOF[i] <- i
    ind_features$RankIndividual[i] <- ind
  }
  ord <- order(ind_features$rankDifference)
  return(ind_features[ord,])
}

filt_rw<-function(prer_cox){
  ft = colnames(prer_cox)
  dd <- matrix(nrow = length(ft),ncol=2)
  dd[,1] <- sapply(strsplit(as.character(ft), split='-', fixed=TRUE), function(x) (x[1]))
  dd[,2] <- sapply(strsplit(as.character(ft), split='-', fixed=TRUE), function(x) (x[2]))
  
  ddd = apply(dd, 1,function(x) (paste(x[2],x[1],sep = "-")))
  
  check = logical(length = length(ft))
  for (i in 2:length(check)) {
    cont = ddd[i]
    check[i] = cont %in% ft[1:(i-1)]
    
  }
  prer_cox[,check] <- 0
  return(prer_cox)
}

