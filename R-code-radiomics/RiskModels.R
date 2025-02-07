options(rgl.useNULL = TRUE)
library(rgl)
library(Rdimtools)
library(stats)
library(boot)
library(ranger)
library(survival)
library(glmnet)
library(mboost)

library(randomForestSRC)

#### Import Libraries ### 
library(survival)
library(purrr)
library(glmnet)
library(dplyr)  
library(ggkm)
library(ggplot2)
library(ggfortify)
'%ni%' <- Negate('%in%')


library(randomForest)
library(randomForestSRC)
library(survival)
library(ranger)

library(praznik)

Spearman_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0) 
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
  c(
    return(cor(x, time, method='s'))
   

  )
}

Pearson_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
  c(
    return(cor(x, time, method='p'))
   

  )
}

Kendall_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
  c(
    return(cor(x, time, method='k'))
   

  )
}

Mifs_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
    res <- do.mifs(as.matrix(x), ndim = ncol(x), time,beta = 1,discretize = "default",preprocess = "null")
    res <- res[1:2]
    res$Y <- colMeans(res$Y)
    res <- as.data.frame(res)
    res <- res[order(res$featidx),]
  c(
    return(res$Y))
}

NJMIM_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
    res <- NJMIM(as.data.frame(x), time, k = ncol(x))
    res <- as.data.frame(res)
    res <- res[order(res$selection),]
  c(
    return(res$score))
}

MRMR_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
    res <- MRMR(as.data.frame(x), time, k = ncol(x))
    res <- as.data.frame(res)
    res <- res[order(res$selection),]
  c(
    return(res$score))
}

JIM_Select <- function(x, indices){
    x <- x[indices,]
    status <- x$status
    drop <- which(status==0)
    x <- x[!seq_len(nrow(x)) %in% drop,]
    time <- x$time
    x[c('time','status')] <- NULL 
    res <- JIM(as.data.frame(x), time, k = ncol(x))
    res <- as.data.frame(res)
    res <- res[order(res$selection),]
  c(
    return(res$score))
}

RFVHimp_Select <- function(x, indices){
    x <- x[indices,]
    x_obj <- rfsrc(Surv(time, status) ~ ., data = as.data.frame(x), bootstap='none', importance = 'permute')
    res <- var.select(Surv(time, status) ~ ., object = x_obj, method='vh.vimp', verbose = FALSE)
    out <- res$varselect 
    res <- out[ order(row.names(out)), ]
    c(
        return(res))
}

RFVH_Select <- function(x, indices){
    x <- x[indices,]
    x_obj <- rfsrc(Surv(time, status) ~ ., data = as.data.frame(x), bootstap='none', importance = 'permute')
    res <- var.select(Surv(time, status) ~ ., object = x_obj, method='vh', verbose = FALSE)
    out <- res$varselect 
    res <- out[ order(row.names(out)), ]
    c(
        return(res[,1]))
}

RFMD_Select <- function(x, indices){
    x <- x[indices,]
    x_obj <- rfsrc(Surv(time, status) ~ ., data = as.data.frame(x), bootstap='none', importance = 'permute')
    res <- var.select(Surv(time, status) ~ ., object = x_obj, method='md', verbose=FALSE)
    out <- res$varselect 
    res <- out[ order(row.names(out)), ]
    c(
        return(res[,1]))
}

PVIRF_Select <- function(x, indices){
    x <- x[indices,]
    x_ob <- ranger(Surv(time,status) ~., data=as.data.frame(x), importance = "permutation",seed = 1)
    res <- importance_pvalues(x_ob,method = c("janitza"),num.permutations = 100,formula = NULL,data = NULL,)
    c(
        return(res[,2]))
}

IMPRF_Select <- function(x, indices){
    x <- x[indices,]
    x_ob <- ranger(Surv(time,status) ~., data=as.data.frame(x), importance = "impurity_corrected",seed = 1)
    res <- importance_pvalues(x_ob,method = c("janitza"),num.permutations = 100,formula = NULL,data = NULL,)
    c(
        return(res[,2]))
}

Cox_Select <- function(df, indices){
    df <- df[indices,]
    df2 <- df 
    df2$time <- NULL 
    df2$status <- NULL
    univ_formulas <- sapply(c(colnames(df2)),function(x)as.formula(paste('Surv(time,status)~',x)))
    #making a list of models
    univ_models <- lapply(univ_formulas, function(x){coxph(x,data=df)})

    #extract data (here I've gone for HR and confint)
    results <- lapply(univ_models,function(x){cbind(Pval=coef(summary(x))[5])})
    results <- as.data.frame(results)
    return(as.matrix(results))}


Cox_Net <- function(xi, y, n, K){ 
    ##hpyerparamaters
   
    #bootstrap:
    #number of bootstrap samples # size of bootstaps
    bootpred <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- NULL
    k.out2 = 0
    bootpred2 <- NULL 
      for (k in K){
        x <- xi[1:k]
        bootpred <- list()
        concordance.tmp2 <- NULL
        for (i in seq_len(n)) {
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            cv <- cv.glmnet(as.matrix(bootdat), y.boot, family="cox",alpha=0.5,  standardize=TRUE)
            bootmod <- glmnet(as.matrix(bootdat), y.boot, family="cox",alpha=0.5, lambda=cv$lambda.min,  standardize=TRUE)
            bootpred[[i]] <- bootmod
            concordance.tmp2 <- append(concordance.tmp2, (min(cv$cvm)))}
        if (is.null(concordance.tmp3)){ concordance.tmp3 = concordance.tmp2}
        if (mean(concordance.tmp2) <= mean(concordance.tmp3)){
            concordance.tmp3 = concordance.tmp2
            bootpred2 = bootpred
            k.out2 = k}}
    out = list(k.out2, bootpred2)
return(bootpred2)}

Cox_Lasso <- function(xi, y, n, K){ 
    ##hpyerparamaters
   
    #bootstrap:
    #number of bootstrap samples # size of bootstaps
    bootpred <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- NULL
    k.out2 = 0
    bootpred2 <- NULL 
      for (k in K){
        x <- xi[1:k]
        bootpred <- list()
        concordance.tmp2 <- NULL
        for (i in seq_len(n)) {
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            cv <- cv.glmnet(as.matrix(bootdat), y.boot, family="cox",alpha=1,  standardize=TRUE)
            bootmod <- glmnet(as.matrix(bootdat), y.boot, family="cox",alpha=1, lambda=cv$lambda.min,  standardize=TRUE)
            bootpred[[i]] <- bootmod
            concordance.tmp2 <- append(concordance.tmp2, (min(cv$cvm)))}
        if (is.null(concordance.tmp3)){ concordance.tmp3 = concordance.tmp2}
        if (mean(concordance.tmp2) <= mean(concordance.tmp3)){
            concordance.tmp3 = concordance.tmp2
            bootpred2 = bootpred
            k.out2 = k}}
    out = list(k.out2, bootpred2)
return(bootpred2)}

### 18/09/10
COX_BOOT <- function(xi, y, n, K){ 
    #bootstrap:
    #number of bootstrap samples # size of bootstaps
    bootpred <- list()
    set.seed(42) #for reproducibility
    #loop over n
    bootmodout <- NULL
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            set.seed(i) 
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            ties=c("efron","breslow","exact")
            concordance.tmp = 0
            concordance.tmp2 <- NULL
            bootmod = list()
            count =1
            for (param in ties){
                for (j in seq_len(n)){
                    smp_size <- smp_size <- floor(0.8 * nrow(bootdat))
                    set.seed(j) 
                    train_ind <- sample(seq_len(nrow(bootdat)), size = smp_size)
                    train <- bootdat[train_ind, ]
                    intval <- bootdat[-train_ind, ]
                    trainY <- y.boot[train_ind, ]
                    Yval <- y.boot[-train_ind, ]
                    bootmod[[count]] <- try(coxph(trainY~.,data =train, ties = param),silent=T)
                    tmpK <- predict(bootmod[[j]], intval)
                    tmpK[!is.finite(tmpK)] <- 0
                    tmp <- coxph(Yval~.,data =as.data.frame(tmpK))              
                    concordance.tmp2 <- append(concordance.tmp2, summary(tmp)$concordance[1])
                    count = count +1}
                    #print(bootmod)
                if (mean(concordance.tmp2) > mean(concordance.tmp)){
                    concordance.tmp = concordance.tmp2
                    bootmodout = bootmod[[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]]}}
            bootpred[[i]] <- bootmodout}
         if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out2 = k}}
        out = list(k.out2, bootpred2)
return(bootpred2)}







### 18/09/10
SURV_REG <- function(xi, y,n, K){ 
    #bootstrap:
    #number of bootstrap samples # size of bootstaps
    bootpred <- list()
    set.seed(42) #for reproducibility
    #loop over n
    bootmodout <- NULL
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            set.seed(i) 
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            concordance.tmp = 0
            dists=c("weibull", "exponential", "gaussian", "logistic","lognormal", "loglogistic")
            bootmod = list()
            for (param in dists){
                concordance.tmp2 <- NULL
                for (j in seq_len(10)){
                    smp_size <- smp_size <- floor(0.8 * nrow(bootdat))
                    set.seed(j) 
                    train_ind <- sample(seq_len(nrow(bootdat)), size = smp_size)
                    train <- bootdat[train_ind, ]
                    intval <- bootdat[-train_ind, ]
                    trainY <- y.boot[train_ind, ]
                    Yval <- y.boot[-train_ind, ]
                    bootmod[[j]] <-survreg(trainY~.,data =train, dist = param)
                    tmpK <- predict(bootmod[[j]], intval)
                    tmpK[!is.finite(tmpK)] <- 0
                    tmp <- coxph(Yval~.,data =as.data.frame(tmpK))              
                    concordance.tmp2 <- append(concordance.tmp2, summary(tmp)$concordance[1])}
            if (mean(concordance.tmp2) > mean(concordance.tmp)){
                concordance.tmp = concordance.tmp2
                bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}
                    bootpred[[i]] <- bootmodout}}
     if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out2 = k}}
        out = list(k.out2, bootpred2)
return(bootpred2)}



BGLM_Cindex<- function(xi, y,n, K){ 
    #bootstrap:
    #number of bootstrap samples defined globally with parameter 'n' # size
    bootpred <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        concordance.tmp3 <- 0
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- glmboost(as.matrix(bootdat[1:ncol(bootdat)]), y.boot, family=Cindex(), control=boost_control(mstop = 200, nu=0.1), center=FALSE)
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- glmboost(y.boot~.,  data =bootdat, family = Cindex(),control = boost_control(mstop = mstop(cvm)), center=FALSE)
            bootpred[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
            concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out2 = k}}
out <- list(k.out2, bootpred2)
return(bootpred2)}



BGLM_Cox <- function(xi, y,n, K){ 
    #bootstrap:
    #number of bootstrap samples defined globally with parameter 'n' # size
    bootpred <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- 0
    for (k in K){
        xi<- as.data.frame(xi)
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- glmboost(as.matrix(bootdat[1:ncol(bootdat)]), y.boot, family=CoxPH(), control=boost_control(mstop = 200, nu=0.1), center=FALSE)
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- glmboost(y.boot~.,  data =bootdat, family = CoxPH(),control = boost_control(mstop = mstop(cvm)), center=FALSE)
            bootpred[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
            concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out = k}}
out <- list(k.out, bootpred2)
return(bootpred2)}


BGLM_Weibull <- function(xi, y,n, K){ 
    #bootstrap:
    #number of bootstrap samples defined globally with parameter 'n' # size
    bootpred <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- glmboost(as.matrix(bootdat[1:ncol(bootdat)]), y.boot, family=Weibull(), control=boost_control(mstop = 200, nu=0.1), center=FALSE)
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- glmboost(y.boot~.,  data =bootdat, family = Weibull(),control = boost_control(mstop = mstop(cvm)), center=FALSE)
            bootpred[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
            concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out = k}}
out <- list(k.out, bootpred2)
return(bootpred2)}


BT_Weibull <- function(xi, y,n, K){ 
    #bootstrap:
    bootout <- list()
    set.seed(42) #for reproducibility
    #loop over n
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- blackboost(y.boot~.,  data =bootdat, family =Weibull() , control = boost_control())
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- blackboost(y.boot~.,  data =bootdat, family = Weibull(), control = boost_control(mstop = mstop(cvm)))
            bootout[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
        concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootout
                k.out2 = k}}
out <- list(k.out2, bootpred2)
return(bootpred2)}

BT_Cox <- function(xi, y,n, K){ 
    #bootstrap:
    bootout <- list()
    set.seed(42) #for reproducibility
    #loop over n
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        concordance.tmp3 <- 0
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- blackboost(y.boot~.,  data =bootdat, family =CoxPH() , control = boost_control())
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- blackboost(y.boot~.,  data =bootdat, family = CoxPH(), control = boost_control(mstop = mstop(cvm)))
            bootout[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
        concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootout
                k.out2 = k}}
out <- list(k.out2, bootpred2)
return(bootpred2)}

BT_Cindex <- function(xi, y,n, K){ 
    #bootstrap:
    bootout <- list()
    set.seed(42) #for reproducibility
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            concordance.tmp2 <- NULL
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            model <- blackboost(y.boot~.,  data =bootdat, family =Cindex() , control = boost_control())
            cv10f <- cv(model.weights(model), type = "kfold")
            cvm <- cvrisk(model, folds = cv10f, papply = lapply)
            bootmod <- blackboost(y.boot~.,  data =bootdat, family = Cindex(), control = boost_control(mstop = mstop(cvm)))
            bootout[[i]] <- bootmod #calculate predictions from this model
            tmpK <- predict(bootmod, bootdat)
            tmpK[!is.finite(tmpK)] <- 0
            tmp <- coxph(y.boot~.,data =as.data.frame(tmpK))
            concordance.tmp2[i] <- summary(tmp)$concordance[1]}
        concordance.tmp2[!is.finite(concordance.tmp2)] <- 0
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootout
                k.out2 = k}}
out <- list(k.out2, bootpred2)
return(bootpred2)}

RSF <- function(xi, y,n, K){ 
    #bootstrap:
    n <- n #number of bootstrap samples # size of bootstaps
    bootout <- list()
    bootpred <- list()
    set.seed(42) #for reproducibility
    #loop over n
    count =0 
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            minnodes <- c(5,10,15,20)
            maxDepth <- c(10,15,20)
            numtrees <- c(100,200,1000,2000)
            concordance.tmp = 0 
            bootmod <- list()
            concordance.tmp2 <- NULL
                    for (node in minnodes){
                        for (depth in maxDepth){
                            for (ntrees in numtrees){
                                for (j in seq_len(10)){
                                    smp_size <- smp_size <- floor(0.8 * nrow(bootdat))
                                    set.seed(j) 
                                    train_ind <- sample(seq_len(nrow(bootdat)), size = smp_size)
                                    train <- bootdat[train_ind, ]
                                    intval <- bootdat[-train_ind, ]
                                    trainY <- y.boot[train_ind, ]
                                    Yval <- y.boot[-train_ind, ]
                                    bootmod[[j]] <-ranger(trainY ~ ., data = train, num.trees = ntrees , min.node.size=node, max.depth=depth, minprop=0.1,splitrule = "logrank", importance = "permutation", local.importance = TRUE)
                                    tmpK <- predict(bootmod[[j]], intval)
                                    tmpK <- rowMeans(tmpK$chf)
                                    tmpK[!is.finite(tmpK)] <- 0
                                    tmp <- coxph(Yval~.,data =as.data.frame(tmpK))
                                    concordance.tmp2[j] <- summary(tmp)$concordance[1]}
                                if (mean(concordance.tmp2) > mean(concordance.tmp)){
                                    concordance.tmp = concordance.tmp2
                                    bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
                            if (mean(concordance.tmp2) > mean(concordance.tmp)){
                                concordance.tmp = concordance.tmp2
                                bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
                        if (mean(concordance.tmp2) > mean(concordance.tmp)){
                            concordance.tmp = concordance.tmp2
                            bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
            bootpred[[i]] <- bootmodout}
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out2 = k}}
        out = list(k.out2, bootpred2)
return(bootpred2)}

ET_RF <- function(xi, y,n, K){ 
    #bootstrap:
    n <- n #number of bootstrap samples # size of bootstaps
    bootout <- list()
    bootpred <- list()
    set.seed(42) #for reproducibility
    #loop over n
    count =0
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)) {
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            minnodes <- c(5,10,15,20)
            maxDepth <- c(10,15,20)
            numtrees <- c(100,200,1000,2000)
            concordance.tmp = 0 
            bootmod <- list()
            concordance.tmp2 <- NULL
                    for (node in minnodes){
                        for (depth in maxDepth){
                            for (ntrees in numtrees){
                                for (j in seq_len(10)){
                                    smp_size <- smp_size <- floor(0.8 * nrow(bootdat))
                                    set.seed(j) 
                                    train_ind <- sample(seq_len(nrow(bootdat)), size = smp_size)
                                    train <- bootdat[train_ind, ]
                                    intval <- bootdat[-train_ind, ]
                                    trainY <- y.boot[train_ind, ]
                                    Yval <- y.boot[-train_ind, ]
                                    bootmod[[j]] <-ranger(trainY ~ ., data = train, num.trees = ntrees , min.node.size=node, max.depth=depth, minprop=0.1,splitrule = "extratrees")
                                    tmpK <- predict(bootmod[[j]], intval)
                                    tmpK <- rowMeans(tmpK$chf)
                                    tmpK[!is.finite(tmpK)] <- 0
                                    tmp <- coxph(Yval~.,data =as.data.frame(tmpK))
                                    concordance.tmp2[j] <- summary(tmp)$concordance[1]}
                                if (mean(concordance.tmp2) > mean(concordance.tmp)){
                                    concordance.tmp = concordance.tmp2
                                    bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
                            if (mean(concordance.tmp2) > mean(concordance.tmp)){
                                concordance.tmp = concordance.tmp2
                                bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
                        if (mean(concordance.tmp2) > mean(concordance.tmp)){
                            concordance.tmp = concordance.tmp2
                            bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]}}
            bootpred[[i]] <- bootmodout}
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
                concordance.tmp3 = concordance.tmp2
                bootpred2 = bootpred
                k.out2 = k}}
        out = list(k.out2, bootpred2)
return(bootpred2)}

MSR_RF <- function(xi, y,n, K){ 
    #bootstrap:
    K=K
    n <- n #number of bootstrap samples # size of bootstaps
    bootout <- list()
    bootpred <- list()
    set.seed(42) #for reproducibility
    #loop over n
    count =0 
    concordance.tmp3 <- 0
    for (k in K){
        x <- xi[1:k]
        bootmod <- list()
        for (i in seq_len(n)){
            set.seed(i)
            sampled <- sample(nrow(x), replace = TRUE)
            bootdat <- x[sampled,] #bootstrap resample of data
            y.boot <- y[sampled]
            minnodes <- c(5,10,15,20)
            maxDepth <- c(10,15,20)
            numtrees <- c(100,200,1000,2000)
            concordance.tmp = 0 
            bootmod <- list()
            concordance.tmp2 <- NULL
                for (node in minnodes){
                    for (depth in maxDepth){
                        for (ntrees in numtrees){
                            for (j in seq_len(10)){
                                smp_size <- smp_size <- floor(0.8 * nrow(bootdat))
                                set.seed(j) 
                                train_ind <- sample(seq_len(nrow(bootdat)), size = smp_size)
                                train <- bootdat[train_ind, ]
                                intval <- bootdat[-train_ind, ]
                                trainY <- y.boot[train_ind, ]
                                Yval <- y.boot[-train_ind, ]
                                bootmod[[j]] <-ranger(trainY ~ ., data = train, num.trees = ntrees , min.node.size=node, max.depth=depth, minprop=0.1,splitrule = "maxstat")
                                tmpK <- predict(bootmod[[j]], intval)
                                tmpK <- rowMeans(tmpK$chf)
                                tmpK[!is.finite(tmpK)] <- 0
                                tmp <- coxph(Yval~.,data =as.data.frame(tmpK))
                                concordance.tmp2[j] <- summary(tmp)$concordance[1]}
                            if (mean(concordance.tmp2) > mean(concordance.tmp)){
                                concordance.tmp = concordance.tmp2
                                bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]
                                k.out =k}}
                        if (mean(concordance.tmp2) > mean(concordance.tmp)){
                            concordance.tmp = concordance.tmp2
                            bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]
                            k.out = k }}
                    if (mean(concordance.tmp2) > mean(concordance.tmp)){
                        concordance.tmp = concordance.tmp2
                        bootmodout = bootmod[which.min(abs(concordance.tmp2-mean(concordance.tmp2)))]
                        k.out=k}}
            bootpred[[i]] <- bootmodout}
        if (mean(concordance.tmp2) > mean(concordance.tmp3)){
            concordance.tmp3 = concordance.tmp2
            bootpred2 = bootpred
            k.out2 = k.out}}
    out = list(k.out2, bootpred2)
return(bootpred2)}


Predict.Risk <- function(test,df){
res.out = list()
for (i in seq_along(test)){
        Model.out <- list()
        count2 =1 
        for (j in test[i]){
            Model<- list()
            count =1
            for (k in j) {
                res <- NULL
                pred.out <- NULL
                for (z in k) {
                    if (class(z)[1] == 'coxnet'){
                        df2 <- df[unlist(z$beta@Dimnames[1])]
                        res <- predict(z, as.matrix(df2))} 
                    if (class(z[1])[1] == 'ranger') {
                        res <- predict(z[1], df)
                        res <- rowMeans(res$chf)} 
                    if (class(z[[1]])[1] == 'ranger') {
                        res <- predict(z[[1]], df)
                        res <- rowMeans(res$chf)} 
                    if (class(z)[1] == 'glmboost' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    if (class(z)[1] == 'blackboost' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    if (class(z)[1] == 'coxph' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    pred.out <- cbind(pred.out, res)}
            Model[[count]] <- rowMeans(pred.out)  # average output of bootrstraps
            count =count +1 } # for k in j-- model bootstraps per a model
            Model.out[[count2]] <- Model # for j in i
            count2 = count2 +1 }
    res.out[[i]] <- Model.out} 
return(res.out)}

Predict.Risk.ALL <- function(test,df){
res.out = list()
for (i in seq_along(test)){
        Model.out <- list()
        count2 =1 
        for (j in test[i]){
            Model<- list()
            count =1
            for (k in j) {
                res <- NULL
                pred.out <- NULL
                for (z in k) {
                    if (class(z)[1] == 'coxnet'){
                        df2 <- df[unlist(z$beta@Dimnames[1])]
                        res <- predict(z, as.matrix(df2))} 
                    if (class(z[1])[1] == 'ranger') {
                        res <- predict(z[1], df)
                        res <- rowMeans(res$chf)} 
                    if (class(z[[1]])[1] == 'ranger') {
                        res <- predict(z[[1]], df)
                        res <- rowMeans(res$chf)} 
                    if (class(z)[1] == 'glmboost' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    if (class(z)[1] == 'blackboost' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    if (class(z)[1] == 'coxph' & (class(z)[1] != 'ranger')){
                        res <- predict(z, df)}
                    pred.out <- cbind(pred.out, res)}
            Model[[count]] <- pred.out  # average output of bootrstraps
            count =count +1 } # for k in j-- model bootstraps per a model
            Model.out[[count2]] <- Model # for j in i
            count2 = count2 +1 }
    res.out[[i]] <- Model.out} 
return(res.out)}

Cindex.Risk <- function(pred, y){
    OUT <- list()
    count3 <- 1 
    for (i in pred){
        count2 <- 1
        Res.out <- list()
        for (j in i){
            count<-1
            Res <- NULL
            for (k in j){
                k[!is.finite(k)] <- 0
                tmp <- coxph(y~.,data =as.data.frame(k))
                Res[count] <- summary(tmp)$concordance[1]
                count <- count +1
            }
            Res.out[[count2]] <- Res 
            count2 <- count2 + 1
        }
        OUT[count3] <- Res.out
        count3 <- count3 +1
    }
return(OUT)}

Risk_R <- function(train, K, n){
    K=K
    n=n
    y <- Surv(train$time, train$status)
    train2 <- train
    train$time <- NULL 
    train$status <- NULL
    ### Order the data frames by colnames and remove outcome vars
    train[ , order(names(train))]
    head(train)
    string_vars <- colnames(train)
    string_vars <- string_vars[string_vars != "time"]
    string_vars <- string_vars[string_vars != "status"]
    ## Creat tmp, colnames, order the n extrapolate colnames
    tmp <- boot(train2,Pearson_Select, R=n,) ## starts i with pearson, then spearman, and so on... 
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select <- NULL 
    feature_select$pearson <- tmp
    tmp <- boot(train2,Spearman_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Spearman <- tmp
    tmp <- boot(train2,Kendall_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Kendall <- tmp
    tmp <- boot(train2,Mifs_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Mifs <- tmp
    tmp <- boot(train2,JIM_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$JIM <- tmp
    tmp <- boot(train2,MRMR_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$MRMR <- tmp
    tmp <- boot(train2,RFVHimp_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFVHimp_Select <- tmp
    tmp <- boot(train2,RFVH_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFVH_Select <- tmp
    tmp <- boot(train2,RFMD_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFMD_Select <- tmp
    tmp <- boot(train2,IMPRF_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$IMPRF_Select <- tmp
    tmp <- boot(train2,PVIRF_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$PVIRF <- tmp
    tmp <- boot(train2,Cox_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Cox <- tmp
    set.seed(47) ## for reproducibility of random features
    feature_select$random <- slice(tmp, sample(1:n()))
    set.seed(NULL)
    count = 1 
    count2 =1 
    out = list()
    res.out = list()
    res.out2 = list()
    ### Models...
    for (i in feature_select){
        res <- BGLM_Cox(train[rownames(i)], y, n, K)
        res.out[["BGLM_CoxPH"]] <- res 
        res <- BGLM_Cindex(train[rownames(i)], y, n, K)
        res.out[["BGLM_Cindex"]] <- res 
        res <- BGLM_Weibull(train[rownames(i)], y, n, K)
        res.out[["BGLM_Weibull"]] <- res 
        res <- BT_Cox(train[rownames(i)], y, n, K)
        res.out[["BT_CoxPH"]] <- res 
        res <- BT_Cindex(train[rownames(i)], y, n, K)
        res.out[["BT_Cindex"]] <- res 
        res <- BT_Weibull(train[rownames(i)], y, n, K)
        res.out[["BT_Weibull"]] <- res 
        res <- Cox_Lasso(train[rownames(i)], y, n, K)
        res.out[["Cox_Lasso"]] <- res 
        res <- Cox_Net(train[rownames(i)], y, n, K)
        res.out[["Cox_Net"]] <- res 
        res <- COX_BOOT(train[rownames(i)], y, n, K)
        res.out[["Cox"]] <- res 
        res <- RSF(train[rownames(i)], y, n, K)
        res.out[["RSF"]] <- res 
        res <- MSR_RF(train[rownames(i)], y, n, K)
        res.out[["MSR_RF"]] <- res 
        res <- ET_RF(train[rownames(i)], y, n, K)
        res.out[["ET_RF"]] <- out 
        res.out2[[i]]<- res.out}
    return(res.out2)
}

Risk_R <- function(train, K, n){
    K=K
    n=n
    y <- Surv(train$time, train$status)
    train2 <- train
    train$time <- NULL 
    train$status <- NULL
    ### Order the data frames by colnames and remove outcome vars
    train[ , order(names(train))]
    head(train)
    string_vars <- colnames(train)
    string_vars <- string_vars[string_vars != "time"]
    string_vars <- string_vars[string_vars != "status"]
    ## Creat tmp, colnames, order the n extrapolate colnames
    tmp <- boot(train2,Pearson_Select, R=n,) ## starts i with pearson, then spearman, and so on... 
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select <- NULL 
    feature_select$pearson <- tmp
    tmp <- boot(train2,Spearman_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Spearman <- tmp
    tmp <- boot(train2,Kendall_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Kendall <- tmp
    tmp <- boot(train2,Mifs_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Mifs <- tmp
    tmp <- boot(train2,JIM_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$JIM <- tmp
    tmp <- boot(train2,MRMR_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$MRMR <- tmp
    tmp <- boot(train2,RFVHimp_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFVHimp_Select <- tmp
    tmp <- boot(train2,RFVH_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFVH_Select <- tmp
    tmp <- boot(train2,RFMD_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$RFMD_Select <- tmp
    tmp <- boot(train2,IMPRF_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$IMPRF_Select <- tmp
    tmp <- boot(train2,PVIRF_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$PVIRF <- tmp
    tmp <- boot(train2,Cox_Select, R=n,)
    tmp <- as.data.frame(colMeans(tmp$t))
    rownames(tmp) <- string_vars
    tmp <- abs(tmp)
    tmp <- tmp[order(-tmp[1]), , drop = FALSE]
    feature_select$Cox <- tmp
    set.seed(47) ## for reproducibility of random features
    feature_select$random <- slice(tmp, sample(1:n()))
    set.seed(NULL)
    count = 1 
    count2 =1 
    out = list()
    res.out = list()
    res.out2 = list()
    ### Models...
     ### Models...
    for (i in feature_select){
        res <- BGLM_Cox(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BGLM_CoxPH"]] <- out 
    out <- NULL
    count=1
     for (i in feature_select){
        res <- BGLM_Cindex(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BGLM_Cindex"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- BGLM_Weibull(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BGLM_Weibull"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- BT_Cox(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BT_CoxPH"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- BT_Cindex(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BT_Cindex"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- BT_Weibull(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["BT_Weibull"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- Cox_Lasso(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["Cox_Lasso"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- Cox_Net(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["Cox_Net"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- COX_BOOT(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["Cox"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- RSF(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["RSF"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- MSR_RF(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["MSR_RF"]] <- out 
    out <- NULL
    count=1
    for (i in feature_select){
        res <- ET_RF(train[rownames(i)], y, n, K)
        out[[count]] <- res
        count = count +1}
    res.out[["ET_RF"]] <- out 
    return(res.out)
}
