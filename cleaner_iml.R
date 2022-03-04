#Imports
{
  defaultW <- getOption("warn") 
  options(warn = -1) 
  
  library(MASS)
  library(caret)
  library(data.table)
  library(doParallel)
  library(rattle)
  library(plyr)
  library(testthat)
  library(checkmate)
  library(ALEPlot)
  library(segmented)
  library(xtable)
  library(quantmod)
  library(lubridate)
  library(GGally)
  library(plm)
  library(parallel)
  library(plotly)
  library(numDeriv)
  
  options(warn = defaultW)
}

{
  #clear env
  rm(list = ls())
  
  #graph output dir
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  output.dir <- getwd()
  fig.dir <- paste(output.dir,"/Figures/",sep = "")
  
  
  #make dirs if needed
  for(m in model_list){
    wd = paste(fig.dir, "/", m, sep = "")
    if ( !(dir.exists(wd)) ) {
      dir.create(wd)
    }else{
      do.call(file.remove, list(list.files(wd, full.names = TRUE)))
    }
  }
  
}



gen_data <- function()
{
  
}

# function to compute number of features
feature_used <- function(pred, feature, sample_size)
{
  dat <- pred$trainingData
  fvalues = dat[,feature] 
  # permute feature
  dat2 = dat[sample(1:nrow(dat), size = sample_size, replace = TRUE),]
  prediction1 = predict(pred,dat2)
  
  sampled_fvalues = sapply(dat2[,feature], function(x){
    sample(setdiff(fvalues, x), size = 1)
  })
  
  dat2 <- data.table(dat2)
  dat2 = dat2[, (feature) := sampled_fvalues]
  dat2$.outcome <- NULL
  prediction2 = predict(pred,dat2)
  
  plot(prediction1~prediction2)
  
  if (any(( prediction1 - prediction2) != 0)) 
    return(TRUE)
  FALSE
}


fit_model <- function(formula, model, train, test, control)
{
  cat("Fitting ", model, "\n")
  
  train_model <- train(formula, data = train, method = model,
                       trControl = control,
                       # preProcess = c("center", "scale"),
                       tuneLength = 10)
  
  y_hat_i <- predict(train_model, test[,x_var])
  y_true <- test[,y_var]
  
  #plot(y_true~y_hat_i)
  #summary(lm(y_true~y_hat_i))$adj
  
  mse = mean((test[,"Y"] - y_hat_i)^2)
  
  return(train_model, mse, y_hat_i)
}


fit_ale <- function(model, train, test, x_var, y_hat_i)
{
  f_0 <- mean(y_hat_i)
  f_ale_1st  <- f_0
  
  
  
  yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata, type = "raw"))
  
  for (v in x_var) {
    ## Do ALE
    j <- which(v == x_var)
    ale_j <- ALEPlot(train, model, pred.fun = yhat, J = j,
                     K = 10^2, NA.plot = TRUE)
    
    loess_j <- loess(ale_j$f.values ~ ale_j$x.values, span = 0.1)
    
    f_j_ale <- predict(loess_j, train[,v]) # based on the ale, we compute the approximation
    
    # f_j_ale <- f_j_ale - mean(f_j_ale)
    f_ale_1st <- f_ale_1st + f_j_ale
    
    {
      file.i <- paste0(fig.dir, model_i,"/ALE","_",j,".pdf")
      pdf(file.i)
      plot(ale_j$f.values~ale_j$x.values,
           ylab = "ALE", xlab = expression(x[j]),
           pch = 20, cex = 0.5)
      lines(predict(loess_j, ale_j$x.values)~ale_j$x.values, col = 2)
      grid(10)
      dev.off()
    }
  }
  #we need to return all the of the ales for each variable from this function
  
}


calc_mec <- function(ale, tol, )
{
  epsilon <- 0.05 # set tolerance
  
  
  MEC_try_function <- function(ale_j) {
    x <- ale_j$x.values
    y <- ale_j$f.values
    K <- 1
    
    lm_j <- lm(y ~ x)
    Rsq_j <- summary(lm_j)$adj
    seg_fit_j <- fitted(lm_j)
    
    ds_plot <- data.frame(seg_fit_j,y,x)
    ds_plot <- ds_plot[order(ds_plot$x),]
    
    {
    file.i <- paste0(fig.dir,model_i,"/MEC","_",j,".pdf")
    pdf(file.i)
    plot(seg_fit_j~x,pch = 20,type = "l", data = ds_plot,
         ylab = "ALE", xlab = expression(x[j]), lwd = 2)
    lines(y ~ x, data = ds_plot,col = 2)
    grid(10)
    dev.off()
    }
    
    
    #why does this need a try catch??
    while(Rsq_j < 1 - epsilon) {
      cat("This is ",K,"\n")
      seg_j <- segmented(lm_j, seg.Z = ~ x, npsi = K)
      seg_fit_j <- fitted(seg_j)
      Rsq_j <- summary(lm(y~seg_fit_j))$adj
      K <- K + 1
      
      ds_plot <- data.frame(seg_fit_j,y,x)
      ds_plot <- ds_plot[order(ds_plot$x),]
      
      {
      pdf(file.i)
      plot(seg_fit_j~x,pch = 20,type = "l", data = ds_plot,
           ylab = "ALE", xlab = expression(x[j]), lwd = 2)
      lines(y ~ x, data = ds_plot,col = 2)
      grid(10)
      dev.off()
      }
    }
    
    return(K)
    
  }
  
  catch_error <- try(MEC_try_function(ale_j),silent  = T)
  try_i <- 1
  
  while(inherits(catch_error,"try-error")) {
    try_i <- try_i + 1
    cat("This is error trial ",try_i,"\n" )
    
    catch_error <- try(MEC_try_function(ale_j), silent  = T)
  }
  
  
  MEC_j <- catch_error
  
  
  MEC <- c(MEC,MEC_j)
  
  V_j <- var(f_j_ale)
  V <- c(V,V_j)
}


main <- function()
{
  cl <- makePSOCKcluster(detectCores())
  registerDoParallel(cl)
  
  #model fitting params
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, allowParallel = T)
  model_list <- c("lm","svmRadial")
  
  #make data
  
  
  #train test split
  ds <- raw_df
  ds_train <- ds
  ds_test <- ds_train
  
  #model formula
  y_var <- "Y"
  x_var <- names(ds)[!names(ds) %in% y_var]
  model_formula <- formula(paste(y_var, " ~ " ,paste(x_var,collapse = " + ")))
  
  #other vars
  MSE_MEC_prop_seq <- MSE_seq <- NF_seq <- IAS_seq <- MEC_model_seq <- train_model_list <- c()
  
  
  for (model_i in model_list) {
    tmp = fit_model(model_formula, model_i, ds_train, ds_test, x_var, y_var, trctrl)
    # model = tmp[1]
    # mse = tmp[2]
    # y_hat_i = tmp[3]
    #train_model_list <- c(train_model_list, list(tmp[1]))
    
    MSE_seq <- c(MSE_seq, tmp[2])
    
    NF <- sapply(x_var, function(x) feature_used(train_model, x, 500) )
    NF_seq <- c(NF_seq, sum(NF))
  
    #doale
    fit_ale(tmp[1], ds_train, ds_test, x_var, tmp[3])
  
    #do mec
  
  }
  
  stopCluster(cl)
  registerDoSEQ()
}


