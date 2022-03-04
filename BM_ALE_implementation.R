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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
output.dir <- getwd()
fig.dir <- paste(output.dir,"/Figures/",sep = "")

simGBM <- function(S = 100, DTintervals = 252, mu = 0.1, sig = 0.2, numSims = 1){
  "
  Generates a Geometric Brownian Motion
  
  Inputs:
  S = starting value
  DTintervals = the number of intervals to generate the bm
  mu = mu variable
  sig = sigma variable
  numSims = number of simulations to run
  
  Returns:
  Matrix
  
the first row is just the starting val
each col is the individual simulation
  "
  dt <- 1/DTintervals#not sure if this is the correct way to do this
  gbmVals <- matrix( rnorm(DTintervals*numSims, dt*(mu - 0.5*sig^2), sig*sqrt(dt)), nrow = DTintervals, ncol = numSims) 
  gbmVals <- apply(gbmVals, 2, cumsum)
  sims <- (matrix(S, nrow = DTintervals, ncol = numSims)) *exp(gbmVals)
  #add the row of initial values
  sims <- rbind(c(rep(S, numSims)),sims )
  return(sims)
}

result = simGBM(numSims = 3)


Sig = matrix(c(0.9, 0.8, 0.8, 0.8), nrow = 2)
sims = mvrnorm(n=100, mu = c(0.14, 0.05), Sigma = Sig)

cov(sims)

eps = rnorm(100, sd = 0.01)

y = 2*sims[,1] + 2*sims[,2]*sims[,1] +  eps

lm(y~sims[,1] + sims[,1]*sims[,2])

# function to compute number of features
feature_used <- function(pred, feature, sample_size){
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


# fit model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3,allowParallel = T)
# model_list <- c("lasso","glmnet","svmLinear","svmRadial","rf","rpart","rpart2")
model_list <- c("lm","svmRadial")

ds <- ds[,!names(ds) %in% names(Filter(is.factor,ds))]


X1 <- sims[,1]
X2 <- sims[,2]
eps <- rnorm(100, sd = 0.1)

non_linear_model <- T
beta_vec <- c(3,-5,3)
Y <- 1 + beta_vec[1]*X1  + beta_vec[2]*X2 + eps
if(non_linear_model) {
  Y <- 1 + beta_vec[1]*X1  +  beta_vec[2] *X2 + beta_vec[3] *(X2^2)  + eps
}

ds <- data.frame(X1,X2,Y)

ds_train <- ds
ds_test <- ds_train
y_var <- "Y"
x_var <- names(ds)[!names(ds) %in% y_var]
model_formula <- formula(paste(y_var, " ~ " ,paste(x_var,collapse = " + ")))

cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

NF_seq <- IAS_seq <- MEC_model_seq <- train_model_list <- c()

for (model_i in model_list) {
  cat("This is model ",model_i,"\n")
  set.seed(13)
  train_model <- train(model_formula, data = ds_train, method = model_i,
                       trControl=trctrl,
                       # preProcess = c("center", "scale"),
                       tuneLength = 10)
  
  train_model_list <- c(train_model_list,list(train_model))
  
  y_hat_i <- predict(train_model,ds_test[,x_var])
  y_true <- ds_test[,y_var]
  plot(y_true~y_hat_i)
  summary(lm(y_true~y_hat_i))$adj
  
  
  
  NF <- sapply(x_var, function(x) feature_used(train_model,x,500) )
  NF_seq <- c(NF_seq,sum(NF))
  
  # compute ALE
  yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata,
                                                        type = "raw"))
  
  {
    f_0 <- mean(y_hat_i)
    f_ale_1st  <- f_0
    
    MEC <- V <- c()
    
    for (v in x_var) {
      j <- which(v == x_var)
      ale_j <- ALEPlot(ds_train, train_model, pred.fun = yhat, J = j,
                       K = 10^2, NA.plot = TRUE)
      
      
      loess_j <- loess(ale_j$f.values ~ ale_j$x.values,span = 0.1)
      
      {
        
        file.i <- paste0(fig.dir,model_i,"/ALE","_",j,".pdf")
        pdf(file.i)
        plot(ale_j$f.values~ale_j$x.values,
             ylab = "ALE", xlab = expression(x[j]),
             pch = 20, cex = 0.5)
        lines(predict(loess_j,ale_j$x.values)~ale_j$x.values,col = 2)
        grid(10)
        dev.off()
        
      }
      x_j <- ds_train[,v]
      f_j_ale <- predict(loess_j,x_j) # based on the ale, we compute the approximation
      
      # in this case, f_j_ale should correspond to
      lm_fit <- beta_vec[j]*(ds_train[,v] - mean(ds_train[,v]))
      plot(lm_fit ~ f_j_ale, pch = 20,cex = 0.5)
      abline(a = 0,b = 1, col = 2) # works for lm model
      
      # f_j_ale <- f_j_ale - mean(f_j_ale)
      f_ale_1st <- f_ale_1st + f_j_ale
      
      # MEC using segmented regression
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
        
        file.i <- paste0(fig.dir,model_i,"/MEC","_",j,".pdf")
        pdf(file.i)
        plot(seg_fit_j~x,pch = 20,type = "l", data = ds_plot,
             ylab = "ALE", xlab = expression(x[j]), lwd = 2)
        lines(y ~ x, data = ds_plot,col = 2)
        grid(10)
        dev.off()
        
        while(Rsq_j < 1 - epsilon) {
          cat("This is ",K,"\n")
          seg_j <- segmented(lm_j,seg.Z = ~ x, npsi = K)
          seg_fit_j <- fitted(seg_j)
          Rsq_j <- summary(lm(y~seg_fit_j))$adj
          K <- K + 1
          
          ds_plot <- data.frame(seg_fit_j,y,x)
          ds_plot <- ds_plot[order(ds_plot$x),]
          
          pdf(file.i)
          plot(seg_fit_j~x,pch = 20,type = "l", data = ds_plot,
               ylab = "ALE", xlab = expression(x[j]), lwd = 2)
          lines(y ~ x, data = ds_plot,col = 2)
          grid(10)
          dev.off()
        }
        
        return(K)
        
      }
      
      
      catch_error <- try(MEC_try_function(ale_j),silent  = T)
      try_i <- 1
      
      while(inherits(catch_error,"try-error")) {
        try_i <- try_i + 1
        cat("This is error trial ",try_i,"\n" )
        catch_error <- try(MEC_try_function(ale_j),silent  = T)
      }
      
      
      MEC_j <- catch_error
      V_j <- var(f_j_ale)
      
      MEC <- c(MEC,MEC_j)
      V <- c(V,V_j)
    }
    
    MEC_model <- sum(MEC*V/sum(V))
    
  }
  
  # compute IAS
  IAS1 <- sum((y_hat_i - f_ale_1st)^2)
  IAS2 <- sum((y_hat_i - f_0)^2)
  IAS <- IAS1/IAS2
  
  # alternatively we can do so using regressing
  1-summary(lm(y_hat_i~f_ale_1st))$adj
  IAS_seq <- c(IAS_seq,IAS)
  
  MEC_model_seq <- c(MEC_model_seq,MEC_model)
  
}


stopCluster(cl)
registerDoSEQ()

sum_df <- data.frame(model = model_list,  MEC = MEC_model_seq, IAS = round(IAS_seq,2),NF = NF_seq)
xtable(sum_df)


