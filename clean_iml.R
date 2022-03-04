
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


##WARNING: THIS WILL DELETE ALL FILES IN MODEL FOLDERS
#configs
{
#clear
rm(list = ls())

#graph output dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
output.dir <- getwd()
fig.dir <- paste(output.dir,"/Figures/",sep = "")

#model fitting params
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, allowParallel = T)
model_list <- c("lm","svmRadial")

#parallel param
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

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

###Data Generation###
{
gen_corr_data <- function(corr_mat, f_sd, f_mu, num_obs){
   d = diag(f_sd)
   cov_mat = d %*% corr_mat %*% d
   return(mvrnorm(n = num_obs, mu = f_mu, Sigma = cov_mat))
}

corr_mat = matrix( c(1.0, 0.8, 0.2,
                     0.8, 1.0, 0.3,
                     0.2, 0.3, 1.0),
                   nrow = 3)

sd = c(2, 1, 0.5)
m = c(0, 0, 0)
num = 1000

raw = gen_corr_data(corr_mat, sd, m, num)

Y = 3*raw[,1] + 2*raw[,2]^3 + 1*raw[,3]
raw_df = data.frame(raw, Y)

#plot
pairs(raw_df)
}


#initialize vars
{
ds <- raw_df
ds_train <- ds
ds_test <- ds_train
y_var <- "Y"
x_var <- names(ds)[!names(ds) %in% y_var]
model_formula <- formula(paste(y_var, " ~ " ,paste(x_var,collapse = " + ")))

#other vars
MSE_MEC_prop_seq <- MSE_seq <- NF_seq <- IAS_seq <- MEC_model_seq <- train_model_list <- c()
}

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
   
   MSE = mean((ds_test[,"Y"] - y_hat_i)^2)
   MSE_seq <- c(MSE_seq, MSE)
   
   NF <- sapply(x_var, function(x) feature_used(train_model,x,500) )
   NF_seq <- c(NF_seq, sum(NF))
   
   # compute ALE
   yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata, type = "raw"))
   
   {
      f_0 <- mean(y_hat_i)
      f_ale_1st  <- f_0
      
      MEC <- V <- c()
      
      for (v in x_var) {
         
         ## Do ALE
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
         
         # # in this case, f_j_ale should correspond to
         # lm_fit <- beta_vec[j]*(ds_train[,v] - mean(ds_train[,v]))
         # plot(lm_fit ~ f_j_ale, pch = 20,cex = 0.5)
         # abline(a = 0,b = 1, col = 2) # works for lm model
         
         # f_j_ale <- f_j_ale - mean(f_j_ale)
         f_ale_1st <- f_ale_1st + f_j_ale
         
         ## MEC using segmented regression
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
      
      #proportion measure
      MSE_MEC_prop_seq <- c(MSE_MEC_prop_seq, MSE * MEC_j)
      MEC_model <- sum(MEC*V/sum(V))
      
   }
   
   # compute IAS
   IAS1 <- sum((y_hat_i - f_ale_1st)^2)
   IAS2 <- sum((y_hat_i - f_0)^2)
   IAS <- IAS1/IAS2
   
   # alternatively we can do so using regressing
   #1-summary(lm(y_hat_i~f_ale_1st))$adj
   IAS_seq <- c(IAS_seq,IAS)
   
   MEC_model_seq <- c(MEC_model_seq, MEC_model)
   
}

stopCluster(cl)
registerDoSEQ()

sum_df <- data.frame(model = model_list,  
                     MEC = MEC_model_seq, 
                     IAS = round(IAS_seq,2), 
                     NF = NF_seq, 
                     MSE = MSE_seq, 
                     MSE_MEC_prop = MSE_MEC_prop_seq)
xtable(sum_df)

sum_df
#introducing a metric for error to (MEC )
#where we have an error function where the best performance is 0, and values at >= 0
#then we will reduce the MSE by the complexity of the model
#the goal is that the model with lower MSE and lower model complexity will have a higher score

##DONE##
#Calibrated data

##TODO##
#slope retreival function
#put majeed code into nice function
   #figure out the try catch thing



###working with the ale stuff
#need the ale plots to get the solpes at each point

