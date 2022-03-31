#Imports
{
   #clear env
   rm(list = ls())
   
   #graph output dir
   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
   output.dir <- getwd()
   fig.dir <- paste(output.dir,"/Figures/",sep = "")
   
   #packages
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
   
   #imports
   source("data_generation.R")
}




clear_figures <- function(model_list)
{
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


fit_model <- function(formula, model, x_var, y_var, train, test, control)
{
   cat("Fitting ", model, "\n")
   
   train_model <- train(formula, data = train, method = model,
                        trControl = control,
                        # preProcess = c("center", "scale"),
                        tuneLength = 10)
   print('Model:')
   print(class(train_model))
   
   y_hat_i <- predict(train_model, test[,x_var])
   y_true <- test[,y_var]
   
   #plot(y_true~y_hat_i)
   #summary(lm(y_true~y_hat_i))$adj
   
   mse = mean((test[,"Y"] - y_hat_i)^2)
   
   NF <- sapply(x_var, function(x) feature_used(train_model, x, 500) )
   num_feat_used = sum(NF)
   
   return(c(train_model, mse, y_hat_i, num_feat_used))
}


#using the predicted y_hats and the model calculate the ale plot for each x_variable
#return the ale values (ig thats all we want)
#

fit_ale <- function(model, x_var, train, y_hat_i)
{
   ale_Js <- c()
   var_Js <- c()
   loess_Js <- c()
   
   f_0 <- mean(y_hat_i)
   f_ale_1st  <- f_0
   
   yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata, type = "raw"))#error here
   
   for (v in x_var) {
      
      print(class(train))
      print(class(v))
      
      j <- which(v == x_var)
      
      ale_j <- ALEPlot(train, model, pred.fun = yhat, J = j, #this is where the error is 
                       K = 10^2, NA.plot = TRUE)#K will need to be changed later
      
      ale_Js <- c(ale_Js, list(ale_j))
      
      loess_j <- loess(ale_j$f.values ~ ale_j$x.values, span = 0.1)
      loess_Js <- c(loess_Js, list(loess_j))
      
      f_j_ale <- predict(loess_j, train[,v]) # based on the ale, we compute the approximation, exists for each factor
      
      # f_j_ale <- f_j_ale - mean(f_j_ale)
      f_ale_1st <- f_ale_1st + f_j_ale
      
      var_Js <- c(var_Js, var(f_j_ale))
      
      {
         file.i <- paste0(fig.dir, model,"/ALE","_",j,".pdf")
         pdf(file.i)
         plot(ale_j$f.values~ale_j$x.values,
              ylab = "ALE", xlab = expression(x[j]),
              pch = 20, cex = 0.5)
         lines(predict(loess_j, ale_j$x.values)~ale_j$x.values, col = 2)
         grid(10)
         dev.off()
      }
   }
   
   # compute IAS
   IAS1 <- sum((y_hat_i - f_ale_1st)^2)
   IAS2 <- sum((y_hat_i - f_0)^2)
   IAS <- IAS1/IAS2
 
   return(list(ale_Js, var_Js, IAS, loess_Js)) #return the ale_j plots, V_j
   #returning in this format flattens the lists
}


calc_mec <- function(model, x_var, ale_Js, tol, V)#ale_j (the plot), tol is epsilon, 
{
   epsilon <- tol # set tolerance
   MEC <- c()
   
   for (v in x_var){
      ale_j = ale_Js[[which(v == x_var)]]
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
            j <- which(v == x_var)
            file.i <- paste0(fig.dir, model,"/MEC","_",j,".pdf")
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
      
      ####
      catch_error <- try(MEC_try_function(ale_j),silent  = T)
      try_i <- 1

      while(inherits(catch_error,"try-error")) {#why is this here
         try_i <- try_i + 1
         cat("This is error trial ", try_i, "\n" )

         catch_error <- try(MEC_try_function(ale_j), silent  = T)

         if(try_i > 1000){
            cat("1000 trials reached, stopping")
            break
         }
      }
      
      ####
      
      MEC_j <- catch_error
      MEC <- c(MEC,MEC_j)
      
   }

   MEC_model <- sum(MEC*V/sum(V))
   
   return(MEC_model)
}

get_xj_slope <- function(loess_model, vec_xval){
   
   sapply(vec_xval, function(xval) return ((predict(loess_model, xval) 
                                                   - predict(loess_model, xval - 0.0001))
                                                  / (xval - (xval-0.0001))) )
}


backtest_hedging <- function(model, loess_model_list, data, x_var, cov_mat){
   ##backtests the hedging for a single model
   
   #TESTING VARS
   loess_model_list = tmp_ale[[4]]
   data = ds
   #
   
   
   ###get the weights for each asset
      #"Y" "X1" ... in this order
         #vars = c("Y", x_var)
   
      #get the returns matrix
         #ds[c("Y", x_var)] #indexing the returns
   
   #get beta matrix for X vars
   betas = sapply(x_var, function(x) return (get_xj_slope(loess_model_list[[which(x == x_var)]], data[x])) )
   
   ###apply function to each row in the testing data
  
   
   ###have a vector of daily returns
      #calculate the evaluation metrics
} 



main <- function()
{
   cl <- makePSOCKcluster(detectCores())
   registerDoParallel(cl)
   
   #model fitting params
   model_list <- c("lm","svmRadial")
   trctrl <- trainControl(method = "repeatedcv", 
                          number = 10, 
                          repeats = 3, 
                          allowParallel = T) # i think this does cross validation
   clear_figures(model_list)
   
   ##################
   #make data
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
      
      Y = 3*raw[,1] + 2*raw[,2] + 1*raw[,3]
      raw_df = data.frame(raw, Y)
      
   }
   
   
   #train test split
   ds <- gen_manual_data()

   ################
   
   #REDO
   ds_train <- ds
   ds_test <- ds_train
   
   #model formula
   y_var <- "Y"
   x_var <- names(ds)[!names(ds) %in% y_var]
   model_formula <- formula(paste(y_var, " ~ " ,paste(x_var,collapse = " + ")))
   
   #info vars
   MSE_MEC_prop_seq <- MSE_seq <- NF_seq <- IAS_seq <- MEC_model_seq <- train_model_list <- c()
   
   #beta extraction var
   loess_model_seq <- c()
   
   
   for (model_i in model_list) {
      
      cat("Fitting ", model_i, "\n")
      
      #### TRAIN MODEL ####
      train_model <- train(model_formula, 
                           data = ds_train, 
                           method = model_i,
                           trControl = trctrl,
                           # preProcess = c("center", "scale"),
                           tuneLength = 10)
      
      y_hat_i <- predict(train_model, ds_test[,x_var])
      y_true <- ds_test[,y_var]
      
      mse = mean((ds_test[,"Y"] - y_hat_i)^2)
      MSE_seq <- c(MSE_seq, mse)
      ####
      
      
      NF <- sapply(x_var, function(x) feature_used(train_model, x, 500) )
      NF_seq <- c(NF_seq, sum(NF))

      
      #### ALE
      tmp_ale = fit_ale(train_model, x_var, ds_train, y_hat_i)
      IAS_seq <- c(IAS_seq, tmp_ale[[3]])
      loess_model_seq <- c(loess_model_seq, list(tmp_ale[[4]]))
      
             
      #### MEC
      tmp_mec = calc_mec(train_model, x_var, tmp_ale[[1]], 0.05, tmp_ale[[2]])
      MEC_model_seq <- c(MEC_model_seq, tmp_mec)
      
   }
   stopCluster(cl)
   registerDoSEQ()
   
   #### SUMMARY MODELS ####
   sum_df <- data.frame(model = model_list,
                        MEC = MEC_model_seq,
                        IAS = round(IAS_seq,2),
                        NF = NF_seq,
                        MSE = MSE_seq,
                        MSE_MEC_prop = MSE_seq * MEC_model_seq)
   print(sum_df)
   ####
   
   #### HEDGING ####
   
   ####
   
   
   
}

main()

