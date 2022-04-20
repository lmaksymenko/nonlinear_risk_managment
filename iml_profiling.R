#Imports

library(profvis)

profvis({
{
   #clear env
   rm(list = ls())
   
   #graph output dir
   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
   output.dir <- getwd()
   fig.dir <- paste(output.dir,"/Figures_Profiling/",sep = "")
   
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
   library(logging)

   
   options(warn = defaultW)
   
   #imports
   source("data_generation.R")
   #source("ale.R")
   #source("mec.R")
   #source("backtesting.R")
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
   
   #plot(prediction1~prediction2)
   
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

fit_ale <- function(model, x_var, train, y_hat_i, plotting)
{
   ale_Js <- c()
   var_Js <- c()
   loess_Js <- c()
   
   f_0 <- mean(y_hat_i)
   f_ale_1st  <- f_0
   
   yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata, type = "raw"))#error here
   
   for (v in x_var) {
      
      j <- which(v == x_var)
      
      ale_j <- ALEPlot(train, model, pred.fun = yhat, J = j,
                       K = 10^2, NA.plot = TRUE)#K will need to be changed later?
      
      ale_Js <- c(ale_Js, list(ale_j))
      
      loess_j <- loess(ale_j$f.values ~ ale_j$x.values, span = 0.1)
      loess_Js <- c(loess_Js, list(loess_j))
      
      f_j_ale <- predict(loess_j, train[,v]) # based on the ale, we compute the approximation, exists for each factor
      
      # f_j_ale <- f_j_ale - mean(f_j_ale)
      f_ale_1st <- f_ale_1st + f_j_ale
      
      var_Js <- c(var_Js, var(f_j_ale))
      
      if(plotting){
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
   sum((as.numeric(y_hat_i) - f_ale_1st)^2)
   IAS2 <- sum((y_hat_i - f_0)^2)
   IAS <- IAS1/IAS2
 
   return(list(ale_Js, var_Js, IAS, loess_Js)) #return the ale_j plots, V_j
   #returning in this format flattens the lists
}


calc_mec <- function(model, x_var, ale_Js, tol, var_Js, plotting)#ale_j (the plot), tol is epsilon, 
{
   #model - current trained model
   #x_var - list of x variables in model
   #ale_Js - list of ALEPlot objects, one for each x_var
   #tol - peicewise regression error tolerance in %
   #var_Js - variance of the loess predicted f_j_ale
   
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
         
         if(plotting){
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
            seg_j <- segmented(lm_j, seg.Z = ~ x, npsi = K)#this could potentially be rewritten, there is a control for segmented that allows us to specify tol
            seg_fit_j <- fitted(seg_j)
            Rsq_j <- summary(lm(y~seg_fit_j))$adj
            K <- K + 1
            
            ds_plot <- data.frame(seg_fit_j,y,x)
            ds_plot <- ds_plot[order(ds_plot$x),]
            
            if(plotting){
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
      catch_error <- try(MEC_try_function(ale_j),silent  = F)
      try_i <- 1

      while(inherits(catch_error,"try-error")) {#why is this here
         try_i <- try_i + 1
         cat("This is error trial ", try_i, "\n" )

         catch_error <- try(MEC_try_function(ale_j), silent  = F)

         if(try_i > 5){
            cat("5 trials reached, stopping")
            catch_error <- length(x_var)
            break
         }
      }
      
      ####
      
      MEC_j <- catch_error
      MEC <- c(MEC,MEC_j)
      
   }

   MEC_model <- sum(MEC*var_Js/sum(var_Js))
   
   return(MEC_model)
}

get_xj_slope <- function(loess_model, vec_xval){
   
   sapply(vec_xval, function(xval) return ((predict(loess_model, xval) 
                                                   - predict(loess_model, xval - 0.0001))
                                                  / (xval - (xval-0.0001))) )
}


iid_backtest_returns <- function(model, loess_model_list, data, x_var){
   ##iid backtest the hedging for a single model
   
   #get beta matrix for X vars
   betas = sapply(x_var, function(x) return (get_xj_slope(loess_model_list[[which(x == x_var)]], data[x])) )
   
   #don't normalize betas yet
   #betas = t(apply(betas, 1, function(x) return(x/sum(x)))) #normalize so that they add to 1
   
   betas = cbind(-betas, Y = rep(1,nrow(betas)))
   
   returns = as.vector(t(apply(betas * data, 1, sum)))
   
   g = data.frame(data[["Y"]],returns)
   
   #standard deviation of the returns
   base = (sd(g$data...Y...))
   hedged = (sd(g$returns, na.rm = TRUE))

   return(c(base,hedged))
} 

iml <- function()
{
   cl <- makePSOCKcluster(detectCores())
   registerDoParallel(cl)
   
   
   #### CONFIGS ####
   model_list <- c("lm", "svmRadial")
   TR_CTRL <- trainControl(method = "repeatedcv", #model fitting params
                          number = 10, 
                          repeats = 3, 
                          allowParallel = T) # i think this does cross validation
   
   #data and train test split for fixed window
   TRAIN_SIZE = 0.5
   PLOTTING = FALSE
   ##################
   clear_figures(model_list)
   
   #ds <- gen_calibrated_data()
   ds <- read.csv("last_data.csv", headers = TRUE)
   train <- ds[0: (nrow(ds) * TRAIN_SIZE),]
   test <- ds[(nrow(ds) * TRAIN_SIZE):nrow(ds),]
   
   #write.csv(ds, "last_data_debug.csv")
   
   #model formula
   y_var <- "Y"
   x_var <- names(ds)[!names(ds) %in% y_var]
   model_formula <- formula(paste(y_var, " ~ " ,paste(x_var,collapse = " + ")))
   
   #info vars
   NF_seq <- IAS_seq <- MEC_model_seq <- train_model_list <- c()
   
   MSE_MEC_prop_seq <- MSE_test_seq <-MSE_train_seq <- c()
   
   #beta extraction var
   loess_model_seq <- c()
   
   
   for (model_i in model_list) {
      
      cat("Fitting ", model_i, "\n")
      
      #### TRAIN MODEL ####
      train_model <- train(model_formula, 
                           data = train, 
                           method = model_i,
                           trControl = TR_CTRL,
                           # preProcess = c("center", "scale"),
                           tuneLength = 10)
      
      y_hat_train <- predict(train_model, train[, x_var])
      y_hat_test <- predict(train_model, test[, x_var])
      y_true <- test[,y_var]
      
      mse = mean((train[,"Y"] - y_hat_train)^2)
      MSE_train_seq <- c(MSE_train_seq, mse)
      
      mse = mean((test[,"Y"] - y_hat_test)^2)
      MSE_test_seq <- c(MSE_test_seq, mse)
      ####
      
      
      NF <- sapply(x_var, function(x) feature_used(train_model, x, 500) )
      NF_seq <- c(NF_seq, sum(NF))

      
      #### ALE
      tmp_ale = fit_ale(train_model, x_var, train, y_hat_train, PLOTTING)
      IAS_seq <- c(IAS_seq, tmp_ale[[3]])
      loess_model_seq <- c(loess_model_seq, list(tmp_ale[[4]]))
      
     
      #### MEC
      tmp_mec = calc_mec(train_model, x_var, tmp_ale[[1]], 0.05, tmp_ale[[2]], PLOTTING) ## change to 90% complexity
      MEC_model_seq <- c(MEC_model_seq, tmp_mec)
      
   }
   stopCluster(cl)
   registerDoSEQ()
   
   base_ret_seq <- c()
   hedge_ret_seq <- c()
   
   #### HEDGING TEST ####
   for (model_i in model_list) {
      cat(model_i, "backtest \n")
      tmp = iid_backtest_returns(model_i, loess_model_seq[[which(model_i == model_list)]], test, x_var)
      base_ret_seq = c(base_ret_seq, tmp[1])
      hedge_ret_seq = c(hedge_ret_seq, tmp[2])
   }
   ####
   
   #### SUMMARY MODELS ####
   sum_df <- data.frame(model = model_list,
                        MEC = MEC_model_seq,
                        IAS = round(IAS_seq,2),
                        NF = NF_seq,
                        MSE_Train = MSE_train_seq,
                        MSE_Test = MSE_test_seq,
                        MSE_Test_MEC_prop = MSE_test_seq * MEC_model_seq,
                        Base_Ret = base_ret_seq,
                        Hedge_Ret = hedge_ret_seq)
   
   return(sum_df)
}

dist_test <- function(){
   NUM_SIMS = 500
   
   #first wirte with colnames
   cat("Trial: ", 1, "\n")
   res = t(unlist(iml()))
   write.table(res ,"test.csv", row.names = FALSE)
   
   #the rest
   for(i in 2:NUM_SIMS){
      cat("Trial: ", i, "\n")
      res = t(unlist(iml()))
      write.table(res ,"test.csv", append = TRUE, row.names = FALSE, col.names = FALSE)
      
   }
}

#dist_test()
iml()

})
