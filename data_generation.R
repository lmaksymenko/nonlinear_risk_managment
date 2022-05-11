####Data Generation###

library(Matrix)
library(MASS)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

gen_corr_data <- function(corr_mat, f_sd, f_mu, num_obs){
   d = diag(f_sd)
   cov_mat = d %*% corr_mat %*% d
   return(mvrnorm(n = num_obs, mu = f_mu, Sigma = cov_mat))
}

gen_2_ord <- function(data, weights){
  
  Y = 0
  i = 1
  #first order
  for(j in len(data)){
    
  }
  #second order
  
}


gen_manual_data <- function(){
  corr_mat = matrix( c(1.0, 0.8, 0.2,
                       0.8, 1.0, 0.3,
                       0.2, 0.3, 1.0),
                     nrow = 3)
  
  sd = c(2, 1, 0.5)
  m = c(0, 0, 0)
  num = 1000
  
  raw = gen_corr_data(corr_mat, sd, m, num)
  
  eps = rnorm(n = num, sd = 0.001)
  
  Y = 1*raw[,1] + 1*raw[,2] + 1*raw[,3] 
      + 1*raw[,1]*raw[,2] + 1*raw[,2]*raw[,3] + 1*raw[,1]*raw[,3]
      + 1*raw[,1]^2 + 1*raw[,2]^2 + 1*raw[,3]^2 + eps
  
  
  raw_df = data.frame(raw, Y)
  
  return(raw_df)
}


### GENERATING CALIBRATED DATA ###

# skip risk free rate
# fix data sampling
#think about nonlinear structure - add more nonlinearity (squared, higher orders)
# what are the correlations between higher orders (x1^2 and so on)
# plot data to view overall complexity (gg plot)


#using https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

gen_calibrated_data <- function(){
  data = read.csv('data/monthlyffdata.csv', header = TRUE)
  colnames(data)
  
  #remove risk free rate
  if(c("RF") %in% colnames(data)){
    data = data[, -which(names(data) == "RF")]
  }
  
  colnames(data) = c(1:length(data))
  
  print(cov(data))
  
  #divide by 20 because the data is monthly 
  #divide by 100 bc data as a percent (eg: 2.58%)
  cov_mat = cov(data)/20/100
  means = unlist(lapply(data, mean))/20/100
  #sds = lapply(data, sd)
  num = 100
  
  
  raw = mvrnorm(n = num, mu = means, Sigma = cov_mat)
  
  
  eps = rnorm(n = num, sd = 0.001)#noise
  
  Y = 1*raw[,1] + 1*raw[,2] + 1*raw[,3] 
  + 1*raw[,1]*raw[,2] + 1*raw[,2]*raw[,3] + 1*raw[,1]*raw[,3]
  + 1*raw[,1]^2 + 1*raw[,2]^2 + 1*raw[,3]^2 + eps
  
  raw_df = data.frame(raw, Y)
  
  return(raw_df)
}



pairs(gen_calibrated_data())



###Recession regime data
combn(4,2)

