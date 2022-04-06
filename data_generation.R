####Data Generation###

library(Matrix)
library(MASS)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

gen_corr_data <- function(corr_mat, f_sd, f_mu, num_obs){
   d = diag(f_sd)
   cov_mat = d %*% corr_mat %*% d
   return(mvrnorm(n = num_obs, mu = f_mu, Sigma = cov_mat))
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
  
  Y = 3*raw[,1] + 2*raw[,2]^3 + 1*raw[,3]
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
  data = read.csv('data/ff_factors.csv', header = TRUE)
  colnames(data) = c(1:length(data))
  
  cov_mat = cov(data)
  means = unlist(lapply(data, mean))
  #sds = lapply(data, sd)
  num = 1000
  
  raw = mvrnorm(n = num, mu = means, Sigma = cov_mat)
  # we need to add noise
  eps = rnorm(n = num, sd = 0.1)
  Y = 0.25*raw[,1] + 0.25*raw[,2]^2 + 0.25*raw[,3] + 0.25*raw[,4] + eps
  
  raw_df = data.frame(raw, Y)
  
  return(raw_df)
}


gen_calibrated_data()

