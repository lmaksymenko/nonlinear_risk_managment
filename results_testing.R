
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)

path = "backtest_results"

cube = read.table(paste(path, "/", "lm_svm_allcube_trials.csv", sep = ""), header = TRUE)
square = read.table(paste(path, "/", "lm_svm_allsquare_trials.csv", sep = ""), header = TRUE)
linear = read.table(paste(path, "/", "lm_svm_linear_trials.csv", sep = ""), header = TRUE, sep = ",")

##Goal: Get Return Closest to Zero?

####
data = data.frame(type = c( rep("Linear", 500), rep("SVM Radial", 500) ),
                  return = c(cube$Hedge_Ret1, cube$Hedge_Ret2) )

ggplot(data = data, aes(x = return, fill = type)) + 
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + 
  labs(title = "All Factors Cubed")


g = cube
cat("CUBE \n")
cat("Linear Return SD: ", sd(g$Hedge_Ret1), "\n")
cat("Linear Mean MEC: ", mean(g$MEC2), "\n")
cat("----------------------------", "\n")
cat("SVM Return SD: ", sd(g$Hedge_Ret2), "\n")
cat("SVM Mean MEC: ", mean(g$MEC2), "\n")
cat("\n\n")


####
data = data.frame(type = c( rep("Linear", 500), rep("SVM Radial", 500) ),
                  return = c(square$Hedge_Ret1, square$Hedge_Ret2) )

ggplot(data = data, aes(x = return, fill = type)) + 
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + 
  labs(title = "All Factors Squared")

g = cube
cat("SQUARED \n")
cat("Linear Return SD: ", sd(g$Hedge_Ret1), "\n")
cat("Linear Mean MEC: ", mean(g$MEC2), "\n")
cat("----------------------------", "\n")
cat("SVM Return SD: ", sd(g$Hedge_Ret2), "\n")
cat("SVM Mean MEC: ", mean(g$MEC2), "\n")
cat("\n\n")

####
data = data.frame(type = c( rep("Linear", 500), rep("SVM Radial", 500) ),
                  return = c(linear$Hedge_Ret1, linear$Hedge_Ret2) )

ggplot(data = data, aes(x = return, fill = type)) + 
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + 
  labs(title = "All Factors Linear")

g = cube
cat("LINEAR \n")
cat("Linear Return SD: ", sd(g$Hedge_Ret1), "\n")
cat("Linear Mean MEC: ", mean(g$MEC2), "\n")
cat("----------------------------", "\n")
cat("SVM Return SD: ", sd(g$Hedge_Ret2), "\n")
cat("SVM Mean MEC: ", mean(g$MEC2), "\n")
cat("\n\n")
