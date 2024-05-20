################################################################################
#####   THIS FILE CONTAINS THE CODES FOR APPLYING GSISd METHOD ON A GIVEN  #####
#####              ON A GIVEN DATA SET AND GROUP STRUCTUE                  ##### 
################################################################################

## 1. PREREQUISITES: If the following R-packages are not installed previously, 
##                   then it will be installed and loaded.  
if(!require(corpcor)){
  install.packages("corpcor")
  library(corpcor)
}
## 2. DATA PROCESSING: The Data should be available as an Rdata file named  
##    "input.RData" containing the following objects:
## (a) X: The design matrix in nXp format, where n is the no. of observations,  
##        and p is the number of covariates.
## (b) y: The n vector of responses
## (c) true_group_structure: A p-length vector indicating the groups indices as  
##                           numbers. For example, for p = 6, if X1, X3, X5 
##                           belong to group 1, X2, X6 belong to group 2, and X4 
##                           belong to group 3, then the true_group_structure 
##                           would be c(1, 2, 1, 3, 1, 2)
rm(list=ls())
set.seed(1947)
load("input.RData")
p <- ncol(X)
X_std <- X
for (j in 1:p) {
  X_std[, j] <- (X[, j] - mean(X[, j])) / sd(X[, j])
}
source("GSIS.R")
d <- 1 # Choose your d 
res_GSIS <- GSISd(y = y, X = X_std, true_group_structure = true_group_structure, 
                 d = d)
group.order <- order(res_GSIS, decreasing = TRUE)
prob.order <- sort(res_GSIS, decreasing = TRUE)
# Printing results
{
  print("The five highest GSISd probability groups are:")
  print(group.order[1:5])
  print("The GSISd probabilities of the above five groups are:")
  print(prob.order[1:5])
}
################################################################################
