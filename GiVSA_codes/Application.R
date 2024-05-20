################################################################################
##### THIS FILE CONTAINS THE CODES FOR APPLYING M-GIVSA OR G-GIVSA METHODS #####
#####                       ON A GIVEN DATA SET                            ##### 
################################################################################

## 1. PREREQUISITES: If the following R-packages are not installed previously, 
##                   then it will be installed and loaded.  
if(!require(corpcor)){
  install.packages("corpcor")
  library(corpcor)
}
if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(lme4)){
  install.packages("lme4")
  library(corpcor)
}
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
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
source("GiVSA.R")
res_MGiVSA <<- Givsa(y = y, X = X_std, groups = true_group_structure, c = 0.025, 
              c0 = 0.01, grp_post = FALSE)
source("Output_GiVSA.r")
Extract.Results.MGiVSA <- Output_GiVSA(res_MGiVSA, X, y)
keep <- c("X", "X_std", "y", "true_group_structure", "Givsa", 
          "res_MGiVSA", "Extract.Results.MGiVSA", "Output_GiVSA", "keep")
rm(list = setdiff(ls(), keep))
res_GGiVSA <<- Givsa(y = y, X = X_std, groups = true_group_structure, c = 0.01, 
                     c0 = 0.0025, grp_post = TRUE)
Extract.Results.GGiVSA <- Output_GiVSA(res_GGiVSA, X, y)
keep <- c(keep, c("res_GGiVSA", "Extract.Results.GGiVSA"))
rm(list = setdiff(ls(), keep))

# Printing results
{
  print("########### Results of M-GiVSA ###########")
  print("The most visited model is:")
  print(Extract.Results.MGiVSA$HPM)
  
  print("Regression parameter of the most visited model is:")
  print(Extract.Results.MGiVSA$HPM_beta[Extract.Results.MGiVSA$HPM])
  
  print("The median probability model is:")
  print(Extract.Results.MGiVSA$MPM)
  
  print("Regression parameter of the median probability model is:")
  print(Extract.Results.MGiVSA$MPM_beta[Extract.Results.MGiVSA$MPM])

  print("The model with highest posterior probability is:")
  print(Extract.Results.MGiVSA$mode)
  
  print("Regression parameter of the highest posterior probability model is:")
  print(Extract.Results.MGiVSA$mode_beta[Extract.Results.MGiVSA$mode])
  print("##########################################")
  
  
  print("########### Results of G-GiVSA ###########")
  print("The most visited model is:")
  print(Extract.Results.GGiVSA$HPM)
  
  print("Regression parameter of the most visited model is:")
  print(Extract.Results.GGiVSA$HPM_beta[Extract.Results.GGiVSA$HPM])
  
  print("The median probability model is:")
  print(Extract.Results.GGiVSA$MPM)
  
  print("Regression parameter of the median probability model is:")
  print(Extract.Results.GGiVSA$MPM_beta[Extract.Results.GGiVSA$MPM])
  
  print("The model with highest posterior probability is:")
  print(Extract.Results.GGiVSA$mode)
  
  print("Regression parameter of the highest posterior probability model is:")
  print(Extract.Results.GGiVSA$mode_beta[Extract.Results.GGiVSA$mode])
  print("##########################################")
  
}


