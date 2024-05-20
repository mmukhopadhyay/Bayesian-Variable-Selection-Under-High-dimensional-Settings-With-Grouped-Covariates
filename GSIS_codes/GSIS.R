################################################################################
##      THIS FILE CONTAINS A FUNCTION TO CALCULATE THE GSIS GROUP INCLUSION   ##
##        PROBABILITIES FOR A GIVEN DATA, GROUP STRUCTURE AND CHOICE OF d     ##
################################################################################
## 1. GSISd function
#' @description
#' This function calculates the GSISd group inclusion probabilities for a given 
#' data, provided the underlying group stucture and choice of d
#' @param y The n-vector of n observations on the response variable.
#' @param X The nXp design matrix of the n observations on the p covariates.
#' @param true_groups_structure  A vector indicating the groups indices as  
#'                               numbers. For example, for p = 6, if X1, X3, X5 
#'                               belong to group 1, X2, X6 belong to group 2, 
#'                               and X4 belong to group 3, then the groups 
#'                               vector would be c(1, 2, 1, 3, 1, 2)
#' @param d The choice of d, 0 <= d <= 1.
#' @return A vector with length as the number of groups containing the GSIS 
#'         group inclusion probabilities. Note that, if the group indices are 1
#'         to k, then the vector will have probabilities of group 1 to group k
#'         in order.
#' @usage GSISd(y = y, X = X_std, true_group_structure = true_group_structure, 
#'        d = d)

GSISd <- function(y, X, truth, true_group_structure, d) {
  library("corpcor")
  try(
    {
      Groups <- 1:max(true_group_structure)
      YTX <- t(y) %*% X
      group_elements <- lapply(Groups, function(x) {
        which(true_group_structure == x)
      })
      group_mu <- numeric(length = length(Groups))
      for (i in Groups) {
        group.i <- group_elements[[i]] 
        X.group <- X[, group.i]
        y_g <- X.group %*% (pseudoinverse(t(X.group) %*% X.group + 
                        1e-3 * diag(1, length(group.i))) %*% (YTX[1, group.i]))
        group_mu[i] <- as.numeric(t(y_g) %*% y_g)/
                                        (length(as.numeric(group.i)))^d
      }
      res <- group_mu / sum(group_mu)
      
    },
    silent = FALSE
  )
  return(res)
}
