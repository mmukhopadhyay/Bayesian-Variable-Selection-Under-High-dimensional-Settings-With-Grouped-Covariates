################################################################################
# THIS FILE CONTAINS CODES FOR EVALUATING THE RESULTS FROM THE Givsa FUNCTION ##
################################################################################
#' @param res An output of the Givsa function. Typically, a list containing two 
#'            elements, "samples" and "mode". The first component "samples", is 
#'            a list of post burnin visited models, while the second component, 
#'            "mode", is the indices of the most frequently visited model after 
#'            burnin.
#' @param X The nXp design matrix of the n observations on the p covariates.
#' @param y The n-vector of n observations on the response variable.

#' @return The output contains a list of six components: 
#'         (a) the indices of the post burnin most visited model (HPM), 
#'         (b) the regression coefficients corresponding to the HPM model (HPM_beta), 
#'         (c) the indices of the median probability model (MPM), 
#'         (d) the regression coefficients corresponding to the MPM model (MPM_beta), 
#'         (e) the model with highest posterior probability (mode), and 
#'         (f) the regression coefficients corresponding to mode model.

#' @usage res <- Givsa(y = y, X = X_std, groups = true_group_structure, 
#'                c = 0.025, c0 = 0.01, grp_post = FALSE)
#'        Extract.Results <- Output_GiVSA(res, X, y)
  

Output_GiVSA <- function(res, X, y){
  Models <- res$sample
  p <- ncol(X)
  chain <- matrix(0, nrow = length(Models), ncol = p)
  Model <- list()
  N <- nrow(chain)
  for (i in 1:length(Models)) {
    chain[i, Models[[i]]] <- 1
    Model[i] <- paste(Models[[i]], collapse = "_")
  }
  Model <- unlist(Model)
  Model <- sort(table(Model), decreasing = TRUE)
  visits <- (apply(chain, MARGIN = 2, sum))
  prob <- visits / nrow(chain)
  coefficients <- rep(0, p)

  HPM_param <- sort(as.numeric(unlist(strsplit(names(Model)[1], "_"))))
  if (length(HPM_param) == 0)  {
      HPM <- coefficients
  } else {
      X_mask <- X[, HPM_param]
      HPM2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
      HPM <- coefficients
      HPM[HPM_param] <- HPM2
  }

  MPM_param <- which(prob >= 0.5)
  if (length(MPM_param) == 0) {
      MPM <- coefficients
  } else {
      X_mask <- X[, MPM_param]
      MPM2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
      MPM <- coefficients
      MPM[MPM_param] <- MPM2
  }

  mode_param <- sort(res$mode)
  if (length(mode_param) == 0) {
      mode <- coefficients
  } else {
      X_mask <- X[, mode_param]
      mode2 <- solve(t(X_mask) %*% X_mask) %*% (t(X_mask) %*% y)
      mode <- coefficients
      mode[mode_param] <- mode2
  }
  ret <- list(HPM_param, HPM, MPM_param, MPM, mode_param, mode)
  names(ret) = c("HPM", "HPM_beta", "MPM", "MPM_beta", "mode", "mode_beta")
  return(ret)

  print("The highest probability model is:")
  print(HPM_param)

  print("The median probability model is:")
  print(MPM_param)

  print("The model with highest posterior probability is:")
  print(mode_param)
}