################################################################################
######### THIS FILE CONTAINS THE CODE FOR THE GIVSA ALGORITHM APPLIED TO #######
#########  THE MODIFIED G-PRIOR AS WELL AS THE GROUP-INFORMED G-PRIOR    #######
################################################################################

## INPUTS OF Givsa

#' @param y The n-vector of n observations on the response variable.
#' @param X The nXp design matrix of the n observations on the p covariates.
#' @param Iter  Total number of iterations.
#' @param burnin  The burnin period.
#'                Must be strictly less than the number of iterations Iter.
#' @param groups  A vector indicating the groups indices as numbers. 
#'                For example, for p = 6, if X1, X3, X5 belong to group 1,
#'                X2, X6 belong to group 2, and X4 belong to group 3, then the
#'                groups vector would be c(1, 2, 1, 3, 1, 2)
#' @param grp_post Logical, if TRUE then the posterior probabilities of the 
#'                 group-informed g-prior will be employed, if FALSE then the
#'                 modified-g prior will be employed.
#' @param min_group_size The minimum size of groups beyond which the groups will 
#'                       be merged together to form a single group for faster 
#'                       implementation   
#' @param c0 The hyperparameter g is set to g = c0 * p^(2 * (1 + delta))) as per
#'           the recommendations of obtained from the theoretical results. 
#'           This parameter is the c0 in the choice of g.
#' @param delta The hyperparameter g is set to g = c0 * p^(2 * (1 + delta))) as 
#'              per the recommendations of obtained from the theoretical results.
#'               This parameter is the delta in the choice of g. 
#' @param c  The hyperparameter kn is set to kn = c * n^r as per the 
#'           recommendations of obtained from the theoretical results. This
#'           parameter is the c in the choice of kn. 
#' @param r  The hyperparameter kn is set to kn = c * n^r as per the 
#'           recommendations of obtained from the theoretical results. This
#'           parameter is the r in the choice of kn.
#' @param al A convex combination of the GSIS group inclusion probability with 
#'           d = 0, and uniform probability on groups is chosen as the required 
#'           group probability in GiVSA. The concerned parameter al is the 
#'           weight of GSIS (d = 0) group inclusion probabilities. 
#' @return  A list containing two elements, "samples" and "mode". The first 
#'          component "samples", is a list of post burnin visited models, while
#'          the second component, "mode", is the indices of the most frequently 
#'          visited model after burnin.
#' @usage Givsa(y = y, X = X_std, groups = true_group_structure, c = 0.025, 
#'        c0 = 0.01, grp_post = FALSE)

Givsa <- function(y, X, Iter = 1e3, burnin = NULL, groups = NULL, 
                  grp_post = FALSE, min_group_size = 2, delta = 1e-2, r = NULL, 
                  c0 = 1e-2, c = 0.025, al = 2 / 3) {
  
  
  ######################
  ##### LIBRARIES ######
  ######################
  library(corpcor)
  library(Rcpp)
  library(lme4)
  library(parallel)
  library(Matrix)
  
  ######################
  ##### FUNCTIONS ######
  ######################
  power <- function(k, n) 
  {
    return(log(k) / log(n))
  }
  
  ######################
  ######## DATA ########
  ######################
  p <- ncol(X)
  n <- nrow(X)
  
  ################################################
  ### IF THE PARAMETER VALUES ARE NOT PROVIDED ###
  ################################################
  if (is.null(Iter)) 
  {   Iter <- 5000  }

  if (is.null(burnin)) 
  {   burnin <- Iter %/% 2  }
  
  if (is.null(groups)) 
  {  groups <- rep(1, p)  }
  
  ###########################
  ###### PREPROCESSING ######
  ###########################
  YTX <- t(y) %*% X
  norm_y <- as.numeric(t(y) %*% y)
  
  # Clubbing independent variables in a single cluster and 
  # assigning the last cluster identity to the merged group
  var2group <- groups
  # No of groups
  group_count <- max(var2group)
  group_sizes <- numeric(group_count)
  for (i in 1:group_count) 
  {
    group_sizes[i] <- length(which(var2group == i))
  }
  Flag <- 0
  if(is.null(group_sizes < min_group_size) == FALSE) {
    Flag <- 1
  }
  
  if(Flag == 1)
  {
    count <- group_count + 1
    count2 <- 1
    no.singleton <- 0
    for (i in 1:group_count) 
    {
      if (group_sizes[i] >= min_group_size) {
        var2group[var2group == i] <- count2
        count2 <- count2 + 1
      } else {
        var2group[var2group == i] <- count
        no.singleton <- no.singleton + 1
      }
    }
    var2group[var2group == count] = count - no.singleton
      
    # Size of each group
    group_sizes <- table(var2group)
    
    # No of groups
    group_count <- length(group_sizes)
  }
  
  ######################
  ### INITIALIZATION ###
  ######################
  # Values of hyper-parameters
  # The variables t, r, b, delta, c0 are not changed, 
  # only c is changed according to the desired effect
  t <- power(p, n)
  b <- 1/3
  if (is.null(r)) 
  {
    r <- 1e-1 + max(b - (t * delta / 2), 1e-1)
  }
  g <- c0 * (p^(2 * (1 + delta)))
  k <- c * (n^r)
  
  count <- rep(0, group_count)
  Groups <- 1:group_count
  group_freq <- rep(1, group_count)
  group_elements <- mclapply(Groups, function(x) { #provides the group identities
    which(var2group == x)
  })
  group_potential <- rep(0, group_count)
  
  GSIS_0 <- numeric(group_count) # GSIS with d = 0
  for (i in Groups) {
    G_k = group_elements[[i]]
    X_group = X[,G_k]
    y_g <- X_group  %*% 
            (pseudoinverse(t(X_group) %*% X_group + 1e-3 * 
                             diag(1, length(G_k))) %*% (YTX[1, G_k]))
    GSIS_0[i] <- max(1e-15, as.numeric(t(y_g) %*% y_g))
  }
  group_potential <- GSIS_0 / max(GSIS_0)
  mode <- numeric()
  mode_val <- 0
  
  ###############################
  ### POSTERIOR PROBABILITIES ###
  ###############################
  ## The gp_posterior function returns the posterior probability of the given 
  ## model up to a proportionality constant under the  modified-g prior setup
  ## The input model is the gamma vector (using the notation of the paper), and
  ## the output is the posterior probability up to a constant of proportionality
  gp_posterior <- function(model) {
    cardinality <- length(model)
    R_squared <- 0
    prob <- (1 + g)^(-cardinality / 2)
    prob <- prob * (k^(-cardinality))
    if (cardinality != 0) {
      X.model <- X[, model]
      xtx <- t(X.model) %*% X.model # X'X
      tryCatch(
        {
          Px <- chol2inv(chol(xtx))  # (X'X)^{-1}
          ytx <- YTX[, model]        # (X'y)' 
          R_squared <- as.numeric(ytx %*% Px %*% t(t(ytx))) / norm_y
        },
        error = function(e) {
          prob <- 0
        }
      )
    }
    prob <- max(0, prob * ((1 - (g / (1 + g)) * 
                              R_squared)^(-n / 2)) / choose(p, cardinality))
    if (prob > mode_val) {
      mode <<- model
      mode_val <<- prob
      group_freq <<- rep(0, group_count)
      group_freq[unique(var2group[mode])] <<- 1
    }
    return(prob)
  }
  
  ## The gigp_posterior function returns the posterior probability of the given 
  ## model up to a proportionality constant under the group informed-g prior 
  ## setup. 
  ## The input model is the gamma vector (using the notation of the paper), and
  ## the output is the posterior probability up to a constant of proportionality
  gigp_posterior <- function(model) {
    cardinality <- length(model)
    R_squared <- 0
    model_groups <- var2group[model]
    permute <- order(model_groups, decreasing = FALSE)
    model <- model[permute]
    mg_count <- cumsum(rle(model_groups[permute])$lengths)
    prob <- g^(-cardinality / 2)
    if (cardinality != 0) {
      X.model = X[, model]
      xtx <- t(X.model) %*% X.model 
      j <- 1
      for (i in 1:length(mg_count)) {
        l <- mg_count[i]
        xtx[j:l, j:l] <- xtx[j:l, j:l] * (1 + 1 / g) # (X'X + G/g)
        j <- l + 1
      }
      tryCatch(
        {
          Px <- chol2inv(chol(xtx))
          ytx <- YTX[, model]
          R_squared <- as.numeric(ytx %*% Px %*% t(t(ytx))) / norm_y
          prob <- prob / sqrt(det(xtx))
        },
        error = function(e) {
          prob <- 0
        }
      )
    }
    prob <- max(0, prob * ((1 - R_squared)^(-n / 2)) * 
                  prod(group_sizes[model_groups] / (k * p) ))
    if (prob > mode_val) {
      mode <<- model
      mode_val <<- prob
      group_freq <<- rep(0, group_count)
      group_freq[unique(var2group[mode])] <<- 1
    }
    return(prob)
  }
  
  if (!grp_post) {
    posterior <- gp_posterior
  } else {
    posterior <- gigp_posterior
  }

  
  ########################
  ### HELPER FUNCTIONS ###
  ########################
  # Returns probability of choosing a move given the size of the current model
  pdf <- function(size, move) { #move 1 for addition, 0 for swap, -1 for removal
    if (move == 1) {
      if (size == 0) { return(1) }
      if (size >= n) { return(0) }
      return(1 / 3)
    }
    if (move == 0) {
      if (size == 0) { return(0) }
      if (size >= n) { return(1 / 2) }
      return(1 / 3)
    }
    if (move == -1) {
      if (size == 0) { return(0) }
      if (size >= n) { return(1 / 2) }
      return(1 / 3)
    }
  }
  
  # An old function, currently it works as follows:
  # Move refer to either add move (1) or remove move (-1)
  # if move = add (1), if the variate is already in the model return 0 else 1
  # if move = remove(-1), if the variate is in the model return 1 else return 0
  
  w <- function(var, block_size, move) 
  {
    var <- as.numeric(var)
    prob <- 0
    if (move == 1) {  
      prob <- (1 - var) }
    if (move == -1) {
      prob <- var }
    prob <- min(1, max(0, prob)) 
    return(c(prob, 1 - prob))
  }
  # Functions for executing the moves
  # Input of these functions are moves and current model gamma, and the 
  # functions return the proposal model
  add <- function(i, block_size, gamma) {
    foo <- NULL
    eta <- sample(c(1, 0), 1, prob = w(0, block_size, 1))
    if (eta == 1) {
      foo <- c(gamma, i)
    }
    return(foo)
  }
  remove <- function(i, block_size, gamma) {
    foo <- NULL
    eta <- sample(c(1, 0), 1, prob = w(1, block_size, -1))
    if (eta == 1) {
      foo <- gamma[gamma != i]
    }
    return(foo)
  }
  swap <- function(j, i, block_size, gamma) {
    foo <- NULL
    card <- length(gamma)
    eta_a <- sample(c(1, 0), 1, prob = w(0, block_size, 1))
    if (eta_a == 1) {
      foo <- c(gamma[gamma != i], j)
    }
    return(foo)
  }
  
    
  #################################
  ### Pushing the for loop to C ###
  #################################
  sourceCpp("./Hash.cpp", env = environment())
  
  ######################
  ### pMTM ALGORITHM ###
  ######################
  MOVES <- c(-1, 0, 1)
  run <- function(gamma, t) {
  #  print(gamma)
  # Setting up values for the current run
    # cardinality of the current model
    card <- length(gamma)
    # move indicates the chosen move based on the cardinality of the current model
    move <- sample(MOVES, 1, prob = c(pdf(card, -1), pdf(card, 0), pdf(card, 1)))
    
    # Code to select a group
    prob <- rep(0, group_count)
    prob <- group_potential
    prob <- prob / sum(prob)
    prob <- al * prob + (1 - al) * (1 / group_count)
    prob <- prob / sum(prob)
    current_block <- sample(Groups, 1, prob = prob) #GiVSA step
    block_size <- group_sizes[[current_block]]
    
    ############ Forward Neighborhood ##########
    N_f <- list()
    # Code for forward addition neighbourhood
    if (move == 1) {
      N_f <- mclapply(setdiff(group_elements[[current_block]], gamma), 
                      add, block_size = block_size, gamma = gamma)
      N_f <- Filter(length, N_f)
    }
    # Code for forward removal neighbourhood
    if (move == -1) {
      N_f <- mclapply(gamma, remove, block_size = block_size, gamma = gamma)
    }
    # Code for forward swap neighbourhood
    if (move == 0) {
      for (i in gamma) {
        eta_r <- sample(c(1, 0), 1, prob = w(1, block_size, -1))
        if (eta_r == 1) {
          temp <- mclapply(setdiff(group_elements[[current_block]], gamma), 
                           swap, i = i, block_size = block_size, gamma = gamma)
          temp <- Filter(length, temp)
          N_f <- append(temp, N_f)
        }
      }
    }
    if (length(N_f) == 0) {
      # Taking care of problematic corner cases
      return(gamma)
    }
    fwd_scores <- unlist(mclapply(N_f, FUN = posterior))
    fwd_scores[is.na(fwd_scores)] <- 0
    
    ############ Forward Proposal ##########
    # Choosing the proposal based on forward scores
    if (all(fwd_scores == 0)) {
      proposal <- sample(N_f, 1)[[1]]
    } else {
      proposal <- (sample(N_f, 1, prob = fwd_scores))[[1]]
    }
    
    ############ Backward Neighborhood ##########
    N_b <- list()
    # Creating backward neighbourhoods
    r <- 0
    removed_block <- 0
    rblock_size <- 0
    # backward addition neighbourhood
    if (move == -1) {
      # r is the removed covariate
      r <- setdiff(gamma, proposal)
      removed_block <- var2group[r]
      rblock_size <- group_sizes[removed_block]
      N_b <- mclapply(setdiff(group_elements[[removed_block]], proposal), 
                      add, block_size = rblock_size, gamma = proposal)
      N_b <- Filter(length, N_b)
    }
    # backward removal neighbourhood
    else if (move == 1) {
      N_b <- mclapply(proposal, remove, block_size = block_size, gamma = proposal)
    }
    # backward swap neighbourhood
    else if (move == 0) {
      r <- setdiff(gamma, proposal)
      removed_block <- var2group[r]
      rblock_size <- group_sizes[var2group[r]]
      for (i in proposal) {
        eta_r <- sample(c(1, 0), 1, prob = w(1, rblock_size, -1))
        if (eta_r == 1) {
          temp <- mclapply(setdiff(group_elements[[removed_block]], proposal), 
                           swap, i = i, block_size = rblock_size, gamma = proposal)
          temp <- Filter(length, temp)
          N_b <- append(N_b, temp)
        }
      }
    }
    
    if (!(list(gamma) %in% N_b)) {
      N_b <- append(N_b, list(gamma))
    }
    bwd_scores <- unlist(mclapply(N_b, FUN = posterior))
    
    # Calculating probability of acceptance
    alpha <- 1
    if (move == 0) {
      alpha <- (prob[removed_block] * w(0, rblock_size, 1)[1] * 
                  w(1, rblock_size, -1)[1] * 
                  sum(fwd_scores)) / (prob[current_block] * 
                  w(0, block_size, 1)[1] * w(1, block_size, -1)[1] * 
                                        sum(bwd_scores))
    }
    if (move == -1) {
      alpha <- (pdf(card - 1, 1) * prob[removed_block] * 
                  w(0, rblock_size, 1)[1] * sum(fwd_scores)) / (pdf(card, -1) * 
                                    w(1, block_size, -1)[1] * sum(bwd_scores))
    }
    if (move == 1) {
      alpha <- (pdf(card + 1, -1) * w(1, block_size, -1)[1] * 
                  sum(fwd_scores)) / (prob[current_block] * pdf(card, 1) * 
                                      w(0, block_size, 1)[1] * sum(bwd_scores))
    }
    alpha <- min(1, alpha)
    
    # Adding the new sample to the chain
    u <- runif(1)
    if (u <= alpha) {
      gamma <- proposal
      # accepted <- TRUE
    }
    return(gamma)
  }
  # The for loop was pushed to C++ for faster computation
  samples <- list()
  start <- numeric()
  samples <- mcmc(start, run, Iter, burnin)
  
  # res contains all the models visited by the chain and the model with highest posterior probability
  res <- list()
  res$samples <- samples
  res$mode <- mode
  return(res)
}