#' Construct variance
#' @keywords internal
#' 
process_variance <- function(num_i, time_new, bio_i, data_predict_all, 
                             long_fit_all, time_variable) {
  
  #LME model fitting
  lfit = long_fit_all$lfit
  #variance-covariance matrix
  Sigma = long_fit_all$Sigma_fit
  long_sub_random = long_fit_all$long_sub_random
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all$long_sub_random[[1]], "|")[[2]])[2]
  ### event type variable name
  #event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  
  # data.long is a list containing all your input data frames
  # data.long is an input list from user
  # data.long must be a list, it can contain n data frames and each element contains one biomarker
  # or data.long can be a list and only contain one data matrix, all biomarkers are contained
  data.long <- data_predict_all
  
  # Convert 'data.long' to a list if it is not a list
  if (!is.list(data.long) || is.data.frame(data.long)) {
    data.long <- list(data.long)
    data.long <- rep(data.long, each = n_longitudinal)
  }
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  #residual variance
  sigma.longitudinal = c()
  for(i in 1:n_longitudinal){
    sigma.longitudinal[i] = lfit[[i]]$sigma
  }
  
  # Initialize lists to store the results for each subject 'num_i'
  # rep_num_i_list will store the repeated ones for each data frame.
  # data_num_i_list will store the filtered data for each patient num_i.
  rep_num_i_list <- list()
  data_num_i_list <- list()
  
  # Iterate over each data frame for different biomarkers
  for (i in 1:n_longitudinal) {
    df <- data.long[[i]]
    # Extract data for patient ID of 'num_i',  where 'num' equals 'num_i'
    selected_data <- df[which(as.character(df[[num]]) == as.character(num_i)), , drop = FALSE]
    
    # the biomarker used to predict,
    if(i == bio_i){
      selected_data <- rbind(selected_data, selected_data[nrow(selected_data), ])
      #replace time variable with predict time
      selected_data[time_variable][nrow(selected_data),] = time_new
    }
    # Store the row length of patient ID of 'num_i', in a vector of repeated 1, 
    # Store in the list for different biomarkers
    rep_num_i_list[[i]] <- rep(1, length(unlist(selected_data[time_variable])))
    # Store the filtered patient ID of 'num_i' data with all variables in the list
    data_num_i_list[[i]] <- selected_data
  }
  #if all biomarkers contained in one data frame
  if(length(data.long) == 1){
    for (i in 1:n_longitudinal) {#
      rep_num_i_list[[i]] <- rep_num_i_list[[1]]
      data_num_i_list[[i]] <- data_num_i_list[[1]]
    }
  }
  
  # Use lapply to check the length of each element, 
  # and then use any to determine whether there is an element with a length of 0
  if(any(sapply(data_num_i_list, nrow) == 0)) {
    return(NA)
  }
  
  A_i <- random_effects_design(data_num_i_list, long_sub_random, Sigma)

  Sigma_vector = c()
  for(i in 1:n_longitudinal){
    Sigma_vector = c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
  }
  ### diag(Sigma_vector) alone is unsafe when Sigma_vector has length 1
  ### (e.g. a patient with a single longitudinal observation to condition
  ### on): diag() then treats that single number as a matrix *dimension*
  ### and returns an n x n identity matrix instead of the intended 1 x 1
  ### diagonal matrix, corrupting Sigma_all's dimensions downstream (a
  ### classic base R diag() gotcha -- see ?diag). Passing the length
  ### explicitly avoids the ambiguity for every length, including 1.
  Sigma_all =  A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector, length(Sigma_vector))
  #Sigma_all <- diag(diag(Sigma_all))
  
  return(Sigma_all) # Return the computed Sigma_all for this iteration
}

#' Multivariate normal log-density for many means sharing one covariance
#' @description Equivalent to calling
#' \code{mvtnorm::dmvnorm(x, mean = means[[k]], sigma = sigma, log = TRUE)}
#' once for every element of \code{means}, but factors \code{sigma} only
#' once. \code{conditionalYTBio()}/\code{conditionalYDTBio()} evaluate the
#' density at every grid point \code{l_i} for the same patient, where only
#' the mean changes (the survival time enters the fixed effects, not the
#' covariance), so re-running the Cholesky decomposition inside
#' \code{dmvnorm()} at each grid point repeated an \code{O(N^3)} step
#' (\code{N} = the patient's total number of observations across all
#' biomarkers) \code{length(l_i)} times. The arithmetic below is the same
#' as \code{mvtnorm::dmvnorm()}'s (including its symmetry check and its
#' \code{-Inf}/\code{Inf} result when \code{sigma} is not positive
#' definite), so results agree with it exactly.
#' @param x Matrix with one row per point at which to evaluate the density.
#' @param means List of mean vectors, each of length \code{ncol(x)}.
#' @param sigma Covariance matrix shared by every element of \code{means}.
#' @return A list the same length as \code{means}; element \code{k} is the
#' vector of log-densities of the rows of \code{x} under \code{means[[k]]}.
#' @keywords internal
dmvnorm_shared_sigma <- function(x, means, sigma) {
  if (is.vector(x)) x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (p != ncol(sigma)) stop("x and sigma have non-conforming size")
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE))
    stop("sigma must be a symmetric matrix")
  dec <- tryCatch(base::chol(sigma), error = function(e) e)
  if (!inherits(dec, "error")) {
    log_const <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi)
  }
  lapply(means, function(mean) {
    mean <- c(mean)
    if (length(mean) != p) stop("x and mean have non-conforming size")
    if (inherits(dec, "error")) {
      x.is.mu <- colSums(t(x) != mean) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
    } else {
      tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
      logretval <- log_const - 0.5 * colSums(tmp^2)
    }
    names(logretval) <- rownames(x)
    logretval
  })
}
