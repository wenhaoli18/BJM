#' Construct variance
#' @keywords internal
#' 
process_variance <- function(num_i, time_new, bio_i, data.predict.all, 
                             long_fit_all, time_variable) {
  
  #LME model fitting
  lfit = long_fit_all[[1]]
  #variance-covariance matrix
  Sigma = long_fit_all[[2]]
  #patient ID
  num <- as.character(nlme::splitFormula(long_fit_all[[4]][[1]], "|")[[2]])[2]
  ### event type variable name
  #event_type_variable = as.character(formula(survival_fit_all[[4]])[[2]])
  
  #number of longitudinal biomarkers
  n_longitudinal <- length(lfit)  #length(data_num_i_list)
  
  # data.long is a list containing all your input data frames
  # data.long is an input list from user
  # data.long must be a list, it can contain n data frames and each element contains one biomarker
  # or data.long can be a list and only contain one data matrix, all biomarkers are contained
  data.long <- data.predict.all
  
  # Convert 'data.long' to a list if it is not a list
  if (!is.list(data.long)) {
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
    selected_data <- df[df[num] == num_i, ]
    
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
  
  A_i_ = list()
  ###random intercept or slope, depend on variance-covariance matrix
  if(dim(Sigma)[1] == n_longitudinal){
    ###random intercept
    for(i in 1:n_longitudinal){
      A_i_[[i]] = rbind(rep_num_i_list[[i]])
    }
  }else{
    ###random slope
    for(i in 1:n_longitudinal){
      A_i_[[i]] = rbind(rep_num_i_list[[i]], unlist(data_num_i_list[[i]][time_variable]))
    }
  }
  
  A_i <- matrix(0, nrow = sum(sapply(A_i_, ncol)), ncol = sum(sapply(A_i_, nrow)))
  length_A = rep(0, n_longitudinal + 1)
  for (i in 1:n_longitudinal) {
    length_A[i + 1] = length_A[i] + dim(A_i_[[i]])[2]
    if(dim(Sigma)[1] == n_longitudinal){
      ###random intercept
      A_i[c((length_A[i] + 1) : length_A[i + 1]), i]  <- t(A_i_[[i]])
    }else{
      ###random slope
      A_i[c((length_A[i] + 1) : length_A[i + 1]), (2*i-1):(2*i)]  <- t(A_i_[[i]])
    }
  }
  
  Sigma_vector = c()
  for(i in 1:n_longitudinal){
    Sigma_vector = c(Sigma_vector, rep(sigma.longitudinal[i]^2, dim(data_num_i_list[[i]])[1]))
  }
  Sigma_all =  A_i %*% Sigma %*% t(A_i) + diag(Sigma_vector)
  #Sigma_all <- diag(diag(Sigma_all))
  
  return(Sigma_all) # Return the computed Sigma_all for this iteration
}
