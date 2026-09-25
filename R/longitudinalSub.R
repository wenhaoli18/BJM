#' The process involves estimating parameters for a multivariate linear mixed-effects 
#' model, which simultaneously analyzes multiple dependent variables that may be 
#' correlated. This approach incorporates both fixed effects, which are consistent 
#' across the population, and random effects, accounting for variations within 
#' groups or subjects. By fitting this model, one can assess the influence of 
#' predictor variables on several longitudinal outcomes while considering the inherent 
#' variability in the data due to random effects.
#' 
#' @param data_fit_all This process requires a set of \code{data.frame} objects 
#' designated for model fitting, with each \code{data.frame} representing a 
#' separate longitudinal outcome. These \code{data.frame} objects must include the 
#' variables identified in \code{long_sub_fixed} and \code{long_sub_random}. 
#' The use of a \code{list} arrangement facilitates the inclusion of 
#' various longitudinal outcomes, which may adhere to different measurement protocols. 
#' When all longitudinal outcomes are measured at the same time points for every patient, 
#' a single \code{data.frame} object can be in a list. 
#' It is assumed that every \code{data.frame} is organized in a long format.
#' 
#' @param long_sub_fixed This refers to a collection of formulas detailing the
#' fixed effects portion for each longitudinal outcome. On the left side of each formula,
#' the response variable is defined, while the right side outlines
#' the fixed effect terms. Should only a single formula be provided - whether
#' as a list with one item or as a standalone formula - it is inferred that
#' a conventional univariate joint model is being constructed.
#' Terms whose basis/contrasts depend on the data they are computed from --
#' \code{poly()} in its default orthogonal mode, \code{splines::ns()}/
#' \code{splines::bs()}, and \code{factor()} -- trigger a warning, because
#' \code{dynamicPrediction()}/\code{dynamicPredictionBio()} rebuild the
#' design matrix from a small, patient-specific slice of data at every point
#' on the prediction grid, so the basis recomputed at prediction time can
#' silently disagree with the one used to fit the model (or fail outright
#' with too few distinct values). Prefer \code{poly(..., raw = TRUE)},
#' \code{I(x^2)}, \code{log()}, \code{sqrt()}, or other terms that do not
#' depend on the surrounding data.
#' 
#' @param long_sub_random A list of one-sided formulas that define the model for the 
#' random effects of each longitudinal outcome. 
#' The number of items in this \code{list} should match the length of \code{formLongFixed}.
#' 
#' @return An object of class \code{"longitudinalSub.BJM"}, a named list with elements:
#' \describe{
#'   \item{lfit}{A list of fitted univariate linear mixed models, one per longitudinal
#'   outcome, each obtained via \code{\link[nlme]{lme}}.}
#'   \item{Sigma_fit}{The estimated variance-covariance matrix of the random effects
#'   in the multivariate linear mixed model.}
#'   \item{long_sub_fixed}{The \code{long_sub_fixed} argument, as supplied.}
#'   \item{long_sub_random}{The \code{long_sub_random} argument, as supplied.}
#'   \item{xlevels}{A list, one element per longitudinal outcome, of the
#'   factor levels observed in the full training data for that outcome.
#'   Used internally at prediction time so that \code{poly()}/
#'   \code{splines::ns()}/\code{splines::bs()}/\code{factor()} terms in
#'   \code{long_sub_fixed} reuse the basis/contrasts fit at training time
#'   instead of recomputing one from a small, patient-specific slice of
#'   data.}
#' }
#' 
#' @examples 
#' \donttest{
#' 
#' long_sub_fixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#' 
#' long_sub_random =list(
#'   "long1" =  ~ year| id,   
#'   "long2" =  ~ year| id,    
#'   "long3" =  ~ year| id,    
#'   "long4" =  ~ year| id,    
#'   "long5" =  ~ year| id,    
#'   "long6" =  ~ year| id)
#' 
#' survival_variable_all = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#' 
#' survival_trans_function = list(
#'   fun1 = function(x){abs(x - 1)}, 
#'   fun2 = function(x){abs(x - 3)}, 
#'   fun3 = function(x){abs(x - 5)}, 
#'   fun4 = function(x){abs(x - 7)}
#' )
#' 
#' # Complete case analysis
#' data_fit_all = list()
#' for(i in seq_len(length(long_sub_fixed))){
#'   data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#' 
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
#'
#' # poly() in its default orthogonal mode, splines::ns()/bs(), and
#' # factor() trigger a warning (see the long_sub_fixed argument above),
#' # but are still safe to use: the terms/xlevels/contrasts fit on the
#' # full training data are cached and reused at prediction time, instead
#' # of being recomputed from each patient's small per-prediction slice.
#' long_fit_poly = longitudinalSub(
#'   pbc3[pbc3$status3 == 1, ],
#'   serBilir ~ year + poly(age, 2) + factor(sex) + years,
#'   ~ year | id)
#'
#' }
#'
#' @export
longitudinalSub <- function(data_fit_all, long_sub_fixed, long_sub_random) {
  long_sub_fixed_check <- if (is.list(long_sub_fixed)) long_sub_fixed else list(long_sub_fixed)
  long_sub_random_check <- if (is.list(long_sub_random)) long_sub_random else list(long_sub_random)
  assert_all_formulas(long_sub_fixed_check, "long_sub_fixed")
  assert_all_formulas(long_sub_random_check, "long_sub_random")
  warn_unsafe_formula_terms(long_sub_fixed_check, "long_sub_fixed")
  if (length(long_sub_fixed_check) != length(long_sub_random_check)) {
    stop(sprintf(
      "`long_sub_fixed` has %d element(s) but `long_sub_random` has %d; they must describe the same number of longitudinal outcomes.",
      length(long_sub_fixed_check), length(long_sub_random_check)
    ), call. = FALSE)
  }

  if (!is.list(long_sub_fixed)) {
    long_sub_fixed <- list(long_sub_fixed)
    long_sub_random <- list(long_sub_random)
    M <- 1
  } else {
    ### number of biomarkers
    M <- length(long_sub_fixed)
  }

  assert_data_list(data_fit_all, "data_fit_all", M, allow_bare_df = TRUE)

  ### Convert 'data.long' to a list if it is not a list
  if (!is.list(data_fit_all) || is.data.frame(data_fit_all)) {
    data_fit_all <- list(data_fit_all)
    data_fit_all <- rep(data_fit_all, each = M)
  }

  ### patient id indicator
  id <- as.character(nlme::splitFormula(long_sub_random[[1]], "|")[[2]])[2]
  for (m in seq_len(M)) {
    assert_vars_in_data(unique(c(all.vars(long_sub_fixed[[m]]), all.vars(long_sub_random[[m]]), id)),
                         data_fit_all[[m]],
                         sprintf("long_sub_fixed[[%d]]/long_sub_random[[%d]]", m, m),
                         sprintf("data_fit_all[[%d]]", m))
  }
  ### number of patients, should be the same for each biomarker list
  ### changed and add unlist
  n <- length(unlist(unique(data_fit_all[[1]][, id])))
  
  lfit <- list()
  lfit_0 <- list()
  mf.fixed <- list()
  xlevels <- list()
  yik <- list()
  Xik <- list()
  nk <- vector(length = M)
  Xik.list <- list()
  nik.list <- list()
  Zik <- list()
  Zik.list <- list()
  
  ### different biomarkers have different missing samples, we have to get the intersect of them
  unique_num = list()
  for (m in 1:M) {
    data.fit.one = data_fit_all[[m]]
    ##exclude NA 
    data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(long_sub_fixed[[m]])])) )),]
    unique_num[[m]] = unique(unlist(data.fit.one[id])) 
  }
  
  if(M != 1){
    ##initiate a vector
    all_biomarker_num <- unlist(unique_num[[1]])
    ## get interaction for different biomarkers
    for (vec in unique_num[-1]) {  
      all_biomarker_num <- intersect(all_biomarker_num, unlist(vec))
    }
  }else{
    all_biomarker_num <- unique_num[[1]]
  }
  
  for (m in 1:M) {
    data.fit.one = data_fit_all[[m]]
    #ctrl <- lmeControl(1000, 1000, opt='optim')
    # List of m separate longitudinal model fits
    lfit[[m]] <- nlme::lme(fixed = long_sub_fixed[[m]], random = long_sub_random[[m]],
                           data = data.fit.one, method = "ML",
                           control = nlme::lmeControl(opt = "optim"), na.action = na.omit)
    lfit[[m]]$call$fixed <- eval(lfit[[m]]$call$fixed)
    
    ##exclude NA 
    data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(formula(lfit[[m]]))])) )),]
    
    ##interaction among different biomarkers
    data.fit.one = data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num),]
    
    # Model frames
    mf.fixed[[m]] <- model.frame(lfit[[m]]$terms,
                                 data.fit.one[, all.vars(long_sub_fixed[[m]])])

    ### factor levels observed in the full training data, cached so that
    ### prediction-time code can rebuild a model.frame() from a small,
    ### patient-specific slice of data without silently dropping levels that
    ### happen not to appear in that slice (which would otherwise make
    ### factor() error with "contrasts can be applied only to factors with 2
    ### or more levels", or -- for poly()/splines::ns()/bs() -- recompute a
    ### different basis than the one the model was fit with).
    xlevels[[m]] <- .getXlevels(lfit[[m]]$terms, mf.fixed[[m]])

    # Longitudinal outcomes by using "model.response" to get the response variable
    yik[[m]] <- by(model.response(mf.fixed[[m]], "numeric"), droplevels(data.fit.one[, id]), as.vector)
    
    # X design matrix, fixed effects design matrix
    Xik[[m]] <- data.frame("id2" = droplevels(data.fit.one[, id]),
                           model.matrix(long_sub_fixed[[m]], data.fit.one))
    
    # n_k (number of observations per each m)
    nk[m] <- nrow(Xik[[m]])
    
    # X design matrix (list)
    #Xik.list[[m]] <- by(Xik[[m]], Xik[[m]][id], function(u) {
    Xik.list[[m]] <- by(Xik[[m]], Xik[[m]]$id2, function(u) {
      as.matrix(u[, -1])
    })
    
    # number of repeated measurements for each subject
    #nik.list[[m]] <- by(Xik[[m]], Xik[[m]][id], nrow)
    nik.list[[m]] <- by(Xik[[m]], Xik[[m]]$id2, nrow)
    
    # Z design matrix, random effects design matrix
    ffk <- nlme::splitFormula(long_sub_random[[m]], "|")[[1]]
    Zik[[m]] <- data.frame("id2" = droplevels(data.fit.one[, id]), model.matrix(ffk, data.fit.one))
    
    # Z design matrix (list), list by subjects
    #Zik.list[[m]] <- by(Zik[[m]], c(unlist(Zik[[m]][id])), function(u) {
    #  as.matrix(u[, -1])
    #})
    Zik.list[[m]] <- by(Zik[[m]], c(Zik[[m]]$id2), function(u) {
      as.matrix(u[, -1])
    })
    
  }
  
  # Flatten lists to length = n
  yi <- sapply(names(yik[[1]]), function(i) {
    unlist(lapply(yik, "[[", i))
  },
  USE.NAMES = TRUE, simplify = FALSE)
  
  Xi <- sapply(names(Xik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Xik.list, "[[", i)))
  },
  USE.NAMES = TRUE, simplify = FALSE)
  
  Zi <- sapply(names(Zik.list[[1]]), function(i) {
    as.matrix(Matrix::bdiag(lapply(Zik.list, "[[", i)))
  },
  USE.NAMES = TRUE, simplify = FALSE)
  
  Zit <- lapply(Zi, t)
  
  # t(X) %*% X [summed over i and inverted]
  XtXi <- lapply(Xi, crossprod)
  XtX.inv <- solve(Reduce("+", XtXi))
  
  # t(X) %*% y [for each i]
  Xtyi <- mapply(function(x, y) {
    crossprod(x, y)
  },
  x = Xi, y = yi,
  SIMPLIFY = FALSE)
  
  # t(X) %*% Z [for each i]
  XtZi <- mapply(function(x, z) {
    crossprod(x, z)
  },
  x = Xi, z = Zi,
  SIMPLIFY = FALSE)
  
  nik <- sapply(names(nik.list[[1]]), function(i) {
    unlist(lapply(nik.list, "[[", i))
  },
  USE.NAMES = TRUE, simplify = FALSE)
  
  # Number of fixed effects
  p <- sapply(1:M, function(i) {
    ncol(Xik[[i]]) - 1
  })
  # Number of random effects
  r <- sapply(1:M, function(i) {
    ncol(Zik[[i]]) - 1
  })
  
  l <- list(yi = yi, Xi = Xi, Zi = Zi, Zit = Zit, nik = nik, yik = yik,
            Xik.list = Xik.list, Zik.list = Zik.list, XtX.inv = XtX.inv,
            Xtyi = Xtyi, XtZi = XtZi, p = p, r = r, M = M, n = n, nk = nk)
  
  #variance-covariance matrices without correlation
  #bdiag stands for block diagonal, and it constructs a block diagonal matrix f
  D <- Matrix::bdiag(lapply(lfit, function(u) matrix(nlme::getVarCov(u),
                                                     dim(nlme::getVarCov(u)))))
  D <- as.matrix(D)
  D.names <- c()
  for (m in 1:M) {
    D.names.k <- paste0(rownames(nlme::getVarCov(lfit[[m]])), "_", m)
    D.names <- c(D.names, D.names.k)
  }
  rownames(D) <- colnames(D) <- D.names
  
  beta.1 <- do.call("c", lapply(lfit, fixef))
  names(beta.1) <- paste0(names(beta.1), "_", rep(1:M, p))
  
  sigma2 <- unlist(lapply(lfit, function(u) u$sigma))^2
  
  result.solve.D <- tryCatch({
    solve(D)
    # code that may produce an error
  }, error = function(e) {
    # code to handle the error, such as printing a message
    message("An error occurred: ", conditionMessage(e))
    NULL  # return NULL to indicate that an error occurred
  })
  
  if (!is.null(result.solve.D)) {
    # code to execute if no error occurs
    # use the 'result' variable here, if necessary
    # ...
    out <- longitudinalSubVar(thetaLong = list("beta" = beta.1, "D" = D, "sigma2" = sigma2),
                              l = l, tol.em = 1e-04, verbose = FALSE)
  } else {
    # code to execute if an error occurs
    # ...
    out <- longitudinalSubVar(thetaLong = list("beta" = beta.1, "D" = diag(1, dim(D)[1]), "sigma2" = sigma2),
                              l = l, tol.em = 1e-04, verbose = FALSE)
  }
  ###**** need to changed back to Sigma_fit = out$D
  Sigma_fit = out$D
  #Sigma_fit = matrix(0,3,3)
  #Sigma_fit = diag(diag(out$D))
  
  long_fit_all = list(lfit = lfit, Sigma_fit = Sigma_fit,
                       long_sub_fixed = long_sub_fixed, long_sub_random = long_sub_random,
                       xlevels = xlevels)

  class(long_fit_all) <- "longitudinalSub.BJM"
  return(long_fit_all)
}
