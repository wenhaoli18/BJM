#' Fit a multivariate longitudinal sub-model
#'
#' @description
#' Fits one linear mixed-effects model per longitudinal biomarker separately
#' (via \code{\link[nlme]{lme}}), keeping each biomarker's own fixed-effects
#' and residual-variance estimates from that separate fit. The \code{M}
#' separate fits are then combined into a single multivariate model by
#' re-estimating the full random-effects variance-covariance matrix jointly
#' across all \code{M} biomarkers, via an EM algorithm initialized from the
#' block-diagonal covariance implied by the separate fits -- so correlation
#' between biomarkers' random effects is captured, even though the fixed
#' effects and residual variances themselves are not re-estimated jointly
#' and remain exactly what each biomarker's own \code{lme()} fit produced.
#' This lets you assess the effect of predictor variables on several
#' longitudinal outcomes at once while accounting for both population-level
#' (fixed) effects and subject-level (random) variability, and for how the
#' outcomes covary within a subject. The result is one of the two sub-models
#' -- together with \code{\link{survivalSub}} -- that
#' \code{\link{predictRisk}}/\code{\link{dynamicPredictionBio}}
#' combine to produce dynamic risk/biomarker predictions.
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
#' \code{splines::bs()}, and \code{factor()} -- are fine in a continuous
#' biomarker's \code{long_sub_fixed} formula: the basis and factor levels the
#' model was fit with are reused when estimating \code{Sigma_fit} and at
#' prediction time. In \code{long_sub_random}, or in an ordinal biomarker's
#' \code{long_sub_fixed}, they trigger a warning, because there the design
#' matrix is still rebuilt from the data at hand (a single patient's rows at
#' prediction time), so the basis can silently disagree with the one used to
#' fit the model. Prefer \code{poly(..., raw = TRUE)}, \code{I(x^2)},
#' \code{log()}, \code{sqrt()}, or other terms that do not depend on the
#' surrounding data there.
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
#' # factor() are safe in a continuous biomarker's long_sub_fixed: the
#' # terms/xlevels/contrasts fit on the full training data are cached and
#' # reused when estimating Sigma_fit and at prediction time, instead of
#' # being recomputed from each patient's small per-prediction slice.
#' long_fit_poly = longitudinalSub(
#'   pbc3[pbc3$status3 == 1, ],
#'   serBilir ~ year + poly(age, 2) + factor(sex) + years,
#'   ~ year | id)
#'
#' # A few more nonlinear-term styles. None of these would trigger the
#' # warning even in long_sub_random, because their basis doesn't depend on the surrounding
#' # data at all -- there's simply nothing that could disagree between the
#' # full training data and the small per-patient slice used at prediction
#' # time.
#' data_fit_extra_terms = pbc3[pbc3$status3 == 1, ]
#' # A transformed-time column can also be built once with survivalTrans()
#' # and merged onto the fitting data ahead of time -- by the time
#' # long_sub_fixed sees it, it is just an ordinary numeric column, same as
#' # any of the other terms below.
#' data_fit_extra_terms$Tyears1 = survivalTrans(c(1, 3, 5, 7))$survival_trans_function[[1]](
#'   data_fit_extra_terms$year)
#'
#' long_fit_nonlinear = list(
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ year + I(year^2) + age + sex,
#'                    ~ year | id),
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ year + I(year^2) * age + sex,
#'                    ~ year | id),
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ year + log(year + 1) + age + sex,
#'                    ~ year | id),
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ year + sqrt(year) + age + sex,
#'                    ~ year | id),
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ poly(year, 2, raw = TRUE) + age + sex,
#'                    ~ year | id),
#'   longitudinalSub(data_fit_extra_terms, serBilir ~ year + Tyears1 + age + sex,
#'                    ~ year | id)
#' )
#'
#' }
#'
#' @param biomarker_type Optional character vector, one entry per
#' longitudinal outcome, each either \code{"continuous"} or \code{"ordinal"}.
#' When supplied, it always takes priority. When \code{NULL} (the default),
#' each biomarker's type is auto-detected from its own response column in
#' \code{data_fit_all}: a \code{factor}/\code{ordered factor} response is
#' treated as \code{"ordinal"}, anything else as \code{"continuous"}. If
#' every biomarker is continuous (the original use case), fitting proceeds
#' exactly as before with no behavior change whatsoever. If at least one
#' biomarker is ordinal, that biomarker is instead fit with
#' \code{ordinal::clmm()} (a probit cumulative link mixed model, requiring
#' the optional \pkg{ordinal} package) and folded into the shared
#' random-effects covariance matrix via a Gaussian-copula extension of the
#' EM algorithm -- see \code{longitudinalSubCopula()} for implementation
#' details.
#'
#' @export
longitudinalSub <- function(data_fit_all, long_sub_fixed, long_sub_random, biomarker_type = NULL) {
  long_sub_fixed_check <- if (is.list(long_sub_fixed)) long_sub_fixed else list(long_sub_fixed)
  long_sub_random_check <- if (is.list(long_sub_random)) long_sub_random else list(long_sub_random)
  assert_all_formulas(long_sub_fixed_check, "long_sub_fixed")
  assert_all_formulas(long_sub_random_check, "long_sub_random")
  if (length(long_sub_fixed_check) != length(long_sub_random_check)) {
    stop(sprintf(
      "`long_sub_fixed` has %d element(s) but `long_sub_random` has %d; they must describe the same number of longitudinal outcomes.",
      length(long_sub_fixed_check), length(long_sub_random_check)
    ), call. = FALSE)
  }
  M <- length(long_sub_fixed_check)

  assert_data_list(data_fit_all, "data_fit_all", M, allow_bare_df = TRUE)
  data_fit_all_norm <- if (!is.list(data_fit_all) || is.data.frame(data_fit_all)) {
    rep(list(data_fit_all), M)
  } else {
    data_fit_all
  }

  biomarker_type_resolved <- resolve_biomarker_type(biomarker_type, data_fit_all_norm, long_sub_fixed_check, M)

  ### Hard requirement: an all-continuous input must behave exactly as
  ### before. longitudinalSubGaussian() is a verbatim copy of this
  ### function's original body, called here with the *pristine, unmodified*
  ### arguments as supplied by the caller -- not the normalized-to-list
  ### copies above -- so the all-continuous code path is byte-for-byte
  ### identical to what ran prior to biomarker_type existing at all.
  if (all(biomarker_type_resolved == "continuous")) {
    return(longitudinalSubGaussian(data_fit_all, long_sub_fixed, long_sub_random))
  }

  ordinal_fixed <- long_sub_fixed_check
  ordinal_fixed[biomarker_type_resolved != "ordinal"] <- list(~ 1)
  warn_unsafe_formula_terms(ordinal_fixed, "long_sub_fixed")
  warn_unsafe_formula_terms(long_sub_random_check, "long_sub_random")

  id <- as.character(nlme::splitFormula(long_sub_random_check[[1]], "|")[[2]])[2]
  for (m in seq_len(M)) {
    assert_vars_in_data(unique(c(all.vars(long_sub_fixed_check[[m]]), all.vars(long_sub_random_check[[m]]), id)),
                         data_fit_all_norm[[m]],
                         sprintf("long_sub_fixed[[%d]]/long_sub_random[[%d]]", m, m),
                         sprintf("data_fit_all[[%d]]", m))
  }

  ### biomarker_type = "ordinal" may be given for a response that is not a
  ### factor (e.g. a 0/1/2 score); clmm() needs a factor, which it used to
  ### fail on. A numeric (or logical) response is ordered by value; a
  ### character one becomes an unordered factor, so fit_marginal_ordinal()
  ### warns that its alphabetical order is taken as the category order.
  for (m in which(biomarker_type_resolved == "ordinal")) {
    resp_name <- all.vars(long_sub_fixed_check[[m]])[1]
    resp <- data_fit_all_norm[[m]][[resp_name]]
    if (is.numeric(resp) || is.logical(resp)) {
      data_fit_all_norm[[m]][[resp_name]] <- factor(resp, levels = sort(unique(resp)), ordered = TRUE)
    } else if (is.character(resp)) {
      data_fit_all_norm[[m]][[resp_name]] <- factor(resp)
    }
  }

  longitudinalSubCopula(data_fit_all_norm, long_sub_fixed_check, long_sub_random_check,
                         biomarker_type_resolved, M)
}

#' Fit a multivariate longitudinal sub-model (all-continuous biomarkers)
#'
#' @description Internal workhorse for \code{\link{longitudinalSub}} when
#' every biomarker is continuous. This is the original \code{longitudinalSub}
#' implementation, kept verbatim and called with the caller's pristine,
#' unmodified arguments, so that the all-continuous case behaves exactly as
#' it always has -- see \code{\link{longitudinalSub}} for the argument
#' documentation, and \code{longitudinalSubCopula()} for the mixed
#' continuous/ordinal extension used when at least one biomarker is
#' categorical.
#' @keywords internal
longitudinalSubGaussian <- function(data_fit_all, long_sub_fixed, long_sub_random) {
  long_sub_fixed_check <- if (is.list(long_sub_fixed)) long_sub_fixed else list(long_sub_fixed)
  long_sub_random_check <- if (is.list(long_sub_random)) long_sub_random else list(long_sub_random)
  assert_all_formulas(long_sub_fixed_check, "long_sub_fixed")
  assert_all_formulas(long_sub_random_check, "long_sub_random")
  warn_unsafe_formula_terms(long_sub_random_check, "long_sub_random")
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

    ### factor levels of the data lme() was fit on (before restricting to
    ### the subjects retained in the joint fit, which may lack some levels)
    xlevels[[m]] <- training_xlevels(lfit[[m]]$terms, data.fit.one, long_sub_fixed[[m]])
    
    ##exclude NA 
    data.fit.one <- data.fit.one[!as.logical(rowSums(data.frame(is.na(data.fit.one[all.vars(formula(lfit[[m]]))])) )),]
    
    ##interaction among different biomarkers
    data.fit.one = data.fit.one[unlist(data.fit.one[id]) %in% unlist(all_biomarker_num),]
    
    # Model frames
    mf.fixed[[m]] <- model.frame(lfit[[m]]$terms,
                                 data.fit.one[, all.vars(long_sub_fixed[[m]])],
                                 xlev = xlevels[[m]])

    ### factor levels observed in the full training data, cached so that
    ### prediction-time code can rebuild a model.frame() from a small,
    ### patient-specific slice of data without silently dropping levels that
    ### happen not to appear in that slice (which would otherwise make
    ### factor() error with "contrasts can be applied only to factors with 2
    ### or more levels", or -- for poly()/splines::ns()/bs() -- recompute a
    ### different basis than the one the model was fit with).

    # Longitudinal outcomes by using "model.response" to get the response variable
    yik[[m]] <- by(model.response(mf.fixed[[m]], "numeric"), factor(data.fit.one[[id]]), as.vector)
    
    # X design matrix, fixed effects design matrix
    ### built from lme()'s own terms (which carry the training-data basis of
    ### poly()/ns()/bs() via predvars) and levels, not re-derived from the
    ### formula on the retained subjects only -- that recomputed a different
    ### basis than the one beta was estimated on whenever subjects were
    ### excluded, corrupting the EM residuals behind Sigma_fit.
    Xik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]),
                           model.matrix(lfit[[m]]$terms, mf.fixed[[m]],
                                        contrasts.arg = lfit[[m]]$contrasts))
    
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
    Zik[[m]] <- data.frame("id2" = factor(data.fit.one[[id]]), model.matrix(ffk, data.fit.one))
    
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
            Xtyi = Xtyi, XtZi = XtZi, p = p, r = r, M = M, n = n_subjects(yi), nk = nk)
  
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
  
  ### starting value for the joint EM: the block-diagonal D from the
  ### separate per-biomarker fits, or the identity if that is singular
  ### (previously reported only as a bare message()).
  D_start <- tryCatch({ solve(D); D }, error = function(e) {
    warning(sprintf(paste0(
      "The random-effects covariance from the separate per-biomarker fits is singular (%s); ",
      "starting the joint EM from the identity matrix instead. Sigma_fit may be unreliable: ",
      "consider simplifying the random-effects structure."), conditionMessage(e)), call. = FALSE)
    diag(1, dim(D)[1])
  })
  out <- longitudinalSubVar(thetaLong = list("beta" = beta.1, "D" = D_start, "sigma2" = sigma2),
                            l = l, tol.em = 1e-04, verbose = FALSE)
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

#' Number of subjects the random-effects covariance is estimated from
#'
#' @description The EM update for \code{Sigma_fit} averages each subject's
#' \eqn{E[b b^T]} over the subjects that are actually in the fit -- those with
#' at least one complete observation of every biomarker (\code{yi} has one
#' element per such subject). It used to divide by the number of subjects
#' in the first biomarker's raw data instead, which underestimated
#' \code{Sigma_fit} (badly, since EM compounds it over iterations) whenever
#' some subject was excluded, e.g. one never measured for some biomarker.
#'
#' @param yi The per-subject list of stacked responses.
#' @return The number of subjects.
#' @keywords internal
n_subjects <- function(yi) length(yi)

#' Factor levels of the data a biomarker's sub-model was fit on
#'
#' @description Used by \code{longitudinalSub()} to cache, per biomarker,
#' the factor levels seen by \code{lme()} -- computed from the full fitting
#' data before it is restricted to the subjects retained in the joint fit,
#' which may lack some levels.
#'
#' @param terms_model The fitted model's \code{terms}.
#' @param data The data the model was fit on.
#' @param fixed_formula The fixed-effects formula (for its variables).
#' @return A named list of factor levels, as from \code{.getXlevels()}.
#' @keywords internal
training_xlevels <- function(terms_model, data, fixed_formula) {
  mf <- model.frame(terms_model, data[, all.vars(fixed_formula), drop = FALSE], na.action = na.omit)
  .getXlevels(terms_model, mf)
}

