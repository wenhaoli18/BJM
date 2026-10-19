#' Impute missing longitudinal biomarker values
#'
#' @description Fills in missing biomarker values before
#' \code{longitudinalSub()} is fit, assuming they are missing at random
#' (MAR): whether a value is missing depends only on what is observed.
#'
#' \strong{It is often not needed.} \code{longitudinalSub()} fits each
#' biomarker's linear mixed model by maximum likelihood on all of that
#' biomarker's observed values, which is already valid under MAR; only the
#' joint covariance of the random effects across biomarkers,
#' \code{Sigma_fit}, is estimated from the subjects with some value of
#' every biomarker. In simulations on \code{pbc3} (see the "Missing data
#' imputation" article on the package website), no imputation method made
#' the predictions closer to those from the complete data than this
#' complete-case analysis when visits were missing values at random; with
#' half of the subjects never measured for a biomarker, multiple
#' imputation with \code{method = "mixed"} was slightly closer. Imputation
#' can also help when the missingness depends on variables outside the
#' model, if they are added to the formulas used here.
#'
#' Three methods are available, selected with \code{method}:
#' \itemize{
#'   \item \code{"mixed"} (default) draws each subject's missing values from
#'   the multivariate linear mixed model \code{longitudinalSub()} fits, given
#'   everything observed for that subject -- all visits and all biomarkers
#'   -- with the model's parameters drawn from an approximation to their
#'   posterior for each completion (proper multiple imputation; see
#'   \code{mixed_model_imputations()}). It needs no extra package, keeps the
#'   variation between subjects, and was the only method in the simulations
#'   that left the random-effects variances and the event-time coefficients
#'   unbiased.
#'   \item \code{"miwae"} fits a MIWAE (Missing data Importance-Weighted
#'   AutoEncoder; Mattei & Frellsen, 2019) -- a variational autoencoder
#'   whose training objective only scores reconstruction on the observed
#'   entries -- jointly across all supplied biomarkers, and uses
#'   self-normalized importance resampling from the fitted decoder to draw
#'   one or more plausible completions of the missing cells.
#'   \item \code{"diffusion"} fits a conditional denoising diffusion
#'   probabilistic model (DDPM; Ho, Jain & Abbeel, 2020) over the same
#'   per-row biomarker vector, trained as in CSDI (Tashiro et al., 2021):
#'   in every mini-batch, each row's observed entries are split at random
#'   into targets and conditions, and the network learns to denoise the
#'   targets given the conditions' values, the covariates and the
#'   observation mask, so the loss only ever involves observed entries.
#'   Missing cells are then filled by running the reverse diffusion with
#'   the missing entries as targets and all observed ones as conditions;
#'   observed entries are returned unchanged.
#' }
#' The two deep generative methods need the \pkg{torch} package. They
#' complete each visit from that visit's covariates and observed biomarkers
#' (and, with \code{history = TRUE}, a summary of the subject's other
#' visits), and in the simulations they shrank the estimated variation
#' between subjects, which every prediction uses; with \code{history = TRUE}
#' they also biased the event-time coefficients when subjects lacked a
#' biomarker entirely. The completed data.frame(s) of any method can be
#' passed directly to \code{longitudinalSub()} in place of the original
#' \code{data_fit_all}.
#'
#' @details Only the biomarker (response) columns named on the left-hand
#' side of \code{long_sub_fixed} are imputed, and they must be numeric
#' (continuous): an ordinal biomarker cannot be imputed with this function,
#' and gives an error. Every covariate referenced on
#' the right-hand side of \code{long_sub_fixed}/\code{long_sub_random} is
#' assumed fully observed (an error is raised if any covariate has missing
#' values). The deep generative methods model
#' \code{p(biomarkers | covariates, missingness mask)} row by row (one row
#' per subject-visit; the covariates include the visit time
#' \code{time_variable} when \code{long_sub_fixed} uses it), and use
#' \code{long_sub_random} only for the id variable; \code{"mixed"} uses the
#' full random-effects structure.
#'
#' Rows are aligned across biomarkers by \code{(id, time_variable)}: when
#' \code{data_fit_all} is a list of different data.frames with different
#' visit schedules per biomarker (rather than one data.frame shared by
#' every biomarker), the union of every \code{(id, time_variable)} pair seen
#' in any of them is used, so a visit at which one biomarker was measured
#' but another was not still gets a value imputed for the missing one.
#'
#' Because the deep generative methods need enough data to fit
#' reliably; the \code{diagnostics} element of the return value also reports
#' a classical single regression imputation (each biomarker regressed on
#' the covariates, fit on its observed rows) as a lightweight baseline, so
#' the two can be compared on a given dataset before trusting the deep
#' model's completions.
#'
#' @param data_fit_all As in \code{\link{longitudinalSub}}: either a single
#' data.frame reused for every biomarker, or a list with one data.frame per
#' biomarker (long format, one row per subject-visit).
#' @param long_sub_fixed As in \code{\link{longitudinalSub}}: a formula, or
#' list of formulas (one per biomarker), whose left-hand side must be a bare
#' column name (no transformation) -- the column that gets imputed.
#' @param long_sub_random As in \code{\link{longitudinalSub}}: used only to
#' read off the subject id variable (\code{~ ... | id}).
#' @param time_variable Name of the visit-time column used, together with
#' the id variable, to align rows across biomarkers.
#' @param n_imputations Number of stochastic completions to draw. Also
#' controls how many completed datasets are returned when
#' \code{impute = "multiple"}.
#' @param latent_dim Dimension of the VAE's latent space. Only used when
#' \code{method = "miwae"}.
#' @param hidden_units Integer vector giving the hidden-layer sizes of the
#' network: the encoder (and, reversed, the decoder) when
#' \code{method = "miwae"}, or the denoising network's body when
#' \code{method = "diffusion"}.
#' @param epochs Number of training epochs: full-batch gradient steps for
#' \code{method = "miwae"}, passes over the data in mini-batches of 128
#' rows for \code{method = "diffusion"}. The default suits data sets of a
#' few hundred to a few thousand rows; many more epochs can overfit, which
#' shows as imputations that vary too little between draws (and so pooled
#' standard errors that are too small) -- in a simulation, MIWAE's spread
#' between draws fell from about the true residual SD at 300 epochs to 80\%
#' of it at 3000.
#' @param importance_samples Number of importance samples (\code{K}) used
#' both in the training objective and at imputation time. Only used when
#' \code{method = "miwae"}.
#' @param diffusion_steps Number of forward/reverse diffusion timesteps.
#' Only used when \code{method = "diffusion"}; more steps generally
#' improve sample quality at proportionally higher training/sampling cost.
#' @param method \code{"mixed"} (default), \code{"miwae"} or
#' \code{"diffusion"}; see Description. \code{latent_dim},
#' \code{hidden_units}, \code{epochs}, \code{importance_samples} and
#' \code{diffusion_steps} only apply to the deep generative methods.
#' @param impute Either \code{"single"} (default; returns one completed
#' \code{data_fit_all}, filling each missing cell with the mean of the
#' \code{n_imputations} draws -- a drop-in replacement for the original
#' \code{data_fit_all}) or \code{"multiple"} (also returns
#' \code{n_imputations} separately-drawn completed datasets in
#' \code{data_fit_all_list}, to fit \code{longitudinalSub()} once per
#' completion and pool the fits with \code{\link{poolLongitudinalSub}}).
#' Prefer \code{"multiple"} unless only a little data is missing. The mean
#' of several draws is less variable than the values it stands in for, so a
#' model fit to the \code{"single"} completion treats the imputed cells as
#' exactly known: its standard errors are too small, and its residual and
#' random-effects variances -- including \code{Sigma_fit}, which
#' \code{\link{predictRisk}} and \code{\link{predictLongitudinal}} use
#' directly -- tend to be underestimated, more so the larger the share of
#' imputed cells.
#' @param history For the deep generative methods: whether each visit's
#' inputs include a summary of the same subject's other visits (see
#' \code{history_features()}). Without it an imputed value does not carry
#' the subject's own level, which shrinks the estimated variation between
#' subjects; with it, that shrinkage is smaller when visits miss values at
#' random, but in simulations it biased the event-time coefficients when
#' subjects lacked a biomarker entirely. Defaults to \code{FALSE}.
#' @param seed Optional integer seed for reproducibility. It seeds both
#' \pkg{torch} and R's random number generator for the duration of the
#' call; R's random number state from before the call is restored
#' afterwards.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{data_fit_all}{A list with one completed data frame per
#'   biomarker (also when \code{data_fit_all} was a single data frame), with
#'   every missing biomarker cell filled. Each holds one row per
#'   (id, \code{time_variable}) pair seen for any biomarker, and only the
#'   columns its model uses: the id, \code{time_variable}, the biomarker and
#'   the covariates of its \code{long_sub_fixed} formula -- other columns of
#'   the input are not carried over. Ready to pass straight to
#'   \code{\link{longitudinalSub}} with the same formulas.}
#'   \item{data_fit_all_list}{When \code{impute = "multiple"}, a list of
#'   \code{n_imputations} completed versions of \code{data_fit_all}, each
#'   from its own draw; \code{NULL} when \code{impute = "single"}.}
#'   \item{diagnostics}{A named list, one element per biomarker (its
#'   \code{deep_generative} entry describes the completions of whichever
#'   \code{method} was used), each with
#'   the number of missing cells and observed-vs-imputed mean/sd for both
#'   the deep-generative and classical-baseline completions.}
#' }
#'
#' @references Mattei, P.-A. and Frellsen, J. (2019). MIWAE: Deep
#' Generative Modelling and Imputation of Incomplete Data Sets.
#' \emph{Proceedings of the 36th International Conference on Machine
#' Learning}, PMLR 97:4413-4423.
#' @references Ho, J., Jain, A., and Abbeel, P. (2020). Denoising Diffusion
#' Probabilistic Models. \emph{Advances in Neural Information Processing
#' Systems}, 33:6840-6851.
#' @references Tashiro, Y., Song, J., Song, Y., and Ermon, S. (2021). CSDI:
#' Conditional Score-based Diffusion Models for Probabilistic Time Series
#' Imputation. \emph{Advances in Neural Information Processing Systems},
#' 34:24804-24816.
#'
#' @examples
#' data(pbc3)
#' data_fit_all <- pbc3[pbc3$status3 == 1, ]
#'
#' # pbc3 has no missing serBilir/albumin values; delete some at random
#' # (MCAR, a special case of MAR) to have something to impute
#' set.seed(1)
#' n <- nrow(data_fit_all)
#' data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
#' data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA
#'
#' long_sub_fixed <- list(
#'   "long1" = serBilir ~ year + age + sex + years,
#'   "long2" = albumin ~ year + age + sex + years)
#' long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
#'
#' # multiple imputation from the mixed model, pooled with Rubin's rules
#' imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
#'                               time_variable = "year", impute = "multiple", seed = 1)
#' fits <- lapply(imputed$data_fit_all_list, longitudinalSub,
#'                long_sub_fixed = long_sub_fixed, long_sub_random = long_sub_random)
#' pooled <- poolLongitudinalSub(fits)
#'
#' \donttest{
#' # a deep generative method instead (needs torch)
#' if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
#'   imputed_miwae <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
#'                                       time_variable = "year", method = "miwae",
#'                                       epochs = 50, seed = 1)
#' }
#' }
#'
#' @export
imputeLongitudinal <- function(data_fit_all, long_sub_fixed, long_sub_random,
                                time_variable, n_imputations = 5,
                                latent_dim = 8, hidden_units = c(64, 32),
                                epochs = 300, importance_samples = 20,
                                diffusion_steps = 100,
                                method = c("mixed", "miwae", "diffusion"),
                                impute = c("single", "multiple"),
                                history = FALSE, seed = NULL) {
  method <- match.arg(method)
  impute <- match.arg(impute)

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

  if (!is.list(long_sub_fixed)) {
    long_sub_fixed <- list(long_sub_fixed)
    long_sub_random <- list(long_sub_random)
  }
  M <- length(long_sub_fixed)

  assert_data_list(data_fit_all, "data_fit_all", M, allow_bare_df = TRUE)
  if (!is.list(data_fit_all) || is.data.frame(data_fit_all)) {
    data_fit_all <- rep(list(data_fit_all), M)
  }

  assert_string(time_variable, "time_variable")
  assert_positive_integer(n_imputations, "n_imputations")
  assert_positive_integer(epochs, "epochs")
  assert_positive_integer(importance_samples, "importance_samples")
  assert_positive_integer(diffusion_steps, "diffusion_steps")
  assert_scalar_numeric(latent_dim, "latent_dim", positive = TRUE)
  if (!is.logical(history) || length(history) != 1 || is.na(history)) {
    stop("`history` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(hidden_units) || length(hidden_units) == 0 ||
      any(hidden_units != round(hidden_units)) || any(hidden_units < 1)) {
    stop("`hidden_units` must be a non-empty vector of positive integers.", call. = FALSE)
  }

  id <- as.character(nlme::splitFormula(long_sub_random[[1]], "|")[[2]])[2]

  for (m in seq_len(M)) {
    lhs <- long_sub_fixed[[m]][[2]]
    if (!is.symbol(lhs)) {
      stop(sprintf(
        "`long_sub_fixed[[%d]]` has a transformed left-hand side (%s); imputeLongitudinal() requires each response to be a bare column name.",
        m, deparse(lhs)
      ), call. = FALSE)
    }
    assert_vars_in_data(unique(c(all.vars(long_sub_fixed[[m]]), all.vars(long_sub_random[[m]]), time_variable)),
                         data_fit_all[[m]],
                         sprintf("long_sub_fixed[[%d]]/long_sub_random[[%d]]/time_variable", m, m),
                         sprintf("data_fit_all[[%d]]", m))
    ### rows are aligned across biomarkers by (id, time_variable), which must
    ### therefore identify a visit
    dup <- duplicated(data_fit_all[[m]][, c(id, time_variable)])
    if (any(dup)) {
      stop(sprintf(paste0(
        "data_fit_all[[%d]] has %d row(s) repeating an (%s, %s) pair; imputeLongitudinal() needs ",
        "one row per subject and visit time."), m, sum(dup), id, time_variable), call. = FALSE)
    }
  }

  resp_vars <- unname(vapply(long_sub_fixed, function(f) as.character(f[[2]]), character(1)))
  if (any(duplicated(resp_vars))) {
    stop(sprintf(
      "`long_sub_fixed` response variables must be distinct across biomarkers; found a repeated response: %s.",
      resp_vars[duplicated(resp_vars)][1]
    ), call. = FALSE)
  }
  ### the imputation model is Gaussian: a factor (e.g. ordinal) response
  ### used to be coerced to numeric, turning the whole column into NaN
  for (m in seq_len(M)) {
    if (!is.numeric(data_fit_all[[m]][[resp_vars[m]]])) {
      stop(sprintf(paste0(
        "Biomarker `%s` is not numeric; imputeLongitudinal() imputes continuous biomarkers only ",
        "(an ordinal biomarker cannot be imputed with it)."), resp_vars[m]), call. = FALSE)
    }
  }
  covariate_vars <- lapply(seq_len(M), function(m) setdiff(all.vars(long_sub_fixed[[m]]), resp_vars[m]))
  all_covariates <- setdiff(unique(unlist(covariate_vars)), c(id, time_variable))

  ### align every biomarker's data onto one (id, time_variable) row set,
  ### pulling each covariate from the first biomarker's data.frame that has
  ### it, so a covariate repeated identically across biomarkers (the common
  ### case) is not merged/duplicated once per biomarker
  keys <- unique(do.call(rbind, lapply(data_fit_all, function(d) d[, c(id, time_variable)])))
  keys <- keys[order(keys[[id]], keys[[time_variable]]), , drop = FALSE]
  rownames(keys) <- NULL

  aligned <- keys
  for (cov in all_covariates) {
    src_m <- which(vapply(covariate_vars, function(v) cov %in% v, logical(1)))[1]
    piece <- unique(data_fit_all[[src_m]][, c(id, time_variable, cov)])
    aligned <- merge(aligned, piece, by = c(id, time_variable), all.x = TRUE, sort = FALSE)
  }
  for (m in seq_len(M)) {
    piece <- unique(data_fit_all[[m]][, c(id, time_variable, resp_vars[m])])
    aligned <- merge(aligned, piece, by = c(id, time_variable), all.x = TRUE, sort = FALSE)
  }
  aligned <- aligned[order(aligned[[id]], aligned[[time_variable]]), , drop = FALSE]
  rownames(aligned) <- NULL

  if (length(all_covariates) > 0) {
    cov_na <- vapply(all_covariates, function(cov) anyNA(aligned[[cov]]), logical(1))
    if (any(cov_na)) {
      stop(sprintf(
        "imputeLongitudinal() only imputes the longitudinal outcome columns (%s); covariate(s) %s have missing values, which is not supported. Impute or drop missing covariates first.",
        paste(resp_vars, collapse = ", "), paste(all_covariates[cov_na], collapse = ", ")
      ), call. = FALSE)
    }
  }

  resp_mat <- as.matrix(aligned[, resp_vars, drop = FALSE])
  mode(resp_mat) <- "double"
  mask <- !is.na(resp_mat)
  if (!any(!mask)) {
    stop("`data_fit_all` has no missing values in the biomarker columns named by `long_sub_fixed`; there is nothing to impute.", call. = FALSE)
  }

  design_covariates <- imputation_covariates(covariate_vars, id)
  if (length(design_covariates) > 0) {
    cov_formula <- as.formula(paste("~", paste(design_covariates, collapse = " + "), "- 1"))
    cov_design <- model.matrix(cov_formula, data = aligned)
  } else {
    cov_design <- matrix(numeric(0), nrow = nrow(aligned), ncol = 0)
  }

  standardize <- function(mat) {
    if (ncol(mat) == 0) return(list(z = mat, mean = numeric(0), sd = numeric(0)))
    mu <- colMeans(mat)
    sdv <- apply(mat, 2, sd)
    sdv[sdv == 0 | is.na(sdv)] <- 1
    z <- sweep(sweep(mat, 2, mu, "-"), 2, sdv, "/")
    list(z = z, mean = mu, sd = sdv)
  }
  cov_std <- standardize(cov_design)

  resp_mean <- vapply(seq_len(M), function(m) mean(resp_mat[mask[, m], m]), numeric(1))
  resp_sd <- vapply(seq_len(M), function(m) sd(resp_mat[mask[, m], m]), numeric(1))
  resp_sd[resp_sd == 0 | is.na(resp_sd)] <- 1
  resp_z <- sweep(sweep(resp_mat, 2, resp_mean, "-"), 2, resp_sd, "/")
  resp_z[!mask] <- 0

  ### classical single-regression-imputation baseline, for the diagnostics
  ### comparison below -- no new dependency, fit on the same covariates
  baseline_mat <- resp_mat
  for (m in seq_len(M)) {
    obs_idx <- mask[, m]
    if (all(obs_idx) || !any(obs_idx)) next
    if (ncol(cov_design) == 0) {
      baseline_mat[!obs_idx, m] <- mean(resp_mat[obs_idx, m])
    } else {
      fit_df <- data.frame(y = resp_mat[obs_idx, m], cov_design[obs_idx, , drop = FALSE])
      fit <- lm(y ~ ., data = fit_df)
      newdata <- data.frame(cov_design[!obs_idx, , drop = FALSE])
      baseline_mat[!obs_idx, m] <- predict(fit, newdata = newdata)
    }
  }

  if (method == "mixed") {
    if (!is.null(seed)) local_r_seed(seed)
    imputations <- mixed_model_imputations(aligned, data_fit_all, long_sub_fixed, long_sub_random,
                                           id, resp_mat, mask, n_imputations)
  } else {
  ### every check above is cheap and needs no optional dependency; only now,
  ### once the arguments are known to be otherwise well-formed, do we require
  ### torch to actually fit the model
  assert_package_installed("torch", "imputeLongitudinal()")
  if (!torch::torch_is_installed()) {
    stop(paste0("imputeLongitudinal() needs the torch backend (LibTorch), which is not installed; ",
                "run torch::install_torch() once to install it."), call. = FALSE)
  }

  ### the MIWAE's importance resampling draws with R's own RNG
  ### (sample.int()/rnorm()), not torch's, so `seed` has to seed both --
  ### seeding torch alone left MIWAE completions irreproducible. The
  ### caller's R RNG state is restored on exit.
  if (!is.null(seed)) {
    torch::torch_manual_seed(seed)
    local_r_seed(seed)
  }

  ### the network's inputs: the covariates and, with `history`, a summary of
  ### the same subject's other visits (see history_features())
  net_std <- if (history) {
    standardize(cbind(cov_design, history_features(aligned, id, time_variable, resp_mat, mask)))
  } else {
    cov_std
  }
  X_cov <- torch::torch_tensor(net_std$z, dtype = torch::torch_float())
  X_resp <- torch::torch_tensor(resp_z, dtype = torch::torch_float())
  X_mask <- torch::torch_tensor(matrix(as.numeric(mask), nrow = nrow(mask)), dtype = torch::torch_float())

  n_cov <- ncol(net_std$z)
  if (method == "miwae") {
    net <- build_miwae_module(n_cov, M, hidden_units, latent_dim)
    train_miwae(net, X_cov, X_resp, X_mask, epochs = epochs, K = importance_samples)
    imputations <- miwae_impute(net, X_cov, X_resp, X_mask, resp_mean, resp_sd,
                                 K = importance_samples, n_imputations = n_imputations)
  } else {
    ### fixed (not user-exposed) width of the sinusoidal timestep embedding
    ### fed into the denoising network -- exposing it as another tunable
    ### argument would add a knob with little practical effect on a small
    ### tabular MLP, so it is kept as an internal implementation constant
    ### instead, matching how e.g. the beta schedule's start/end are not
    ### user-exposed either.
    time_emb_dim <- 16L
    net <- build_diffusion_module(n_cov, M, hidden_units, time_emb_dim)
    trained <- train_diffusion(net, X_cov, X_resp, X_mask, epochs = epochs,
                                diffusion_steps = diffusion_steps, time_emb_dim = time_emb_dim)
    imputations <- diffusion_impute(net, X_cov, X_resp, X_mask, resp_mean, resp_sd,
                                     schedule = trained$schedule, diffusion_steps = diffusion_steps,
                                     time_emb_dim = time_emb_dim, n_imputations = n_imputations)
  }
  }

  build_output_list <- function(resp_filled_mat) {
    aligned_out <- aligned
    aligned_out[, resp_vars] <- resp_filled_mat
    lapply(seq_len(M), function(m) {
      cols <- unique(c(id, time_variable, resp_vars[m], covariate_vars[[m]]))
      aligned_out[, cols, drop = FALSE]
    })
  }

  filled_mean <- Reduce(`+`, imputations) / n_imputations
  resp_mat_single <- resp_mat
  resp_mat_single[!mask] <- filled_mean[!mask]
  data_fit_all_single <- build_output_list(resp_mat_single)

  data_fit_all_list <- NULL
  if (impute == "multiple") {
    data_fit_all_list <- lapply(seq_len(n_imputations), function(imp) {
      resp_mat_m <- resp_mat
      resp_mat_m[!mask] <- imputations[[imp]][!mask]
      build_output_list(resp_mat_m)
    })
  }

  summarize_marker <- function(observed_vals, imputed_vals) {
    list(n_observed = length(observed_vals), n_imputed = length(imputed_vals),
         observed_mean = mean(observed_vals), observed_sd = sd(observed_vals),
         imputed_mean = if (length(imputed_vals) > 0) mean(imputed_vals) else NA_real_,
         imputed_sd = if (length(imputed_vals) > 1) sd(imputed_vals) else NA_real_)
  }
  diagnostics <- lapply(seq_len(M), function(m) {
    obs_idx <- mask[, m]
    list(biomarker = resp_vars[m], n_missing = sum(!obs_idx),
         deep_generative = summarize_marker(resp_mat[obs_idx, m], resp_mat_single[!obs_idx, m]),
         classical_baseline = summarize_marker(resp_mat[obs_idx, m], baseline_mat[!obs_idx, m]))
  })
  names(diagnostics) <- resp_vars

  list(data_fit_all = data_fit_all_single, data_fit_all_list = data_fit_all_list,
       diagnostics = diagnostics)
}

#' Covariates of the imputation model
#'
#' @description The right-hand-side variables of \code{long_sub_fixed}
#' other than the subject id, including the visit time
#' \code{time_variable} when a formula uses it. The visit time is also the
#' key that aligns rows across biomarkers, and it used to be left out of the
#' model along with the id, so values were imputed without knowing when the
#' measurement was taken (pulling a trending biomarker towards its overall
#' level).
#' @param covariate_vars List of each biomarker's right-hand-side variables.
#' @param id Name of the id column.
#' @return Character vector of covariate names.
#' @keywords internal
imputation_covariates <- function(covariate_vars, id) {
  setdiff(unique(unlist(covariate_vars)), id)
}

#' Build the MIWAE encoder/decoder network
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}. The
#' encoder input is (standardized covariates, zero-filled standardized
#' biomarkers, missingness mask); the decoder input is (latent draw,
#' standardized covariates), matching the "condition the generative model
#' on fully-observed covariates too" design used throughout.
#' @keywords internal
build_miwae_module <- function(n_cov, M, hidden_units, latent_dim) {
  encoder_input_dim <- n_cov + 2 * M
  decoder_input_dim <- latent_dim + n_cov

  miwae_module <- torch::nn_module(
    "miwae_module",
    initialize = function() {
      enc_layers <- list()
      in_dim <- encoder_input_dim
      for (h in hidden_units) {
        enc_layers[[length(enc_layers) + 1]] <- torch::nn_linear(in_dim, h)
        enc_layers[[length(enc_layers) + 1]] <- torch::nn_relu()
        in_dim <- h
      }
      self$encoder_body <- do.call(torch::nn_sequential, enc_layers)
      self$encoder_mean <- torch::nn_linear(in_dim, latent_dim)
      self$encoder_logvar <- torch::nn_linear(in_dim, latent_dim)

      dec_layers <- list()
      in_dim <- decoder_input_dim
      for (h in rev(hidden_units)) {
        dec_layers[[length(dec_layers) + 1]] <- torch::nn_linear(in_dim, h)
        dec_layers[[length(dec_layers) + 1]] <- torch::nn_relu()
        in_dim <- h
      }
      self$decoder_body <- do.call(torch::nn_sequential, dec_layers)
      self$decoder_mean <- torch::nn_linear(in_dim, M)
      self$decoder_logvar <- torch::nn_linear(in_dim, M)
    },
    encode = function(enc_input) {
      h <- self$encoder_body(enc_input)
      list(mean = self$encoder_mean(h), logvar = self$encoder_logvar(h))
    },
    decode = function(z, cov) {
      dec_input <- torch::torch_cat(list(z, cov), dim = -1)
      h <- self$decoder_body(dec_input)
      list(mean = self$decoder_mean(h), logvar = self$decoder_logvar(h))
    }
  )
  miwae_module()
}

#' Draw K importance-weighted latent samples and decode them
#'
#' @description Internal helper shared by \code{\link{train_miwae}} and
#' \code{\link{miwae_impute}}: encodes each row once, draws K samples from
#' q(z|x_o) (reusing the same encoder mean/logvar for all K, as in the
#' original MIWAE), decodes all of them, and returns everything needed to
#' compute the importance-weighted log-likelihood terms.
#' @keywords internal
miwae_forward <- function(net, X_cov, X_resp, X_mask, K) {
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  enc_input <- torch::torch_cat(list(X_cov, X_resp, X_mask), dim = -1)
  q <- net$encode(enc_input)
  q_mean_k <- q$mean$unsqueeze(2)$expand(c(n, K, q$mean$shape[2]))
  q_logvar_k <- q$logvar$unsqueeze(2)$expand(c(n, K, q$logvar$shape[2]))
  eps <- torch::torch_randn_like(q_mean_k)
  z <- q_mean_k + eps * torch::torch_exp(0.5 * q_logvar_k)

  cov_k <- X_cov$unsqueeze(2)$expand(c(n, K, X_cov$shape[2]))
  latent_dim <- z$shape[3]
  dec <- net$decode(z$reshape(c(n * K, latent_dim)), cov_k$reshape(c(n * K, X_cov$shape[2])))
  dec_mean <- dec$mean$reshape(c(n, K, M))
  dec_logvar <- dec$logvar$reshape(c(n, K, M))

  resp_k <- X_resp$unsqueeze(2)$expand(c(n, K, M))
  mask_k <- X_mask$unsqueeze(2)$expand(c(n, K, M))
  log2pi <- log(2 * pi)

  log_px <- -0.5 * (log2pi + dec_logvar + (resp_k - dec_mean)^2 / torch::torch_exp(dec_logvar))
  log_px_obs <- (log_px * mask_k)$sum(dim = 3)
  log_pz <- (-0.5 * (log2pi + z^2))$sum(dim = 3)
  log_qz <- (-0.5 * (log2pi + q_logvar_k + (z - q_mean_k)^2 / torch::torch_exp(q_logvar_k)))$sum(dim = 3)
  log_w <- log_px_obs + log_pz - log_qz

  list(log_w = log_w, dec_mean = dec_mean, dec_logvar = dec_logvar)
}

#' Train the MIWAE network with a full-batch importance-weighted ELBO
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}.
#' @keywords internal
train_miwae <- function(net, X_cov, X_resp, X_mask, epochs, K, learning_rate = 1e-3) {
  optimizer <- torch::optim_adam(net$parameters, lr = learning_rate)
  for (epoch in seq_len(epochs)) {
    optimizer$zero_grad()
    fwd <- miwae_forward(net, X_cov, X_resp, X_mask, K)
    iwae_bound <- torch::torch_logsumexp(fwd$log_w, dim = 2) - log(K)
    loss <- -iwae_bound$mean()
    loss$backward()
    optimizer$step()
  }
  invisible(net)
}

#' Draw self-normalized-importance-resampled completions from a fitted MIWAE
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}:
#' implements the MIWAE multiple-imputation procedure (Mattei & Frellsen,
#' 2019, Section 3.2) -- draw K decoder samples per row, importance-weight
#' them by their (masked-to-observed-entries) likelihood, then resample one
#' index per completion and use that draw's decoder mean plus Gaussian
#' noise as the imputed value. Each completion uses its own K samples, so
#' completions are independent draws. Already-observed cells are left
#' untouched by the caller; this only needs to return sensible values for
#' missing cells.
#' @keywords internal
miwae_impute <- function(net, X_cov, X_resp, X_mask, resp_mean, resp_sd, K, n_imputations) {
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  imputations <- vector("list", n_imputations)

  torch::with_no_grad({
    for (imp in seq_len(n_imputations)) {
      ### a fresh importance-sampling pass for every imputation, so the
      ### completions are independent draws (they used to resample one shared
      ### set of K latent samples)
      fwd <- miwae_forward(net, X_cov, X_resp, X_mask, K)
      w_arr <- torch::as_array(torch::nnf_softmax(fwd$log_w, dim = 2))
      dec_mean_arr <- torch::as_array(fwd$dec_mean)
      dec_logvar_arr <- torch::as_array(fwd$dec_logvar)

      completed_z <- matrix(nrow = n, ncol = M)
      for (i in seq_len(n)) {
        k_draw <- sample.int(K, size = 1, prob = w_arr[i, ])
        mu <- dec_mean_arr[i, k_draw, ]
        sdv <- exp(0.5 * dec_logvar_arr[i, k_draw, ])
        completed_z[i, ] <- mu + rnorm(M) * sdv
      }
      imputations[[imp]] <- sweep(sweep(completed_z, 2, resp_sd, "*"), 2, resp_mean, "+")
    }
  })

  imputations
}

#' Build a linear variance-preserving diffusion (DDPM) noise schedule
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend. Returns the standard Ho et al.
#' (2020) linear \code{beta} schedule and the derived \code{alpha}/
#' \code{alpha_bar} sequences, as plain numeric vectors (not torch tensors)
#' so \code{\link{diffusion_impute}} can index them by a plain R integer
#' timestep inside its reverse-sampling loop without repeated
#' tensor<->R round trips.
#' @keywords internal
diffusion_beta_schedule <- function(diffusion_steps, beta_start = 1e-4, beta_end = 0.02) {
  beta <- seq(beta_start, beta_end, length.out = diffusion_steps)
  alpha <- 1 - beta
  alpha_bar <- cumprod(alpha)
  list(beta = beta, alpha = alpha, alpha_bar = alpha_bar)
}

#' Sinusoidal diffusion-timestep embedding
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend: the standard Transformer/DDPM
#' sinusoidal position embedding, applied to a (possibly per-row-varying)
#' integer timestep instead of a sequence position. The frequency basis is
#' computed once in plain R (not with \code{torch::torch_arange()}, whose
#' inclusive/exclusive endpoint convention is a frequent source of
#' off-by-one bugs) and then multiplied against the timestep tensor, so the
#' output width is always exactly \code{dim} regardless of parity.
#' @param t A 1-D torch float tensor of timesteps, one per row.
#' @param dim Output embedding width.
#' @keywords internal
sinusoidal_time_embedding <- function(t, dim) {
  half <- max(dim %/% 2, 1)
  freqs <- exp(-log(10000) * (seq_len(half) - 1) / max(half - 1, 1))
  freqs_t <- torch::torch_tensor(freqs, dtype = torch::torch_float())
  args <- t$unsqueeze(2) * freqs_t$unsqueeze(1)
  emb <- torch::torch_cat(list(torch::torch_sin(args), torch::torch_cos(args)), dim = -1)
  if (emb$shape[2] < dim) {
    pad <- torch::torch_zeros(c(t$shape[1], dim - emb$shape[2]))
    emb <- torch::torch_cat(list(emb, pad), dim = -1)
  }
  emb
}

#' Build the conditional diffusion denoising network
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend. A single MLP predicts the noise in
#' the target entries of a row from: the noisy target entries (zero
#' elsewhere), the conditioning entries' values (zero elsewhere), the
#' conditioning mask, the covariates and an embedding of the diffusion step
#' -- a tabular version of CSDI's conditional denoiser (Tashiro et al.,
#' 2021).
#' @param n_cov Number of (standardized) covariate columns.
#' @param M Number of biomarkers.
#' @param hidden_units Hidden-layer sizes of the MLP body.
#' @param time_emb_dim Width of the timestep embedding.
#' @return An instantiated \code{torch::nn_module}.
#' @keywords internal
build_diffusion_module <- function(n_cov, M, hidden_units, time_emb_dim) {
  input_dim <- 3 * M + n_cov + time_emb_dim

  diffusion_module <- torch::nn_module(
    "diffusion_module",
    initialize = function() {
      self$time_mlp <- torch::nn_sequential(
        torch::nn_linear(time_emb_dim, time_emb_dim), torch::nn_relu()
      )
      body_layers <- list()
      in_dim <- input_dim
      for (h in hidden_units) {
        body_layers[[length(body_layers) + 1]] <- torch::nn_linear(in_dim, h)
        body_layers[[length(body_layers) + 1]] <- torch::nn_relu()
        in_dim <- h
      }
      self$body <- do.call(torch::nn_sequential, body_layers)
      self$out <- torch::nn_linear(in_dim, M)
    },
    forward = function(x_t, x_cond, cond_mask, cov, t_emb) {
      te <- self$time_mlp(t_emb)
      h <- self$body(torch::torch_cat(list(x_t, x_cond, cond_mask, cov, te), dim = -1))
      self$out(h)
    }
  )
  diffusion_module()
}

#' Train the conditional diffusion network
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend, following CSDI's self-supervised
#' training (Tashiro et al., 2021): in every mini-batch, each row's observed
#' entries are split at random into targets (each with a probability drawn
#' uniformly per row) and conditions; the targets are noised to a random
#' diffusion step and the network learns to predict that noise from the
#' conditions' true values. Only observed entries ever serve as targets, so
#' the loss never touches a missing cell. An earlier version noised whole
#' rows without giving the network the other observed values, so imputations
#' made little use of the row's other biomarkers.
#' @param net Network from \code{build_diffusion_module()}.
#' @param X_cov,X_resp,X_mask Covariates, standardized responses (0 where
#'   missing) and the observation mask, as \code{torch} tensors.
#' @param epochs Passes over the data.
#' @param diffusion_steps Number of diffusion steps.
#' @param time_emb_dim Width of the timestep embedding.
#' @param batch_size Rows per mini-batch.
#' @param learning_rate Adam learning rate.
#' @return A list with the trained \code{net} and its noise \code{schedule}.
#' @keywords internal
train_diffusion <- function(net, X_cov, X_resp, X_mask, epochs, diffusion_steps,
                            time_emb_dim, batch_size = 128, learning_rate = 1e-3) {
  schedule <- diffusion_beta_schedule(diffusion_steps)
  alpha_bar <- torch::torch_tensor(schedule$alpha_bar, dtype = torch::torch_float())
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  optimizer <- torch::optim_adam(net$parameters, lr = learning_rate)

  for (epoch in seq_len(epochs)) {
    order <- sample.int(n)
    for (start in seq(1, n, by = batch_size)) {
      idx <- torch::torch_tensor(order[start:min(n, start + batch_size - 1)], dtype = torch::torch_long())
      x0 <- X_resp[idx, ]
      observed <- X_mask[idx, ]
      b <- x0$shape[1]
      share <- torch::torch_rand(c(b, 1))
      target <- observed * (torch::torch_rand(c(b, M)) < share)$to(dtype = torch::torch_float())
      cond <- observed - target
      t_idx <- torch::torch_randint(1, diffusion_steps + 1, c(b), dtype = torch::torch_long())
      ab_t <- alpha_bar[t_idx]$unsqueeze(2)
      eps <- torch::torch_randn_like(x0)
      x_t <- (torch::torch_sqrt(ab_t) * x0 + torch::torch_sqrt(1 - ab_t) * eps) * target
      t_emb <- sinusoidal_time_embedding(t_idx$to(dtype = torch::torch_float()), time_emb_dim)
      eps_hat <- net(x_t, x0 * cond, cond, X_cov[idx, ], t_emb)
      loss <- ((eps_hat - eps)^2 * target)$sum() / torch::torch_clamp(target$sum(), min = 1)
      optimizer$zero_grad()
      loss$backward()
      optimizer$step()
    }
  }

  list(net = net, schedule = schedule)
}

#' Draw conditional reverse-diffusion completions
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend. The missing entries of every row are
#' the targets and all its observed entries the conditions: starting from
#' pure noise on the targets, the DDPM reverse process (Ho et al., 2020) is
#' run with the network conditioned on the observed values at every step.
#' Observed entries are returned unchanged.
#' @param net Trained network.
#' @param X_cov,X_resp,X_mask As in \code{train_diffusion()}.
#' @param resp_mean,resp_sd Means and standard deviations that undo the
#'   responses' standardization.
#' @param schedule,diffusion_steps,time_emb_dim As used in training.
#' @param n_imputations Number of completions to draw.
#' @return A list of \code{n_imputations} completed response matrices on the
#'   original scale.
#' @keywords internal
diffusion_impute <- function(net, X_cov, X_resp, X_mask, resp_mean, resp_sd,
                             schedule, diffusion_steps, time_emb_dim, n_imputations) {
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  target <- 1 - X_mask
  x_cond <- X_resp * X_mask
  imputations <- vector("list", n_imputations)

  torch::with_no_grad({
    for (imp in seq_len(n_imputations)) {
      x <- torch::torch_randn(c(n, M)) * target
      for (t in diffusion_steps:1) {
        t_emb <- sinusoidal_time_embedding(
          torch::torch_tensor(rep(t, n), dtype = torch::torch_float()), time_emb_dim)
        eps_hat <- net(x, x_cond, X_mask, X_cov, t_emb)
        alpha_bar_t <- schedule$alpha_bar[t]
        beta_t <- schedule$beta[t]
        mu <- (x - (beta_t / sqrt(1 - alpha_bar_t)) * eps_hat) / sqrt(schedule$alpha[t])
        if (t > 1) {
          # variance of the true reverse step q(x_{t-1} | x_t, x_0)
          var_t <- beta_t * (1 - schedule$alpha_bar[t - 1]) / (1 - alpha_bar_t)
          x <- (mu + sqrt(var_t) * torch::torch_randn_like(x)) * target
        } else {
          x <- mu * target
        }
      }
      x_arr <- torch::as_array(X_mask * X_resp + target * x)
      imputations[[imp]] <- sweep(sweep(x_arr, 2, resp_sd, "*"), 2, resp_mean, "+")
    }
  })

  imputations
}

#' Seed R's RNG until the calling function returns
#'
#' @description Calls \code{set.seed(seed)} and registers, in the calling
#' function's frame, an exit handler that puts back the
#' \code{.Random.seed} that existed before (or removes it if there was
#' none), so a \code{seed} argument does not change the caller's random
#' number stream.
#'
#' @param seed Integer seed.
#' @param frame The frame whose exit restores the RNG state.
#' @keywords internal
local_r_seed <- function(seed, frame = parent.frame()) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  restore <- function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
  do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = frame)
  set.seed(seed)
}

#' What a subject's other visits say about each biomarker
#'
#' @description Inputs for the deep imputation models (\code{history = TRUE}
#' in \code{\link{imputeLongitudinal}}): without them each visit is
#' completed from that visit alone, so an imputed value does not carry the
#' subject's own level, which shrinks the estimated variation between
#' subjects. For each visit and biomarker: the mean of the subject's
#' \emph{other} observed values of the biomarker, how many there are
#' (\code{log1p}), and the value at the nearest other observed visit and its
#' distance in time (\code{log1p}). The visit's own value is left out, so
#' the model cannot simply copy it while training. A subject with no other
#' observed value gets the biomarker's overall mean, a count of 0 and the
#' largest distance.
#' @return A matrix with four columns per biomarker, one row per visit.
#' @keywords internal
history_features <- function(aligned, id, time_variable, resp_mat, mask) {
  n <- nrow(aligned)
  M <- ncol(resp_mat)
  times <- aligned[[time_variable]]
  overall <- vapply(seq_len(M), function(m) mean(resp_mat[mask[, m], m]), numeric(1))
  far <- log1p(diff(range(times)) + 1)
  out <- matrix(NA_real_, n, 4 * M)
  colnames(out) <- paste0(rep(colnames(resp_mat), each = 4),
                          c(".other_mean", ".other_count", ".nearest", ".nearest_distance"))
  for (rows in split(seq_len(n), as.character(aligned[[id]]))) {
    t <- times[rows]
    for (m in seq_len(M)) {
      observed <- mask[rows, m]
      y <- resp_mat[rows, m]
      for (j in seq_along(rows)) {
        others <- observed
        others[j] <- FALSE
        cols <- 4 * (m - 1) + 1:4
        if (!any(others)) {
          out[rows[j], cols] <- c(overall[m], 0, overall[m], far)
        } else {
          k <- which(others)[which.min(abs(t[others] - t[j]))]
          out[rows[j], cols] <- c(mean(y[others]), log1p(sum(others)), y[k], log1p(abs(t[k] - t[j])))
        }
      }
    }
  }
  out
}

#' Multiple imputation from the multivariate linear mixed model
#'
#' @description \code{method = "mixed"} in \code{\link{imputeLongitudinal}}.
#' Fits the model \code{\link{longitudinalSub}} fits -- each biomarker's
#' linear mixed model, with the random effects of all biomarkers jointly
#' normal with covariance \code{Sigma_fit} -- and draws each subject's
#' missing values from their normal distribution given everything observed
#' for that subject: all of its visits and all biomarkers. For each
#' completion the model's parameters are first drawn from an approximation
#' to their posterior, so that the completions also reflect the uncertainty
#' in the fit (proper multiple imputation): the fixed effects from their
#' estimated sampling distribution, each residual variance from a scaled
#' inverse chi-squared distribution, and \code{Sigma_fit} from an inverse
#' Wishart distribution with as many degrees of freedom as subjects.
#' @return A list of \code{n_imputations} matrices like \code{resp_mat},
#'   with the missing cells filled.
#' @keywords internal
mixed_model_imputations <- function(aligned, data_fit_all, long_sub_fixed, long_sub_random, id,
                                    resp_mat, mask, n_imputations) {
  M <- length(long_sub_fixed)
  fit <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)
  lfit <- fit$lfit
  X <- lapply(seq_len(M), function(m) {
    xlev <- if (!is.null(fit$xlevels)) fit$xlevels[[m]] else NULL
    mf <- model.frame(lfit[[m]]$terms, aligned, xlev = xlev, na.action = stats::na.pass)
    model.matrix(lfit[[m]]$terms, mf, contrasts.arg = lfit[[m]]$contrasts)
  })
  ### subjects the joint covariance was estimated from: some value of every biomarker
  ids <- as.character(aligned[[id]])
  n_subjects <- sum(tapply(seq_len(nrow(aligned)), ids, function(r) all(colSums(mask[r, , drop = FALSE]) > 0)))
  Sigma_hat <- fit$Sigma_fit
  rows_by_id <- split(seq_len(nrow(aligned)), ids)
  Z_by_id <- lapply(rows_by_id, function(rows) {
    random_effects_design(rep(list(aligned[rows, , drop = FALSE]), M), long_sub_random, Sigma_hat)
  })

  lapply(seq_len(n_imputations), function(k) {
    beta <- lapply(lfit, function(f) {
      as.vector(mvtnorm::rmvnorm(1, nlme::fixef(f), as.matrix(f$varFix)))
    })
    sigma2 <- vapply(seq_len(M), function(m) {
      df <- lfit[[m]]$dims$N - length(beta[[m]])
      lfit[[m]]$sigma^2 * df / stats::rchisq(1, df)
    }, numeric(1))
    Sigma <- solve(stats::rWishart(1, max(n_subjects, nrow(Sigma_hat) + 1),
                                   solve(Sigma_hat) / max(n_subjects, nrow(Sigma_hat) + 1))[, , 1])
    filled <- resp_mat
    for (s in names(rows_by_id)) {
      rows <- rows_by_id[[s]]
      obs <- as.vector(mask[rows, , drop = FALSE])           # stacked by biomarker, then visit
      if (all(obs)) next
      mu <- unlist(lapply(seq_len(M), function(m) X[[m]][rows, , drop = FALSE] %*% beta[[m]]))
      Z <- Z_by_id[[s]]
      V <- Z %*% Sigma %*% t(Z) + diag(rep(sigma2, each = length(rows)), length(mu))
      y <- as.vector(resp_mat[rows, , drop = FALSE])
      mis <- !obs
      if (any(obs)) {
        gain <- V[mis, obs, drop = FALSE] %*% solve(V[obs, obs, drop = FALSE])
        cond_mean <- mu[mis] + as.vector(gain %*% (y[obs] - mu[obs]))
        cond_cov <- V[mis, mis, drop = FALSE] - gain %*% V[obs, mis, drop = FALSE]
      } else {
        cond_mean <- mu[mis]
        cond_cov <- V[mis, mis, drop = FALSE]
      }
      cond_cov <- (cond_cov + t(cond_cov)) / 2
      draw <- cond_mean + as.vector(t(chol(cond_cov + diag(1e-10, nrow(cond_cov)))) %*% stats::rnorm(sum(mis)))
      y[mis] <- draw
      filled[rows, ] <- matrix(y, length(rows), M)
    }
    filled
  })
}
