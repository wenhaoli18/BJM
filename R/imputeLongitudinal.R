#' Impute missing longitudinal biomarker values with a deep generative model
#'
#' @description \code{longitudinalSub()} fits each biomarker's mixed model on
#' the complete cases for that biomarker, and then keeps only the subjects
#' who have at least one non-missing observation of \emph{every} biomarker
#' (see \code{\link{longitudinalSub}}). With interrupted/irregular follow-up this
#' throws away information: a subject missing just one of several biomarkers
#' at a visit -- or missing a biomarker's measurements entirely -- is dropped
#' from every biomarker's fit, not just the one it is missing.
#'
#' \code{imputeLongitudinal()} is an optional preprocessing step that fills
#' these gaps \emph{before} \code{longitudinalSub()} runs, assuming the
#' missingness is at random (MAR) given the observed covariates and
#' biomarkers. Two deep generative backends are available, selected with
#' \code{method}:
#' \itemize{
#'   \item \code{"miwae"} (default) fits a MIWAE (Missing data
#'   Importance-Weighted AutoEncoder; Mattei & Frellsen, 2019) -- a
#'   variational autoencoder whose training objective only scores
#'   reconstruction on the observed entries -- jointly across all supplied
#'   biomarkers, and uses self-normalized importance resampling from the
#'   fitted decoder to draw one or more plausible completions of the
#'   missing cells.
#'   \item \code{"diffusion"} fits a conditional denoising diffusion
#'   probabilistic model (DDPM; Ho, Jain & Abbeel, 2020) over the same
#'   per-row biomarker vector, trained with a masked denoising
#'   score-matching loss that -- like the MIWAE objective above -- is only
#'   evaluated on observed entries (the same masking principle used by
#'   Ouyang et al.'s MissDiff, 2023, for training diffusion models on
#'   incomplete tabular data). Missing cells are then filled by a
#'   RePaint-style (Lugmayr et al., 2022) reverse-diffusion sampler: at
#'   every denoising step, observed dimensions are re-noised directly from
#'   their true value and only the unobserved dimensions are drawn from the
#'   learned reverse process, so the final sample matches the data exactly
#'   wherever it was actually observed.
#' }
#' Both backends draw one or more plausible completions of the missing
#' cells, and the completed data.frame(s) can be passed directly to
#' \code{longitudinalSub()} in place of the original \code{data_fit_all}.
#'
#' @details Only the biomarker (response) columns named on the left-hand
#' side of \code{long_sub_fixed} are imputed; every covariate referenced on
#' the right-hand side of \code{long_sub_fixed}/\code{long_sub_random} is
#' assumed fully observed (an error is raised if any covariate has missing
#' values). Each biomarker's own within-subject serial correlation continues
#' to be modeled downstream by \code{nlme::lme()} inside
#' \code{longitudinalSub()}, exactly as today; this function only models
#' \code{p(biomarkers | covariates, missingness mask)} row by row (one row
#' per subject-visit), so it does not need to know the random-effects
#' structure in \code{long_sub_random} beyond the id variable it names.
#'
#' Rows are aligned across biomarkers by \code{(id, time_variable)}: when
#' \code{data_fit_all} is a list of different data.frames with different
#' visit schedules per biomarker (rather than one data.frame shared by
#' every biomarker), the union of every \code{(id, time_variable)} pair seen
#' in any of them is used, so a visit at which one biomarker was measured
#' but another was not still gets a value imputed for the missing one.
#'
#' Because this is a deep generative model, it needs enough data to fit
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
#' @param epochs Number of full-batch training epochs.
#' @param importance_samples Number of importance samples (\code{K}) used
#' both in the training objective and at imputation time. Only used when
#' \code{method = "miwae"}.
#' @param diffusion_steps Number of forward/reverse diffusion timesteps.
#' Only used when \code{method = "diffusion"}; more steps generally
#' improve sample quality at proportionally higher training/sampling cost.
#' @param method Which deep generative backend to fit: \code{"miwae"}
#' (default) or \code{"diffusion"}. See Details.
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
#' imputed cells. Note also that each row (subject-visit) is imputed from
#' that row's covariates and observed biomarkers only, not from the same
#' subject's other visits, so imputed values do not carry a subject's own
#' level or trend; this too pulls the estimated between-subject
#' (random-effects) variation towards zero.
#' @param seed Optional integer seed for reproducibility. It seeds both
#' \pkg{torch} and R's random number generator for the duration of the
#' call; R's random number state from before the call is restored
#' afterwards.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{data_fit_all}{A completed version of \code{data_fit_all}, in the
#'   same list/data.frame shape, with every missing biomarker cell filled.
#'   Ready to pass straight to \code{\link{longitudinalSub}}.}
#'   \item{data_fit_all_list}{When \code{impute = "multiple"}, a list of
#'   \code{n_imputations} independently completed versions of
#'   \code{data_fit_all}; \code{NULL} when \code{impute = "single"}.}
#'   \item{diagnostics}{A named list, one element per biomarker, each with
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
#' @references Lugmayr, A., Danelljan, M., Romero, A., Yu, F., Timofte, R.,
#' and Van Gool, L. (2022). RePaint: Inpainting using Denoising Diffusion
#' Probabilistic Models. \emph{Proceedings of the IEEE/CVF Conference on
#' Computer Vision and Pattern Recognition}, 11461-11471.
#' @references Ouyang, Y., Xie, L., Li, C., and Cheng, G. (2023). MissDiff:
#' Training Diffusion Models on Tabular Data with Missing Values.
#' \emph{ICML 2023 Workshop on Structured Probabilistic Inference &
#' Generative Modeling}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE)) {
#'   data(pbc3)
#'   data_fit_all <- pbc3[pbc3$status3 == 1, ]
#'
#'   # pbc3's complete-case subset has no missing serBilir/albumin values;
#'   # artificially delete some at random (MCAR, a special case of MAR) to
#'   # simulate the interrupted-follow-up gaps imputeLongitudinal() targets
#'   set.seed(1)
#'   n <- nrow(data_fit_all)
#'   data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
#'   data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA
#'
#'   long_sub_fixed <- list(
#'     "long1" = serBilir ~ year + age + sex + years,
#'     "long2" = albumin ~ year + age + sex + years)
#'   long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)
#'
#'   imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed,
#'                                  long_sub_random, time_variable = "year",
#'                                  epochs = 50, seed = 1)
#'   long_fit_all <- longitudinalSub(imputed$data_fit_all, long_sub_fixed,
#'                                    long_sub_random)
#'
#'   # the diffusion backend fills the same gaps via a conditional DDPM
#'   # instead of a VAE -- same call shape, just a different `method`
#'   imputed_diffusion <- imputeLongitudinal(data_fit_all, long_sub_fixed,
#'                                            long_sub_random, time_variable = "year",
#'                                            method = "diffusion", epochs = 50,
#'                                            diffusion_steps = 30, seed = 1)
#' }
#' }
#'
#' @export
imputeLongitudinal <- function(data_fit_all, long_sub_fixed, long_sub_random,
                                time_variable, n_imputations = 5,
                                latent_dim = 8, hidden_units = c(64, 32),
                                epochs = 300, importance_samples = 20,
                                diffusion_steps = 100,
                                method = c("miwae", "diffusion"),
                                impute = c("single", "multiple"),
                                seed = NULL) {
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
  }

  resp_vars <- unname(vapply(long_sub_fixed, function(f) as.character(f[[2]]), character(1)))
  if (any(duplicated(resp_vars))) {
    stop(sprintf(
      "`long_sub_fixed` response variables must be distinct across biomarkers; found a repeated response: %s.",
      resp_vars[duplicated(resp_vars)][1]
    ), call. = FALSE)
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

  if (length(all_covariates) > 0) {
    cov_formula <- as.formula(paste("~", paste(all_covariates, collapse = " + "), "- 1"))
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

  ### every check above is cheap and needs no optional dependency; only now,
  ### once the arguments are known to be otherwise well-formed, do we require
  ### torch to actually fit the model
  assert_package_installed("torch", "imputeLongitudinal()")

  ### the MIWAE's importance resampling draws with R's own RNG
  ### (sample.int()/rnorm()), not torch's, so `seed` has to seed both --
  ### seeding torch alone left MIWAE completions irreproducible. The
  ### caller's R RNG state is restored on exit.
  if (!is.null(seed)) {
    torch::torch_manual_seed(seed)
    local_r_seed(seed)
  }

  X_cov <- torch::torch_tensor(cov_std$z, dtype = torch::torch_float())
  X_resp <- torch::torch_tensor(resp_z, dtype = torch::torch_float())
  X_mask <- torch::torch_tensor(matrix(as.numeric(mask), nrow = nrow(mask)), dtype = torch::torch_float())

  n_cov <- ncol(cov_std$z)
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
#' noise as the imputed value. Already-observed cells are left untouched by
#' the caller; this only needs to return sensible values for missing cells.
#' @keywords internal
miwae_impute <- function(net, X_cov, X_resp, X_mask, resp_mean, resp_sd, K, n_imputations) {
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  imputations <- vector("list", n_imputations)

  torch::with_no_grad({
    fwd <- miwae_forward(net, X_cov, X_resp, X_mask, K)
    w <- torch::nnf_softmax(fwd$log_w, dim = 2)

    dec_mean_arr <- torch::as_array(fwd$dec_mean)
    dec_logvar_arr <- torch::as_array(fwd$dec_logvar)
    w_arr <- torch::as_array(w)

    for (imp in seq_len(n_imputations)) {
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

#' Build the diffusion denoising network
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend. A single MLP predicts the noise
#' \code{eps} added to the (standardized) biomarker vector, conditioned on
#' the noisy input itself, the observed-data mask, the fully-observed
#' covariates, and an embedding of the current timestep -- the same
#' "condition the generative model on fully-observed covariates and the
#' mask" design used by \code{\link{build_miwae_module}}.
#' @keywords internal
build_diffusion_module <- function(n_cov, M, hidden_units, time_emb_dim) {
  input_dim <- M + M + n_cov + time_emb_dim

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
    forward = function(x_t, mask, cov, t_emb) {
      te <- self$time_mlp(t_emb)
      input <- torch::torch_cat(list(x_t, mask, cov, te), dim = -1)
      h <- self$body(input)
      self$out(h)
    }
  )
  diffusion_module()
}

#' Train the diffusion network with a masked denoising score-matching loss
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend. At every epoch, each row is noised
#' to an independently-sampled random timestep, the network predicts the
#' noise that was added, and the squared-error loss is averaged only over
#' \emph{observed} entries (\code{X_mask == 1}) -- missing entries have no
#' real \code{x0} to noise in the first place, so scoring them would just
#' train the network to fit whatever placeholder (zero) value they were
#' filled with. This mirrors \code{\link{train_miwae}}'s
#' observed-entries-only reconstruction term, and matches the masked
#' training objective used by MissDiff (Ouyang et al., 2023) for diffusion
#' models on incomplete tabular data.
#' @keywords internal
train_diffusion <- function(net, X_cov, X_resp, X_mask, epochs, diffusion_steps,
                             time_emb_dim, learning_rate = 1e-3) {
  schedule <- diffusion_beta_schedule(diffusion_steps)
  n <- X_resp$shape[1]
  optimizer <- torch::optim_adam(net$parameters, lr = learning_rate)

  for (epoch in seq_len(epochs)) {
    optimizer$zero_grad()

    t_idx <- sample.int(diffusion_steps, size = n, replace = TRUE)
    ab_t <- torch::torch_tensor(schedule$alpha_bar[t_idx], dtype = torch::torch_float())$unsqueeze(2)
    eps <- torch::torch_randn_like(X_resp)
    x_t <- torch::torch_sqrt(ab_t) * X_resp + torch::torch_sqrt(1 - ab_t) * eps

    t_emb <- sinusoidal_time_embedding(
      torch::torch_tensor(t_idx, dtype = torch::torch_float()), time_emb_dim)
    eps_hat <- net(x_t, X_mask, X_cov, t_emb)

    sq_err <- (eps_hat - eps)^2 * X_mask
    loss <- sq_err$sum() / torch::torch_clamp(X_mask$sum(), min = 1)
    loss$backward()
    optimizer$step()
  }

  list(net = net, schedule = schedule)
}

#' Draw RePaint-style conditional reverse-diffusion completions
#'
#' @description Internal helper for \code{\link{imputeLongitudinal}}'s
#' \code{method = "diffusion"} backend: implements the RePaint (Lugmayr et
#' al., 2022) inpainting sampler, adapted from image patches to a per-row
#' biomarker vector. Starting from pure noise, at every reverse step
#' \code{t}, observed dimensions are overwritten by re-noising the true
#' (standardized) observed value directly to noise level \code{t - 1} (the
#' "known" branch, sampled from the forward process posterior, not the
#' network), while unobserved dimensions are drawn from the network's
#' learned reverse step (the "unknown" branch); the two are recombined with
#' the mask before the next step. This guarantees the final sample matches
#' the data exactly at every previously-observed cell, and only genuinely
#' extrapolates the missing ones. Already-observed cells are left untouched
#' by the caller, matching \code{\link{miwae_impute}}'s contract.
#' @keywords internal
diffusion_impute <- function(net, X_cov, X_resp, X_mask, resp_mean, resp_sd,
                              schedule, diffusion_steps, time_emb_dim, n_imputations) {
  n <- X_resp$shape[1]
  M <- X_resp$shape[2]
  imputations <- vector("list", n_imputations)

  torch::with_no_grad({
    for (imp in seq_len(n_imputations)) {
      x <- torch::torch_randn(c(n, M))

      for (t in diffusion_steps:1) {
        t_emb <- sinusoidal_time_embedding(
          torch::torch_tensor(rep(t, n), dtype = torch::torch_float()), time_emb_dim)
        eps_hat <- net(x, X_mask, X_cov, t_emb)

        alpha_t <- schedule$alpha[t]
        alpha_bar_t <- schedule$alpha_bar[t]
        beta_t <- schedule$beta[t]
        mu <- (x - (beta_t / sqrt(1 - alpha_bar_t)) * eps_hat) / sqrt(alpha_t)

        if (t > 1) {
          z <- torch::torch_randn_like(x)
          x_unknown <- mu + sqrt(beta_t) * z

          alpha_bar_prev <- schedule$alpha_bar[t - 1]
          eps_known <- torch::torch_randn_like(x)
          x_known <- sqrt(alpha_bar_prev) * X_resp + sqrt(1 - alpha_bar_prev) * eps_known
        } else {
          x_unknown <- mu
          x_known <- X_resp
        }

        x <- X_mask * x_known + (1 - X_mask) * x_unknown
      }

      ### numerical hygiene: guarantee observed cells are exactly preserved
      ### regardless of any residual sampling noise accumulated above
      x <- X_mask * X_resp + (1 - X_mask) * x
      x_arr <- torch::as_array(x)
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
