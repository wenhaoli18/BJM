#' Plot of risk using dynamic prediction
#' 
#' @description This function gives the risk prediction plot.
#' 
#' @param data_predict_all_pre This involves a collection of \code{data.frame} objects for
#' dynamic prediction, each corresponding to a distinct longitudinal outcome. These data
#' frames should contain the variables specified in \code{long_sub_fixed} and
#' \code{long_sub_random}. Utilizing a list structure allows for the incorporation of
#' multiple longitudinal outcomes, each potentially following different measurement
#' protocols. In instances where all longitudinal outcomes are recorded at identical time
#' points across patients, a singular \code{data.frame} object may be used in a \code{list}.
#' Alternatively, a single bare \code{data.frame} (not wrapped in a list) may be supplied
#' directly; it is then reused for every longitudinal outcome. It is presumed that each
#' data frame is structured in a long format.
#' @param long_fit_all Outputs from the model fitting process using the \code{nlme} package,
#' encompassing the results and parameters obtained from the analysis.
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param prediction_time Time used to make the prediction.
#' @param bio_i Biomarker used to do prediction. 
#' @param horizon Prediction horizon.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes, 
#' in the order of \code{survival_variable_all}.
#' @param bandcount1 The number of grid points spanning the prediction window,
#' from \code{prediction_time} to \code{prediction_time + horizon}. Larger
#' values give a more accurate but slower estimate. Defaults to \code{"auto"},
#' which resolves it once, before looping over the landmark times, (using
#' the first landmark time as a representative probe) by doubling from a
#' built-in starting value until the predicted risk stabilizes; see
#' \code{\link{dynamicPrediction}}'s \code{bandcount1} for details of that
#' search. The resolved value is then reused, fixed, for every landmark
#' time -- it is not re-searched on every iteration.
#' @param bandcount2 The number of grid points used to approximate
#' integrating out to infinity when normalizing the predicted risk. A wider
#' follow-up range needs a larger \code{bandcount2} to keep the grid
#' spacing comparable. Defaults to \code{"auto"}; resolved the same way as
#' \code{bandcount1} (jointly with it, when both are \code{"auto"}).
#'
#' Pass explicit numbers instead of \code{"auto"} for full manual control, or
#' use \code{checkBandcountConvergence()} (applied to \code{dynamicPrediction()}
#' directly) to inspect the convergence behavior yourself. See also
#' \code{vignette("BJM-intro", package = "BJM")} for further guidance on
#' choosing \code{bandcount1}/\code{bandcount2}.
#'
#' @param n_cores Number of CPU cores to use for computing the prediction at
#' each landmark time in \code{prediction_time}. Each landmark time is
#' computed independently, so this loop can be dispatched across cores.
#' Defaults to \code{1} (serial execution; identical behavior/output to
#' versions of this function without this argument). Values greater than
#' \code{1} use \code{parallel::mclapply()}, which relies on forking and is
#' therefore only actually parallel on Unix-like systems (Linux, macOS); on
#' Windows, \code{mclapply()} silently runs the iterations serially
#' regardless of \code{n_cores} (a limitation of R's fork-based parallelism,
#' not of this package). Parallel execution produces exactly the same
#' numeric result as serial execution -- only the order in which iterations
#' are computed (not the order results are assembled in) changes.
#' @return Plot of risk using dynamic prediction.
#' @export
riskPlot = function(data_predict_all_pre, long_fit_all, survival_fit_all,
                       prediction_time = NULL, bio_i = NULL,
                       horizon, time_variable,
                       survival_variable_all, survival_trans_function,
                       bandcount1 = "auto", bandcount2 = "auto", n_cores = 1){

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_data_list(data_predict_all_pre, "data_predict_all_pre", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all_pre) || is.data.frame(data_predict_all_pre)) {
    data_predict_all_pre <- rep(list(data_predict_all_pre), each = length(long_fit_all$lfit))
  }
  if (!is.null(bio_i)) {
    assert_index(bio_i, length(long_fit_all$lfit), "bio_i", "longitudinal outcomes in long_fit_all")
  }
  assert_string(time_variable, "time_variable")
  assert_scalar_numeric(horizon, "horizon")
  assert_bandcount(bandcount1, "bandcount1")
  assert_bandcount(bandcount2, "bandcount2")
  assert_positive_integer(n_cores, "n_cores")
  # prediction_time may be NULL (landmark defaults to each patient's first
  # observed time_variable value) or a vector of landmark times; probe with
  # the first usable value, or skip the probe entirely if none is available yet.
  survival_trans_probe <- if (is.numeric(prediction_time) && length(prediction_time) >= 1 && !anyNA(prediction_time[1])) {
    prediction_time[1]
  } else {
    NULL
  }
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = survival_trans_probe)
  for (i in seq_along(data_predict_all_pre)) {
    assert_vars_in_data(time_variable, data_predict_all_pre[[i]], "time_variable",
                         sprintf("data_predict_all_pre[[%d]]", i))
  }

  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]
  
  ### event type variable name
  if(length(survival_fit_all$form_conditional_cr) != 0){
    event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  }
  
  #name of biomarker
  if(is.null(bio_i)){
    bio_i_name = as.character(formula(long_fit_all$long_sub_fixed[[1]])[[2]]) 
    
    DP_data_bio = data.frame(time = unlist(data_predict_all_pre[[1]][time_variable]),
                             longitudinal = unlist(data_predict_all_pre[[1]][bio_i_name]))
  }else{
    bio_i_name = as.character(formula(long_fit_all$long_sub_fixed[[bio_i]])[[2]])

    DP_data_bio = data.frame(time = unlist(data_predict_all_pre[[bio_i]][time_variable]),
                             longitudinal = unlist(data_predict_all_pre[[bio_i]][bio_i_name]))
    
  }
  
  
  if(is.null(prediction_time)){
    landmark.time = unlist(getFirst(data_predict_all_pre)[time_variable])
  }else if(length(prediction_time) == 1){
    landmark.time = c(prediction_time, 1.5 * prediction_time, 2 * prediction_time)
  }else{
    landmark.time = prediction_time
  }

  ### bandcount1/bandcount2 = "auto" (the default): resolve them once here,
  ### using the first landmark time as a representative probe, rather than
  ### re-running the auto-tuning search on every landmark time in the loop
  ### below (see auto_tune_bandcount()).
  auto_names_1_2 <- c("bandcount1", "bandcount2")[c(identical(bandcount1, "auto"), identical(bandcount2, "auto"))]
  if (length(auto_names_1_2) > 0) {
    probe_data_predict_all <- list()
    for (i in seq_len(length(long_fit_all$long_sub_fixed))) {
      probe_data_predict_all[[i]] = data_predict_all_pre[[i]][data_predict_all_pre[[i]][time_variable] <= landmark.time[1], ]
    }
    probe_args <- list(data_predict_all = probe_data_predict_all, long_fit_all = long_fit_all,
                        survival_fit_all = survival_fit_all, prediction_time = landmark.time[1],
                        horizon = horizon, time_variable = time_variable,
                        survival_variable_all = survival_variable_all,
                        survival_trans_function = survival_trans_function,
                        bandcount1 = bandcount1, bandcount2 = bandcount2)
    resolved_1_2 <- auto_tune_bandcount(dynamicPrediction, probe_args, auto_names_1_2)$bandcount
    if (!is.null(resolved_1_2$bandcount1)) bandcount1 <- resolved_1_2$bandcount1
    if (!is.null(resolved_1_2$bandcount2)) bandcount2 <- resolved_1_2$bandcount2
  }

  ### each landmark time only reads data_predict_all_pre/long_fit_all/
  ### survival_fit_all/bandcount1/bandcount2 (all fixed above), so it is an
  ### independent unit of work. It is dispatched below via lapply() (serial,
  ### n_cores == 1, the default -- identical to the old sequential for()
  ### loop) or parallel::mclapply() (n_cores > 1), and the per-landmark-time
  ### results are reassembled afterwards, in the original landmark.time
  ### order, into exactly the vectors the old loop built with c(...)
  ### accumulation.
  has_cr <- length(survival_fit_all$form_conditional_cr) != 0

  compute_one_landmark <- function(time.cutoff) {
    data_predict_all = list()
    for(i in seq_len(length(long_fit_all$long_sub_fixed))){
      data_predict_all[[i]] = data_predict_all_pre[[i]][data_predict_all_pre[[i]][time_variable] <= time.cutoff,]
    }

    risk.prob = dynamicPrediction(data_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time = time.cutoff,
                                  horizon, time_variable,
                                  survival_variable_all, survival_trans_function,
                                  bandcount1, bandcount2)

    list(time.cutoff = time.cutoff,
         risk_prob_1 = risk.prob$risk_prob_1,
         risk_prob_2 = if (has_cr) risk.prob$risk_prob_2 else NULL,
         keep = length(risk.prob$risk_prob_1) != 0)
  }

  if (n_cores > 1) {
    results_per_landmark <- parallel::mclapply(landmark.time, compute_one_landmark, mc.cores = n_cores)
    # Unlike lapply(), mclapply() does not propagate an error raised inside a
    # worker: it catches it and returns a "try-error" object in that slot of
    # the result list instead, leaving the other slots unaffected. Detect
    # that here and re-raise the original error, so a failure behaves the
    # same way (stops riskPlot() with the same message) regardless of
    # n_cores.
    failed <- vapply(results_per_landmark, function(r) inherits(r, "try-error"), logical(1))
    if (any(failed)) {
      stop(conditionMessage(attr(results_per_landmark[[which(failed)[1]]], "condition")), call. = FALSE)
    }
  } else {
    results_per_landmark <- lapply(landmark.time, compute_one_landmark)
  }

  risk.prob.1 = unlist(lapply(results_per_landmark, function(r) r$risk_prob_1))
  risk.prob.2 = if (has_cr) unlist(lapply(results_per_landmark, function(r) r$risk_prob_2)) else c()
  landmark.time.new = unlist(lapply(results_per_landmark, function(r) if (r$keep) r$time.cutoff else NULL))

  if(length(survival_fit_all$form_conditional_cr) != 0){
    # with competing risks
    DP_data = data.frame(time = landmark.time.new, 
                         probType1 = risk.prob.1, 
                         probType2 = risk.prob.2)
    
    if(is.null(bio_i)){
      ## do not plot longitudinal biomarker information
      dp_risk = ggplot() +
        geom_line(data = DP_data, aes(x = time, y = probType1, color = "Event type1"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = probType1, color = "Event type1"), size = 3) + 
        geom_line(data = DP_data, aes(x = time, y = probType2, color = "Event type2"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = probType2, color = "Event type2"), size = 3) +  
        
        scale_color_manual(name = "Lines",
                           values = c("Event type1" = "black", "Event type2" = "red"))   +
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Predicted risk probability") + xlab("Follow-up time") +   
        #0 - 1 black, 1 - 2 red
        geom_vline(xintercept = unlist(getFirst(data_predict_all_pre)[survival_variable])[1], 
                   linetype = "solid", color = unlist(getFirst(data_predict_all_pre)[event_type_variable])[1] + 1, size = 2) + 
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }else{
      ## plot longitudinal biomarker information with risk prediction
      
      scale_prob = 2 * max(na.omit(DP_data_bio$longitudinal))
      dp_risk = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), size = 3) + 
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), size = 3) +  
        
        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",
                                      "Event type1" = "black", "Event type2" = "red"))   +
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        geom_vline(xintercept = unlist(getFirst(data_predict_all_pre)[survival_variable])[1], 
                   linetype = "solid", color = unlist(getFirst(data_predict_all_pre)[event_type_variable])[1] + 1, size = 2) + 
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }

  }else{
    # without competing risks
    DP_data = data.frame(time = landmark.time.new, 
                         probType1 = risk.prob.1)
    
    if(is.null(bio_i)){
      ## do not plot longitudinal biomarker information
      dp_risk = ggplot() +
        geom_line(data = DP_data, aes(x = time, y = probType1, color = "black"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = probType1, color = "black"), size = 3) + 
        
        ylab("Predicted risk probability") + xlab("Follow-up time") +   
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank(),
              legend.position = "none")
      
    }else{
      ## plot longitudinal biomarker information
      scale_prob = 2 * max(na.omit(DP_data_bio$longitudinal))
      dp_risk = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), size = 3) + 

        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",
                                      "Risk probability" = "black"))   +
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        geom_vline(xintercept = unlist(getFirst(data_predict_all_pre)[survival_variable])[1], 
                   linetype = "solid", color = "red", size = 2) + 
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }

    
  }
 return(dp_risk)
}
