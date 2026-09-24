#' Plot of risk and future biomarker with density using dynamic prediction
#' 
#' @description This function gives the risk and biomarker prediction plot.
#' 
#' @param data_predict_all_one This involves a collection of \code{data.frame} one object for
#' dynamic prediction and making plots, each corresponding to a distinct longitudinal outcome. 
#' These data frames should contain the variables specified in \code{long_sub_fixed} 
#' and \code{long_sub_random}. Utilizing a list structure 
#' allows for the incorporation of multiple longitudinal outcomes, 
#' each potentially following different measurement protocols. 
#' In instances where all longitudinal outcomes are recorded at identical
#' time points across patients, a singular \code{data.frame} object may
#' be used in a \code{list}. Alternatively, a single bare \code{data.frame}
#' (not wrapped in a list) may be supplied directly; it is then reused for
#' every longitudinal outcome. It is presumed that each data frame is
#' structured in a long format.
#'
#' @param long_fit_all Outputs from the model fitting process using the \code{nlme} package,
#' encompassing the results and parameters obtained from the analysis.
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param prediction_time Time used to make the prediction.
#' @param horizon Prediction horizon.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable_all The name of the transformed time-to-event outcomes variable.
#' @param survival_trans_function The transformation function used for time-to-event outcomes, 
#' in the order of \code{survival_variable_all}.
#' @param bandcount1 The number of grid points spanning the prediction window,
#' from \code{prediction_time} to \code{prediction_time + horizon}. Larger
#' values give a more accurate but slower estimate. Defaults to \code{"auto"},
#' which resolves it once, before looping over \code{horizon} (using the
#' largest requested horizon as a representative probe), by doubling from a
#' built-in starting value until the predicted risk stabilizes; see
#' \code{\link{dynamicPrediction}}'s \code{bandcount1} for details of that
#' search. The resolved value is then reused, fixed, for every point in
#' \code{horizon} -- it is not re-searched on every iteration.
#' @param bandcount2 The number of grid points used to approximate
#' integrating out to infinity when normalizing the predicted risk/density.
#' A wider follow-up range needs a larger \code{bandcount2} to keep the
#' grid spacing comparable. Defaults to \code{"auto"}; resolved the same way
#' as \code{bandcount1} (jointly with it, when both are \code{"auto"}).
#' @param bandcount3 The number of points in the candidate-biomarker-value
#' grid used to build the predicted density curve; controls the resolution
#' of the density, not a time integral. Defaults to \code{"auto"}; resolved
#' the same way, but only when \code{bio_pred} is non-\code{NULL} (it is
#' unused otherwise).
#'
#' Pass explicit numbers instead of \code{"auto"} for full manual control, or
#' use \code{checkBandcountConvergence()} (applied to \code{dynamicPrediction()}/
#' \code{dynamicPredictionBio()} directly) to inspect the convergence behavior
#' yourself. See also \code{vignette("BJM-intro", package = "BJM")} for
#' further guidance on choosing \code{bandcount1}/\code{bandcount2}/
#' \code{bandcount3}.
#'
#' @param bio_his Which biomarker history will be plotted
#' @param bio_pred Indicator, predict future biomarker or not, if NULL do not predict
#' @param density Indicator, plot future biomarker density or not, if NULL do not plot
#' @return Plot of risk and future biomarker with density using dynamic prediction.
#' 
#' @examples 
#' \donttest{
#' data(pbc3)
#' 
#' data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]
#' 
#' form_marginal_surv = Surv(years, status3) ~ age + sex
#' form_conditional_cr = NULL
#' 
#' survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv, 
#'                                form_conditional_cr)
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
#' i_PID = 2
#' data.raw.predict.plot = pbc3[pbc3$id == i_PID, ]
#' data_predict_all_one = list(data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot,
#'                             data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot)
#'
#' # plot biomarker 1 history,  predict future biomarker
#'
#' predictPlot(data_predict_all_one, long_fit_all, survival_fit_all,
#'             prediction_time = 5, bio_his = 1, bio_pred = 1,
#'             horizon = seq(0.5, 3.0, 0.5), time_variable = "year",
#'             survival_variable_all, survival_trans_function,
#'            bandcount1 = 10, bandcount2 = 10, bandcount3 = 200)
#'        
#' }
#'     
#' @export
predictPlot = function(data_predict_all_one, long_fit_all, survival_fit_all,
                    prediction_time = 4, horizon = seq(0.0, 3.0, 0.5), time_variable,
                    survival_variable_all, survival_trans_function,
                    bandcount1 = "auto", bandcount2 = "auto", bandcount3 = "auto",
                    bio_his = 1, bio_pred = 1, density = 1){

  assert_class(long_fit_all, "longitudinalSub.BJM", "long_fit_all", "longitudinalSub")
  assert_class(survival_fit_all, "survivalSub.BJM", "survival_fit_all", "survivalSub")
  assert_index(bio_his, length(long_fit_all$lfit), "bio_his", "longitudinal outcomes in long_fit_all")
  assert_data_list(data_predict_all_one, "data_predict_all_one", length(long_fit_all$lfit), allow_bare_df = TRUE)
  if (!is.list(data_predict_all_one) || is.data.frame(data_predict_all_one)) {
    data_predict_all_one <- rep(list(data_predict_all_one), each = length(long_fit_all$lfit))
  }
  assert_scalar_numeric(prediction_time, "prediction_time")
  assert_string(time_variable, "time_variable")
  assert_bandcount(bandcount1, "bandcount1")
  assert_bandcount(bandcount2, "bandcount2")
  assert_bandcount(bandcount3, "bandcount3")
  assert_survival_trans(survival_variable_all, survival_trans_function, probe_value = prediction_time)
  for (i in seq_along(data_predict_all_one)) {
    assert_vars_in_data(time_variable, data_predict_all_one[[i]], "time_variable",
                         sprintf("data_predict_all_one[[%d]]", i))
  }

  coxph_fit = survival_fit_all$coxph_fit
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]
  ### event type variable name
  if(length(survival_fit_all$form_conditional_cr) != 0)  event_type_variable = as.character(formula(survival_fit_all$form_conditional_cr)[[2]])
  
  #name of biomarker
  bio_i_name = as.character(formula(long_fit_all$long_sub_fixed[[bio_his]])[[2]])
  
  DP_data_bio = data.frame(time = unlist(data_predict_all_one[[bio_his]][time_variable]),
                           longitudinal = unlist(data_predict_all_one[[bio_his]][bio_i_name]))
  
  DP_data_bio = DP_data_bio[DP_data_bio$time <= prediction_time, ]

  ### data before the prediction time -- this does not depend on
  ### prediction.horizon, so it is built once here rather than inside the
  ### loop below.
  data_predict_all = list()
  for(i in seq_len(length(long_fit_all$long_sub_fixed))){
    data_predict_all[[i]] = data_predict_all_one[[i]][data_predict_all_one[[i]][time_variable] <= (prediction_time + 1e-8),]
  }

  ### bandcount1/bandcount2/bandcount3 = "auto" (the default): resolve them
  ### once here, using the largest requested horizon as a representative
  ### probe, rather than re-running the auto-tuning search on every horizon
  ### in the loop below (see auto_tune_bandcount()).
  auto_names_1_2 <- c("bandcount1", "bandcount2")[c(identical(bandcount1, "auto"), identical(bandcount2, "auto"))]
  if (length(auto_names_1_2) > 0) {
    probe_args <- list(data_predict_all = data_predict_all, long_fit_all = long_fit_all,
                        survival_fit_all = survival_fit_all, prediction_time = prediction_time,
                        horizon = max(horizon), time_variable = time_variable,
                        survival_variable_all = survival_variable_all,
                        survival_trans_function = survival_trans_function,
                        bandcount1 = bandcount1, bandcount2 = bandcount2)
    resolved_1_2 <- auto_tune_bandcount(dynamicPrediction, probe_args, auto_names_1_2)$bandcount
    if (!is.null(resolved_1_2$bandcount1)) bandcount1 <- resolved_1_2$bandcount1
    if (!is.null(resolved_1_2$bandcount2)) bandcount2 <- resolved_1_2$bandcount2
  }
  if (!is.null(bio_pred) && identical(bandcount3, "auto")) {
    probe_args <- list(bio_i = bio_his, data_predict_all = data_predict_all, long_fit_all = long_fit_all,
                        survival_fit_all = survival_fit_all, prediction_time = prediction_time,
                        horizon = max(horizon), time_variable = time_variable,
                        survival_variable_all = survival_variable_all,
                        survival_trans_function = survival_trans_function,
                        bandcount2 = bandcount2, bandcount3 = bandcount3)
    resolved_3 <- auto_tune_bandcount(dynamicPredictionBio, probe_args, "bandcount3")$bandcount
    bandcount3 <- resolved_3$bandcount3
  }

  ### risk predicted probability
  risk.prob.1 = c(); risk.prob.2 = c()
  ### mode prediction
  Y_predict_mode = c()
  ### quantiles prediction
  Y_predict_quantile_1_10 = c()
  Y_predict_quantile_2_10 = c()
  Y_predict_quantile_3_10 = c()
  Y_predict_quantile_4_10 = c()
  Y_predict_quantile_5_10 = c()
  Y_predict_quantile_6_10 = c()
  Y_predict_quantile_7_10 = c()
  Y_predict_quantile_8_10 = c()
  Y_predict_quantile_9_10 = c()
  tt = 0
  for(prediction.horizon in horizon){
    tt = tt + 1

    risk.prob = dynamicPrediction(data_predict_all, long_fit_all, survival_fit_all,
                                  prediction_time, 
                                  horizon = prediction.horizon, time_variable,
                                  survival_variable_all, survival_trans_function,
                                  bandcount1, bandcount2)
    
    if(!is.null(bio_pred)){
    Y_predict_all = dynamicPredictionBio(bio_i = bio_his, data_predict_all, long_fit_all, 
                                     survival_fit_all, 
                                     prediction_time, 
                                     horizon = prediction.horizon, time_variable,
                                     survival_variable_all, survival_trans_function,
                                     bandcount2, bandcount3)
    
    Y_predict_mode = c(Y_predict_mode, Y_predict_all$Y_predict)
    Y_all = unlist(Y_predict_all$Y_all)
    Y_all_diff = Y_all[2] - Y_all[1]
    my_vector = Y_predict_all$Y_density[,1]
    index <- which(cumsum(my_vector) >= 0.1 * 1/Y_all_diff)[1]
    Y_predict_quantile_1_10 <- c(Y_predict_quantile_1_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.2 * 1/Y_all_diff)[1]
    Y_predict_quantile_2_10 <- c(Y_predict_quantile_2_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.3 * 1/Y_all_diff)[1]
    Y_predict_quantile_3_10 <- c(Y_predict_quantile_3_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.4 * 1/Y_all_diff)[1]
    Y_predict_quantile_4_10 <- c(Y_predict_quantile_4_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.5 * 1/Y_all_diff)[1]
    Y_predict_quantile_5_10 <- c(Y_predict_quantile_5_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.6 * 1/Y_all_diff)[1]
    Y_predict_quantile_6_10 <- c(Y_predict_quantile_6_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.7 * 1/Y_all_diff)[1]
    Y_predict_quantile_7_10 <- c(Y_predict_quantile_7_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.8 * 1/Y_all_diff)[1]
    Y_predict_quantile_8_10 <- c(Y_predict_quantile_8_10, Y_all[index])
    index <- which(cumsum(my_vector) >= 0.9 * 1/Y_all_diff)[1]
    Y_predict_quantile_9_10 <- c(Y_predict_quantile_9_10, Y_all[index])
    }
    
    if(length(survival_fit_all$form_conditional_cr) != 0){
      # with competing risks
      risk.prob.1 = c(risk.prob.1, risk.prob$risk_prob_1)
      risk.prob.2 = c(risk.prob.2, risk.prob$risk_prob_2)   
    }else{
      # without competing risks
      risk.prob.1 = c(risk.prob.1, risk.prob$risk_prob_1)
    }
    
  }
  
  ### plot figure
  scale_prob = 2 * max(na.omit(DP_data_bio$longitudinal))
  if(length(survival_fit_all$form_conditional_cr) != 0 & is.null(bio_pred)){
      ## with competing risks, without longitudinal biomarker information
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1, 
                           probType2 = risk.prob.2)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), size = 3) + 
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), size = 3)   + 
        
        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Event type1" = "black", "Event type2" = "red"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
      
    } else if (length(survival_fit_all$form_conditional_cr) != 0 & !is.null(bio_pred) & !is.null(density)){
      
      ## with competing risks, with longitudinal biomarker information and density plots
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1, 
                           probType2 = risk.prob.2,
                           predMode = Y_predict_mode,
                           predQuan1 = Y_predict_quantile_1_10,
                           predQuan2 = Y_predict_quantile_2_10,
                           predQuan3 = Y_predict_quantile_3_10,
                           predQuan4 = Y_predict_quantile_4_10,
                           predQuan5 = Y_predict_quantile_5_10,
                           predQuan6 = Y_predict_quantile_6_10,
                           predQuan7 = Y_predict_quantile_7_10,
                           predQuan8 = Y_predict_quantile_8_10,
                           predQuan9 = Y_predict_quantile_9_10)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), linetype = "solid", size = 4) + 
        geom_point(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal")) +
        geom_text(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), label = "L", size = 6, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predQuan1, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan2, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan3, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan4, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan5, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan6, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan7, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan8, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan9, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        
        # Add the shaded area between the two lines
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan1, ymax = predQuan2), fill = "#00FFCC", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan2, ymax = predQuan3), fill = "#33CC99", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan3, ymax = predQuan4), fill = "#009956", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan4, ymax = predQuan5), fill = "#003300", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan5, ymax = predQuan6), fill = "#006600", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan6, ymax = predQuan7), fill = "#009956", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan7, ymax = predQuan8), fill = "#33CC99", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan8, ymax = predQuan9), fill = "#00FFCC", alpha = 0.5) + 
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), size = 3) + 
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), size = 3)   + 
        
        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Event type1" = "black", "Event type2" = "red"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +  
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }else if (length(survival_fit_all$form_conditional_cr) != 0 & !is.null(bio_pred) & is.null(density)){
      ## with competing risks, with longitudinal biomarker information without density plots
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1, 
                           probType2 = risk.prob.2,
                           predMode = Y_predict_mode,
                           predQuan1 = Y_predict_quantile_1_10,
                           predQuan2 = Y_predict_quantile_2_10,
                           predQuan3 = Y_predict_quantile_3_10,
                           predQuan4 = Y_predict_quantile_4_10,
                           predQuan5 = Y_predict_quantile_5_10,
                           predQuan6 = Y_predict_quantile_6_10,
                           predQuan7 = Y_predict_quantile_7_10,
                           predQuan8 = Y_predict_quantile_8_10,
                           predQuan9 = Y_predict_quantile_9_10)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), linetype = "solid", size = 4) + 
        geom_point(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal")) +
        geom_text(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), label = "L", size = 6, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Event type1"), size = 3) + 
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType2, color = "Event type2"), size = 3)   + 
        
        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Event type1" = "black", "Event type2" = "red"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }else if (length(survival_fit_all$form_conditional_cr) == 0 & is.null(bio_pred) ){
      ## without competing risks, without longitudinal biomarker information
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), size = 3) + 

        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Risk probability" = "black"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
    }else if (length(survival_fit_all$form_conditional_cr) == 0 & !is.null(bio_pred) & !is.null(density)){
      ## without competing risks, with longitudinal biomarker information with density plots
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1,
                           predMode = Y_predict_mode,
                           predQuan1 = Y_predict_quantile_1_10,
                           predQuan2 = Y_predict_quantile_2_10,
                           predQuan3 = Y_predict_quantile_3_10,
                           predQuan4 = Y_predict_quantile_4_10,
                           predQuan5 = Y_predict_quantile_5_10,
                           predQuan6 = Y_predict_quantile_6_10,
                           predQuan7 = Y_predict_quantile_7_10,
                           predQuan8 = Y_predict_quantile_8_10,
                           predQuan9 = Y_predict_quantile_9_10)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), linetype = "solid", size = 4) + 
        geom_point(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal")) +
        geom_text(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), label = "L", size = 6, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predQuan1, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan2, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan3, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan4, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan5, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan6, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan7, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan8, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        geom_line(data = DP_data, aes(x = time, y = predQuan9, color = "Longitudinal"), linetype = "dashed", size = 0.4) + 
        
        # Add the shaded area between the two lines
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan1, ymax = predQuan2), fill = "#00FFCC", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan2, ymax = predQuan3), fill = "#33CC99", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan3, ymax = predQuan4), fill = "#009956", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan4, ymax = predQuan5), fill = "#003300", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan5, ymax = predQuan6), fill = "#006600", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan6, ymax = predQuan7), fill = "#009956", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan7, ymax = predQuan8), fill = "#33CC99", alpha = 0.5) + 
        geom_ribbon(data = DP_data, aes(x = time, ymin = predQuan8, ymax = predQuan9), fill = "#00FFCC", alpha = 0.5) + 
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), size = 3) + 

        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Risk probability" = "black"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
    }else if (length(survival_fit_all$form_conditional_cr) == 0 & !is.null(bio_pred) & is.null(density)){
      ## without competing risks, with longitudinal biomarker information without density plots
      
      DP_data = data.frame(time = prediction_time + horizon, 
                           probType1 = risk.prob.1,
                           predMode = Y_predict_mode,
                           predQuan1 = Y_predict_quantile_1_10,
                           predQuan2 = Y_predict_quantile_2_10,
                           predQuan3 = Y_predict_quantile_3_10,
                           predQuan4 = Y_predict_quantile_4_10,
                           predQuan5 = Y_predict_quantile_5_10,
                           predQuan6 = Y_predict_quantile_6_10,
                           predQuan7 = Y_predict_quantile_7_10,
                           predQuan8 = Y_predict_quantile_8_10,
                           predQuan9 = Y_predict_quantile_9_10)
      
      dp_plot = ggplot() +
        geom_line(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) + 
        geom_point(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal")) +
        geom_text(data = DP_data_bio, aes(x = time, y = longitudinal, color = "Longitudinal"), label = "L", size = 4, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), linetype = "solid", size = 4) + 
        geom_point(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal")) +
        geom_text(data = DP_data, aes(x = time, y = predMode, color = "Longitudinal"), label = "L", size = 6, vjust = -0.5)    +
        
        geom_line(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), linetype = "solid", size = 1) +
        geom_point(data = DP_data, aes(x = time, y = scale_prob * probType1, color = "Risk probability"), size = 3) + 
        
        scale_color_manual(name = "Lines",
                           values = c("Longitudinal" = "green",  
                                      "Risk probability" = "black"),
                           guide = guide_legend(title = "Linetype"))  + 
        scale_y_continuous(sec.axis = sec_axis(~./scale_prob, name="Risk Probabilities")) + 
        ylab("Longitudinal biomarker") + xlab("Follow-up time") +   
        scale_x_continuous(breaks = seq(0, 15, 1))  +
        geom_vline(xintercept = prediction_time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
    }
    
  return(dp_plot)
}
