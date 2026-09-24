#' Plot of risk using dynamic prediction
#' 
#' @description This function gives the risk prediction plot.
#' 
#' @param data_predict_all_pre This involves a collection of \code{data.frame} objects for
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
#' @param bandcount1 The number of points used to perform the numerical integral, 
#' from the prediction time to the prediction time plus the horizon.
#' @param bandcount2 The number of points used to perform the numerical integral,
#'  from the prediction time to infinity.
#' 
#' @return Plot of risk using dynamic prediction.
#' @export
riskPlot = function(data_predict_all_pre, long_fit_all, survival_fit_all,
                       prediction_time = NULL, bio_i = NULL,
                       horizon, time_variable,
                       survival_variable_all, survival_trans_function,
                       bandcount1 = 10, bandcount2 = 10){
  
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
  
  landmark.time.new = c(); risk.prob.1 = c(); risk.prob.2 = c()
  tt = 0
  for(time.cutoff in landmark.time){
    tt = tt + 1

    data_predict_all = list()
    for(i in seq_len(length(long_fit_all$long_sub_fixed))){
      data_predict_all[[i]] = data_predict_all_pre[[i]][data_predict_all_pre[[i]][time_variable] <= time.cutoff,]
    }
    
    risk.prob = dynamicPrediction(data_predict_all, long_fit_all, survival_fit_all, 
                                  prediction_time = time.cutoff, 
                                  horizon, time_variable,
                                  survival_variable_all, survival_trans_function,
                                  bandcount1, bandcount2)
    
    if(length(survival_fit_all$form_conditional_cr) != 0){
      risk.prob.1 = c(risk.prob.1, risk.prob$risk_prob_1)
      risk.prob.2 = c(risk.prob.2, risk.prob$risk_prob_2)   
    }else{
      risk.prob.1 = c(risk.prob.1, risk.prob$risk_prob_1)
    }
    
    if(length(risk.prob$risk_prob_1) !=0 ) landmark.time.new = c(landmark.time.new, time.cutoff)
    
  }
  
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
