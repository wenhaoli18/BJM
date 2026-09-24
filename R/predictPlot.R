#' Plot of risk and future biomarker with density using dynamic prediction
#' 
#' @description This function gives the risk and biomarker prediction plot.
#' 
#' @param data.predict.all.one This involves a collection of \code{data.frame} one object for 
#' dynamic prediction and making plots, each corresponding to a distinct longitudinal outcome. 
#' These data frames should contain the variables specified in \code{LongSubFixed} 
#' and \code{LongSubRandom}. Utilizing a list structure 
#' allows for the incorporation of multiple longitudinal outcomes, 
#' each potentially following different measurement protocols. 
#' In instances where all longitudinal outcomes are recorded at identical 
#' time points across patients, a singular \code{data.frame} object may 
#' be used in a \code{list}. It is presumed that each data frame is 
#' structured in a long format.
#' 
#' @param long_fit_all Outputs from the model fitting process using the \code{nlme} package, 
#' encompassing the results and parameters obtained from the analysis.
#' @param survival_fit_all Results and parameters generated from the model fitting 
#' procedure, utilizing the \code{coxph} function. These outputs include the comprehensive 
#' findings and variables derived from the analysis.
#' @param prediction.time Time used to make the prediction.
#' @param horizon Prediction horizon.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survivalVariableAll The name of the transformed time-to-event outcomes variable.
#' @param survivalTransFunction The transformation function used for time-to-event outcomes, 
#' in the order of \code{survivalVariableAll}.
#' @param bandcount1 The number of points used to perform the numerical integral, 
#' from the prediction time to the prediction time plus the horizon.
#' @param bandcount2 The number of points used to perform the numerical integral,
#'  from the prediction time to infinity.
#' @param bandcount3 The number of points used to calculate the probability density function.
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
#' data.survival.fitting =  pbc3[!duplicated(pbc3$id), ]
#' 
#' formMarginalSurv = Surv(years, status3) ~ age + sex
#' formConditionalCR = NULL
#' 
#' survival_fit_all = survivalSub(data.survival.fitting, formMarginalSurv, 
#'                                formConditionalCR)
#' 
#' LongSubFixed = list(
#'   "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
#'   "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
#'   "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
#'   "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
#'   "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
#'   "long6" = platelets ~ year + age + sex + (years)  + (years) * year)
#' 
#' LongSubRandom =list(
#'   "long1" =  ~ year| id,   
#'   "long2" =  ~ year| id,    
#'   "long3" =  ~ year| id,    
#'   "long4" =  ~ year| id,    
#'   "long5" =  ~ year| id,    
#'   "long6" =  ~ year| id)
#' 
#' survivalVariableAll = list(
#'   "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
#' )
#' 
#' survivalTransFunction = list(
#'   fun1 = function(x){abs(x - 1)}, 
#'   fun2 = function(x){abs(x - 3)}, 
#'   fun3 = function(x){abs(x - 5)}, 
#'   fun4 = function(x){abs(x - 7)}
#' )
#' 
#' # Complete case analysis
#' data.fit.all = list()
#' for(i in seq_len(length(LongSubFixed))){
#'   data.fit.all[[i]] = pbc3[pbc3$status3 == 1, ]
#' }
#' 
#' # fitting longitudinal submodel
#' long_fit_all = longitudinalSub(data.fit.all, LongSubFixed, LongSubRandom)
#' 
#' i_PID = 2
#' data.raw.predict.plot = pbc3[pbc3$id == i_PID, ]
#' data.predict.all.pre = list(data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot,
#'                             data.raw.predict.plot, data.raw.predict.plot, data.raw.predict.plot)
#' 
#' # plot biomarker 1 history,  predict future biomarker
#' 
#' predictPlot(data.predict.all.pre, long_fit_all, survival_fit_all, 
#'             prediction.time = 5, bio_his = 1, bio_pred = 1,
#'             horizon = seq(0.5, 3.0, 0.5), time_variable = "year",
#'             survivalVariableAll, survivalTransFunction,
#'            bandcount1 = 10, bandcount2 = 10, bandcount3 = 200)
#'        
#' }
#'     
#' @export
predictPlot = function(data.predict.all.one, long_fit_all, survival_fit_all, 
                    prediction.time = 4, horizon = seq(0.0, 3.0, 0.5), time_variable,
                    survivalVariableAll, survivalTransFunction,
                    bandcount1 = 10, bandcount2 = 10, bandcount3 = 200,
                    bio_his = 1, bio_pred = 1, density = 1){
  
  coxph_fit = survival_fit_all[[1]]
  survival_variable = as.character(formula(coxph_fit)[[2]])[2]
  ### event type variable name
  if(length(survival_fit_all[[4]]) != 0)  event_type_variable = as.character(formula(survival_fit_all[[4]])[[2]])
  
  #name of biomarker
  bio_i_name = as.character(formula(long_fit_all[[3]][[bio_his]])[[2]])
  
  DP_data_bio = data.frame(time = unlist(data.predict.all.one[[bio_his]][time_variable]), 
                           longitudinal = unlist(data.predict.all.one[[bio_his]][bio_i_name]))
  
  DP_data_bio = DP_data_bio[DP_data_bio$time <= prediction.time, ]
  
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
    
    ### data before the prediction time
    data.predict.all = list()
    for(i in seq_len(length(long_fit_all[[3]]))){
      data.predict.all[[i]] = data.predict.all.one[[i]][data.predict.all.one[[i]][time_variable] <= (prediction.time + 1e-8),]
    }
    
    risk.prob = dynamicPrediction(data.predict.all, long_fit_all, survival_fit_all, 
                                  prediction.time, 
                                  horizon = prediction.horizon, time_variable,
                                  survivalVariableAll, survivalTransFunction,
                                  bandcount1, bandcount2)
    
    if(!is.null(bio_pred)){
    Y_predict_all = dynamicPredictionBio(bio_i = bio_his, data.predict.all, long_fit_all, 
                                     survival_fit_all, 
                                     prediction.time, 
                                     horizon = prediction.horizon, time_variable,
                                     survivalVariableAll, survivalTransFunction,
                                     bandcount2, bandcount3)
    
    Y_predict_mode = c(Y_predict_mode, Y_predict_all[[1]])
    Y_all = unlist(Y_predict_all[[3]])
    Y_all_diff = Y_all[2] - Y_all[1]
    my_vector = Y_predict_all[[2]][,1]
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
    
    if(length(survival_fit_all[[4]]) != 0){
      # with competing risks
      risk.prob.1 = c(risk.prob.1, risk.prob[[1]])
      risk.prob.2 = c(risk.prob.2, risk.prob[[2]])   
    }else{
      # without competing risks
      risk.prob.1 = c(risk.prob.1, risk.prob[[1]])
    }
    
  }
  
  ### plot figure
  scale_prob = 2 * max(na.omit(DP_data_bio$longitudinal))
  if(length(survival_fit_all[[4]]) != 0 & is.null(bio_pred)){
      ## with competing risks, without longitudinal biomarker information
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
      
    } else if (length(survival_fit_all[[4]]) != 0 & !is.null(bio_pred) & !is.null(density)){
      
      ## with competing risks, with longitudinal biomarker information and density plots
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }else if (length(survival_fit_all[[4]]) != 0 & !is.null(bio_pred) & is.null(density)){
      ## with competing risks, with longitudinal biomarker information without density plots
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
    }else if (length(survival_fit_all[[4]]) == 0 & is.null(bio_pred) ){
      ## without competing risks, without longitudinal biomarker information
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
    }else if (length(survival_fit_all[[4]]) == 0 & !is.null(bio_pred) & !is.null(density)){
      ## without competing risks, with longitudinal biomarker information with density plots
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
        geom_hline(yintercept = c(0, scale_prob/5, scale_prob/5*2, scale_prob/5*3, 
                                  scale_prob/5*4, scale_prob), 
                   linetype = "dotted", color = "pink", size = 1.2) +
        theme_bw(base_size = 25)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank()) 
      
      
    }else if (length(survival_fit_all[[4]]) == 0 & !is.null(bio_pred) & is.null(density)){
      ## without competing risks, with longitudinal biomarker information without density plots
      
      DP_data = data.frame(time = prediction.time + horizon, 
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
        geom_vline(xintercept = prediction.time, linetype = "solid", color = "brown", size = 1) + 
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
