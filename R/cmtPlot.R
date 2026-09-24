#' Plot conditional mean trajectories (CMT)
#' 
#' @description This function generates the Conditional Mean Trajectories (CMT) plot. 
#' All patients in this plot experience events at the same time point, 
#' specified by \code{condi_time2event}. Several evenly spaced time points between 
#' the baseline and \code{condi_time2event} are selected for plotting. Each point is 
#' calculated using the mean value of all patients' biomarker values at that time point. 
#' The interval between two time points is defined by \code{interval_time}
#' 
#' @param data_plot_all A \code{data.frame} that includes the biomarker used for plotting. 
#' It is utilized to plot conditional mean trajectories (CMT).
#' @param condi_time2event Conditional event time, indicating that all patients should
#' have events at this time in the plot. If \code{NULL}, it defaults to the midpoint
#' of the observed range of \code{time_variable} in \code{data_plot_all}.
#' @param event_type_variable Competing risks variable indicator name. Set to NULL if 
#' there are no competing risks.
#' @param event_type A vector containing the names of all event types.
#' @param bio_variable Name of the biomarker variable used for plotting.
#' @param time_variable The name of time variable in linear mixed model.
#' @param survival_variable Name of the time-to-event outcomes variable.
#' @param interval_time The time interval between two time points. Time points are 
#' plotted within the baseline to event time.
#' @param id_variable Name of the patient ID column in \code{data_plot_all}. Default is \code{"id"}.
#' 
#' @return Conditional mean trajectories plot.
#' 
#' @examples 
#' 
#' # example without competing risks
#' 
#' data(pbc3)
#' 
#' pbc.cmt <- cmtPlot(data_plot_all = pbc3, condi_time2event = 5, 
#'    event_type_variable = NULL, event_type = NULL,
#'    bio_variable = "serBilir", time_variable = "year", 
#'    survival_variable = "years", 
#'    interval_time = 1/12
#' )
#' 
#' pbc.cmt
#' 
#' @examples 
#' 
#' # example with competing risks
#' 
#' data(pbc3)
#' 
#' data_plot_all = pbc3[!is.na(pbc3$status4),]
#' 
#' pbc.cmt.cr <- cmtPlot(data_plot_all, condi_time2event = 5, 
#'    event_type_variable = 'status4', event_type = c("0", "1"),
#'    bio_variable = "albumin", time_variable = "year", 
#'    survival_variable = "years", 
#'    interval_time = 1/4
#' )
#' 
#' pbc.cmt.cr
#'
#' @export
cmtPlot = function(data_plot_all, condi_time2event, event_type_variable, event_type,
                   bio_variable, time_variable, survival_variable, interval_time = 1/12,
                   id_variable = "id"){

  assert_data_frame(data_plot_all, "data_plot_all")
  assert_string(bio_variable, "bio_variable")
  assert_string(time_variable, "time_variable")
  assert_string(survival_variable, "survival_variable")
  assert_string(id_variable, "id_variable")
  assert_vars_in_data(c(bio_variable, time_variable, survival_variable, id_variable), data_plot_all,
                       "bio_variable/time_variable/survival_variable/id_variable", "data_plot_all")
  if (!is.null(event_type_variable)) {
    assert_string(event_type_variable, "event_type_variable")
    if (length(event_type) < 2) {
      stop("`event_type` must have at least 2 elements (one per competing event type) when `event_type_variable` is supplied.", call. = FALSE)
    }
  }
  if (!is.null(condi_time2event)) {
    assert_scalar_numeric(condi_time2event, "condi_time2event")
  }
  assert_scalar_numeric(interval_time, "interval_time", positive = TRUE)

  ### A sequence of conditional time-to-event times
  if(!is.null(condi_time2event)){
    ### if condi_time2event is pre-defined
    condi_time2event_seq = condi_time2event
  }else{
    ### Pick the middel point of the time duration
    condi_time2event_seq = ceiling(1/2 * min(unlist(data_plot_all[time_variable])) +
                                     1/2 * max(unlist(data_plot_all[time_variable])))
    ### otherwise generate a sequence automatically
    #condi_time2event_seq = seq(ceiling(min(unlist(plot_data[time_variable]))),
    #                           ceiling(max(unlist(plot_data[time_variable]))),
    #                           by = interval_time_1 )[-1]
  }
  
  if(!is.null(event_type_variable)){
  ### with competing risks
    
  ### Event type 1
  plot_data = data_plot_all[data_plot_all[event_type_variable] == event_type[1],]
  cluster_event_value = list() 
    
  #### number of patients has events at certain period
  event1_number = c()
  event_year_i = 0
  for(event_year in condi_time2event_seq){
    event_year_i = event_year_i + 1
    ### sample data who has events in the interval of [event_year - 1, event_year] 
    plot_data_event_year = plot_data[plot_data[survival_variable] <= event_year & 
                                     plot_data[survival_variable] >= event_year - 1,]
    event1_number[event_year_i] = dim(plot_data_event_year[!duplicated(plot_data_event_year[[id_variable]]), ])[1]
    if(dim(plot_data_event_year)[1] == 0) {
      ### if no data in plot_data_event_year, use previous event year
      cluster_event_value[[event_year_i]] = cluster_event_value[[event_year_i - 1]]
    }else{
      cluster_event_value[[event_year_i]] = vector()
      ### calculate the mean of biomarker value at a sequence of points
      for(year in seq(from = ceiling(min(plot_data_event_year[time_variable])), 
                      to = event_year, by = interval_time) ){ 
        plot_data_interval = plot_data_event_year[plot_data_event_year[time_variable] <= year & 
                                        plot_data_event_year[time_variable] >= year - interval_time,]
        if(dim(plot_data_interval)[1] == 0) {
          mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
        }else{
          mean_value  = mean(unlist(plot_data_interval[bio_variable]))
          if(is.na(mean_value)){
            mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
          }
        }
        
        cluster_event_value[[event_year_i]] =  c(cluster_event_value[[event_year_i]], mean_value)
      }
    }
    
  }
  cluster_event_value_event1 = cluster_event_value
  
  #### Event type 2
  plot_data = data_plot_all[data_plot_all[event_type_variable] == event_type[2],]
  cluster_event_value = list() 
  
  #### number of patients has events at certain period
  event2_number = c()
  event_year_i = 0
  for(event_year in condi_time2event_seq){
    event_year_i = event_year_i + 1
    ### sample data who has events in the interval of [event_year - 1, event_year] 
    plot_data_event_year = plot_data[plot_data[survival_variable] <= event_year & 
                                       plot_data[survival_variable] >= event_year - 1,]
    event2_number[event_year_i] = dim(plot_data_event_year[!duplicated(plot_data_event_year[[id_variable]]), ])[1]
    if(dim(plot_data_event_year)[1] == 0) {
      ### if no data in plot_data_event_year, use previous event year
      cluster_event_value[[event_year_i]] = cluster_event_value[[event_year_i - 1]]
    }else{
      cluster_event_value[[event_year_i]] = vector()
      ### calculate the mean of biomarker value at a sequence of points
      for(year in seq(from = ceiling(min(plot_data_event_year[time_variable])), 
                      to = event_year, by = interval_time) ){ 
        plot_data_interval = plot_data_event_year[plot_data_event_year[time_variable] <= year & 
                                                    plot_data_event_year[time_variable] >= year - interval_time,]
        if(dim(plot_data_interval)[1] == 0) {
          mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
        }else{
          mean_value  = mean(unlist(plot_data_interval[bio_variable]))
          if(is.na(mean_value)){
            mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
          }
        }
        
        cluster_event_value[[event_year_i]] =  c(cluster_event_value[[event_year_i]], mean_value)
      }
    }
    
  }
  cluster_event_value_event2 = cluster_event_value
  
  for(event_year_i in seq_len(length(condi_time2event_seq))){
    DP_event1 = data.frame(time = rep(1: length(cluster_event_value_event1[[event_year_i]])), 
                         probEvent = cluster_event_value_event1[[event_year_i]])
    DP_event2 = data.frame(time = rep(1: length(cluster_event_value_event2[[event_year_i]])), 
                          probEvent = cluster_event_value_event2[[event_year_i]])
    
    nam <- paste("Plot_p", event_year_i, sep = "") 
    assign(nam,  ggplot() +
             geom_point(data = DP_event1, aes(x = time, y = probEvent, color = "Event1"), size = 3) +  
             geom_smooth(data = DP_event1, aes(x = time, y = probEvent, color = "Event1"), method = "lm") +
             geom_point(data = DP_event2, aes(x = time, y =  probEvent, color = "Event2"), size = 3) +
             geom_smooth(data = DP_event2, aes(x = time, y = probEvent, color = "Event2"), method = "lm") +
             
             scale_color_manual(name = "Lines",
                                values = c("Event1" = "red", "Event2" = "black"))   +
             ylab("Biomarker value") + xlab("Follow-up time") +   
             theme_bw(base_size = 22)+
             theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   plot.background = element_blank()) )
  }
  
  }else{
    # without competing risks
    
    plot_data = data_plot_all
    cluster_event_value = list() 
    #### number of patients has events at certain period
    event1_number = c()
    event_year_i = 0
    for(event_year in condi_time2event_seq){
      event_year_i = event_year_i + 1
      ### sample data who has events in the interval of [event_year - 1, event_year] 
      plot_data_event_year = plot_data[plot_data[survival_variable] <= event_year & 
                                         plot_data[survival_variable] >= event_year - 1,]
      event1_number[event_year_i] = dim(plot_data_event_year[!duplicated(plot_data_event_year[[id_variable]]), ])[1]
      if(dim(plot_data_event_year)[1] == 0) {
        ### if no data in plot_data_event_year, use previous event year
        cluster_event_value[[event_year_i]] = cluster_event_value[[event_year_i - 1]]
      }else{
        cluster_event_value[[event_year_i]] = vector()
        ### calculate the mean of biomarker value at a sequence of points
        for(year in seq(from = ceiling(min(plot_data_event_year[time_variable])), 
                        to = event_year, by = interval_time) ){ 
          plot_data_interval = plot_data_event_year[plot_data_event_year[time_variable] <= year & 
                                                      plot_data_event_year[time_variable] >= year - interval_time,]
          if(dim(plot_data_interval)[1] == 0) {
            mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
          }else{
            mean_value  = mean(unlist(plot_data_interval[bio_variable]))
            if(is.na(mean_value)){
              mean_value  = cluster_event_value[[event_year_i]][length(cluster_event_value[[event_year_i]])]
            }
          }
          
          cluster_event_value[[event_year_i]] =  c(cluster_event_value[[event_year_i]], mean_value)
        }
      }
      
    }
    cluster_event_value_event1 = cluster_event_value
    
    for(event_year_i in seq_len(length(condi_time2event_seq))){
      DP_event1 = data.frame(time = rep(1: length(cluster_event_value_event1[[event_year_i]])), 
                             probEvent = cluster_event_value_event1[[event_year_i]])
      
      nam <- paste("Plot_p", event_year_i, sep = "") 
      assign(nam,  ggplot() +
               geom_point(data = DP_event1, aes(x = time, y = probEvent, color = "Event1"), size = 3) +  
               geom_smooth(data = DP_event1, aes(x = time, y = probEvent, color = "Event1"), method = "lm") +
               
               scale_color_manual(name = "Lines",
                                  values = c("Event1" = "red"))   +
               ylab("Biomarker value") + xlab("Follow-up time") +   
               theme_bw(base_size = 22)+
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_blank(),
                     legend.position = "none") )
    }
    
    }
  
  return(Plot_p1)
}

