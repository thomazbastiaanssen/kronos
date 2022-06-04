#' Fit a cosinor model and extract relevant parameters
#' @description This is the main workhorse function in the kronos package. It manages the individual functionalities of kronos, including rhythmicity analysis and differential rhythmicity.
#' @param x input data
#' @param measurement A string. Should be the column name containing the meansurement data
#' @param groups A string. Should be the column name containing the group names. If this is kept NULL, kronos will fit a single model without groups.  
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param pairwise A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A kronosOut S4 object containing coefficients and all operations.
#' @export
#' @examples 
#' \dontrun{
#' library(tidyverse)
#' 
#' dat = read.delim("~/Documents/PhD/circadian/diff_rhythm/TPH1_colon.csv", sep = ",")
#' 
#' out = kronos(x = dat, measurement = "Value", time = "Time", groups = "Treatment", period = 24)
#' 
#' out %>% {
#'   ggplot() +
#'     #Light gray rectangle for dark phase
#'     geom_rect(data=NULL,aes(xmin=12,xmax=Inf,ymin=-Inf,ymax=Inf),
#'               fill="lightgray") +
#'     #Whiskers for SEM
#'     geom_errorbar(data = .@input %>%
#'                     group_by(Time, Treatment) %>%
#'                     summarise(stder = std(Value),
#'                               Value = mean(Value, na.rm = T)),
#'                   aes(x = Time,
#'                       ymin = Value - stder,
#'                       ymax = Value + stder,
#'                       group = interaction(Time, Treatment)),
#'                   width = 2, position = position_dodge(0.8)) +
#'     #Big average point per group
#'     geom_point(data = .@input %>%
#'                  group_by(Time, Treatment) %>%
#'                  summarise(Value = mean(Value, na.rm = T)),
#'                aes(x = Time, y = Value, fill = Treatment),
#'                size = 4, shape = 21,
#'                position = position_dodge(0.8))+
#'     #Smaller points per observation
#'     geom_point(data = .@input, aes(x = Time, y = Value, fill = Treatment),
#'                shape = 21,
#'                position = position_dodge(0.8),
#'                alpha = 3/4) +
#'     
#'     #add the line information
#'     geom_line(data = .@to_plot, aes(x = zt, y = pred_value, colour = Treatment)) +
#'     
#'     #Separate by treatment?
#'     #facet_wrap(~Treatment) +
#'     #Fix scales and general layout
#'     scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
#'     theme_bw() +
#'     xlab("Time (h)") +
#'     ylab("Gene expression")
#' }
#' } 
#' 
kronos <- function(x, measurement, groups = NULL, time, period = 24, pairwise = F, verbose = T){
  
  
  stopifnot("The 'time' argument needs to be exactly the same name as one of the columns in input data." = 
              time %in% colnames(x))
  stopifnot("The 'groups' argument needs to be exactly the same name as one of the columns in input data." = 
              groups %in% colnames(x))
  stopifnot("The 'measurement' argument needs to be exactly the same name as one of the columns in input data." = 
              measurement %in% colnames(x))
  
  if(is.null(groups)){x = cbind(x, dummy_groups = 1); groups = "dummy_groups"}
  
  allgroups = sort(unique(x[,groups]))
  
  
  x = cbind(x, get_cos_sine(x = x[,time], paste0(time, "_"), period = period))
  
  fit = fit_cosinor_model(x = x, measurement = measurement, groups = groups, time = time,  verbose = verbose)
  
  vals = get_rhythmic_trace(fit = fit, groups = groups, allgroups = allgroups, time = time, period = period, verbose = verbose)
  
  groupwise = list()
  for(g in 1:length(allgroups)){
    groupwise[[g]] = fit_groupwise_model(x = x, measurement = measurement, 
                                         groups = groups, group = allgroups[g], time = time, period = period, verbose = T)
  }
  bygroup = do.call(rbind, groupwise)
  row.names(bygroup) <- NULL
  
  pairwise_t <- list("pairwise was not activated.")
  if(pairwise){  
    pairwise_t <- list()
    if(verbose){print("Fitting pairwise models")}
    pairwise_t = pairwise_cosinor_model(x = x, measurement = measurement, groups = groups, time = time, pairwise_t = pairwise_t, verbose = verbose)
  }
  
  
  #initialize output container object
  output = new("kronosOut", 
               input      = x,
               fit        = fit,
               to_plot    = vals, 
               ind_fit    = bygroup, 
               pairwise_t = pairwise_t)
  
  return(output)
}

#' Calculate standard error
#' @description calculate standard error for plotting purposes
#' @param x input vector
#' @return standard error of the input vector 
#' @export
std <- function(x, na.rm = TRUE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}

#' Get sine and cosine components
#' @description Based on cosinor and limorhyde packages
get_cos_sine <- function(x, period, colnamePrefix = NULL){
  b = cbind(cos(x/period * 2 * pi), 
            sin(x/period * 2 * pi))
  colnames(b) = c("cos", "sin")
  colnames(b) = paste0(colnamePrefix, colnames(b))
  return(b)
}

#' Fit cosinor model
#' @description Fit cosinor model for totality of data
fit_cosinor_model <- function(x, measurement, groups, time,  verbose = T){

  
  form = as.formula(paste0(measurement, " ~ ", groups, " * (", time, "_cos + ", time, "_sin) -1" ))
  
  #in case of singular group
  if(groups == "dummy_groups"){ form = as.formula(paste0(measurement, " ~ (", time, "_cos + ", time, "_sin)" ))}
  
  if(verbose){print(paste0("Using the following model: ", c(form) ))}
  
  fit = lm(form, data = x)
  
  return(fit)
}


#' Fit cosinor model
#' @description Fit cosinor model for one aspect of data
fit_groupwise_model <- function(x, measurement, groups, group, time, period, verbose = T){
  
  x = x[x[,groups] == group,]
  form = as.formula(paste0(measurement, " ~ (", time, "_cos + ", time, "_sin)" ))
  
  if(verbose){print(paste0("Using the following model: ", c(form) ))}
  
  fit = lm(form, data = x)
  
  fstat = summary(fit)$fstatistic
  p.val = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  r.sq  = summary(fit)$r.squared
  
  vals = expand.grid(groups = group, zt = seq(0, period, 0.25))
  colnames(vals)[1] <- group
  
  vals$zt_cos = cos(vals$zt * 2 * pi / period)
  vals$zt_sin = sin(vals$zt * 2 * pi / period)
  
  
  vals$base     = coef(fit)[1]
  vals$cos_part = coef(fit)[paste0(time, "_cos")] * vals$zt_cos
  vals$sin_part = coef(fit)[paste0(time, "_sin")] * vals$zt_sin
  
  #sum them all up
  vals$pred_value = vals$base + vals$cos_part + vals$sin_part
  
  
  
  ###############This part needs to be fixed
  avg       = mean(vals$pred_value)
  highpoint = which(vals$pred_value == max(vals$pred_value))
  which(vals$pred_value >= avg)
  sin_beta  = coef(fit)[paste0(time, "_sin")]
  cos_beta  = coef(fit)[paste0(time, "_cos")]
  acrophase = (atan2(sin_beta,cos_beta) / 2 / pi * period + period) %% period
  amplitude = sqrt(sin_beta^2 + cos_beta^2)
  
  out_df = data.frame(
    group     = group, 
    p.val     = p.val, 
    r.sq      = r.sq, 
    avg       = avg, 
    acro      = acrophase, 
    amplitude = amplitude)
  
  colnames(out_df)[1] <- groups
  
  return(out_df)
}

#' Get coordinates for plotting rhythmic traces
#' @description Calculates coordinates to plot traces for rhythmic data based on /code{fit_groupwise_model}
#' 
get_rhythmic_trace <- function(fit, groups, allgroups, time, period, verbose = T){
  
  vals = expand.grid(groups = allgroups, zt = seq(0, period, 0.25))
  colnames(vals)[1] <- groups
  
  vals$zt_cos = cos(vals$zt * 2 * pi / period)
  vals$zt_sin = sin(vals$zt * 2 * pi / period)
  
  
  vals$base     = rep(coef(fit)[1:length(allgroups)])
  vals$cos_part = coef(fit)[paste0(time, "_cos")] * vals$zt_cos
  vals$sin_part = coef(fit)[paste0(time, "_sin")] * vals$zt_sin
  vals$cos_int  = rep(c(0, coef(fit)[paste0(groups, allgroups[-1], ":", time, "_cos")])) * vals$zt_cos
  vals$sin_int  = rep(c(0, coef(fit)[paste0(groups, allgroups[-1], ":", time, "_sin")])) * vals$zt_sin
  #sum them all up
  
  #In the case of a singular group; kill the interaction terms:
  if(groups == "dummy_groups"){vals$cos_int = 0; vals$sin_int = 0}
  
  vals$pred_value = vals$base + vals$cos_part + vals$sin_part + vals$cos_int + vals$sin_int
  
  return(vals)
}

#' Fit pairwise cosinor models as some sort of TukeyHSD. 
#' @description Fit cosinor model for subset of data.
#' 
pairwise_cosinor_model <- function(x = x, measurement = measurement, groups = groups, time = time, pairwise_t = pairwise_t, verbose = verbose){
  stopifnot("You need more that two groups in order to sensibly do pairwise comparisons. " = (unique(x[,groups])) > 2)
  combos <- combn(c(unique(x[,groups])), m = 2)
  
  for(c in 1:ncol(combos)){
    group_pair = combos[,c]
    x_pair     = x[x[,groups] %in% group_pair,]
    
    pairwise_t[[c]] <- anova(fit_cosinor_model(x = x_pair, measurement = measurement, groups = groups, time = time,  verbose = T))
    names(pairwise_t)[c] <- paste(group_pair, collapse = " vs ")
  }
  return(pairwise_t)
}



