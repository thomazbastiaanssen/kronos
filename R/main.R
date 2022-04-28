#' Fit a cosinor model and extract relevant parameters
#' @description This is the main workhorse function in the kronos package. It manages the individual functionalities of kronos, including rhythmicity analysis and differential rhythmicity.
#' @param x input data
#' @param measurement A string. Should be the column name containing the meansurement data
#' @param groups A string. Should be the column name containing the group names.  
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A list containing coefficients and all operations.
#' @export 
#' 
pairwise_rhythm <- function(x, measurement, groups, time, period = 24, verbose = T){
  stopifnot("The 'time' argument needs to be exactly the same name as one of the columns in input data." = 
              time %in% colnames(x))
  stopifnot("The 'groups' argument needs to be exactly the same name as one of the columns in input data." = 
              groups %in% colnames(x))
  stopifnot("The 'measurement' argument needs to be exactly the same name as one of the columns in input data." = 
              measurement %in% colnames(x))
  
  allgroups = sort(unique(x[,groups]))
  
  
  x = cbind(x, get_cos_sine(x = x[,time], paste0(time, "_"), period = period))
  
  fit = fit_cosinor_model(x = x, measurement = measurement, groups = groups, time = time,  verbose = verbose)
  
  #tidyfit = broom::tidy(fit)
  vals = get_rhythmic_trace(fit = fit, groups = groups, allgroups = allgroups, time = time, period = period, verbose = verbose)
  
  groupwise = list()
  for(g in 1:length(allgroups)){
    groupwise[[g]] = fit_groupwise_model(x = x, measurement = measurement, groups = groups, group = allgroups[g], time = time, verbose = T)
  }
  bygroup = do.call(rbind, groupwise)
  row.names(bygroup) <- NULL
  ###############Plot
  
  
  return(list(input   = x, 
              fit     = fit, 
              to_plot = vals, 
              ind_fit = bygroup))
}
#' Calculate standard error
#' @description calculate standard error for plotting purposes
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
  
  if(verbose){print(paste0("Using the following model: ", c(form) ))}
  
  fit = lm(form, data = x)
  
  return(fit)
}


#' Fit cosinor model
#' @description Fit cosinor model for one aspect of data
fit_groupwise_model <- function(x, measurement, groups, group, time,  verbose = T){
  
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
  vals$pred_value = vals$base + vals$cos_part + vals$sin_part + vals$cos_int + vals$sin_int
  
  return(vals)
}



