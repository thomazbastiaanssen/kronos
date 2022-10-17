#' Fit a cosinor model and extract relevant parameters
#' @description This is the main workhorse function in the kronos package. It manages the individual functionalities of kronos, including rhythmicity analysis and differential rhythmicity.
#' @param data input data
#' @param formula A formula. Use the time() function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param pairwise A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A kronosOut S4 object containing coefficients and all operations.
#' @export
kronos_test <- function(formula, data, time = NULL, period = 24, verbose = T, pairwise = T){
  
#Identify time component
time <- get_time_test(formula = formula, time = time, data = data, verbose = verbose)  

#Set up sine and cosine component
data <- cbind(data, get_cos_sine_test(data = data[,time], colnamePrefix = paste0(time, "_"), period = period))

#Fit the general model
fit  <- fit_cosinor_model_test(formula = formula, data = data, time = time, verbose = verbose, for_pw = F)

#Predict values for plotting
vals <- kronos_predict_test(fit = fit, period = period, time = time, verbose = verbose)

fit$model$unique_name <- paste(fit$model[,names(fit$xlevels)])

groupwise = vector(mode = "list", length = length(unique(fit$model$unique_name)))

for(g in 1:length(groupwise)){
  groupwise[[g]] = fit_groupwise_model_test(data = fit$model, 
                                       group = unique(fit$model$unique_name)[g], 
                                       time = time, period = period, verbose = T)
}
bygroup = do.call(rbind, groupwise)
row.names(bygroup) <- NULL

pairwise_t <- list("pairwise was not activated.")
if(pairwise){  
  pairwise_t <- list()
  if(verbose){print("Fitting pairwise models")}
  pairwise_t = pairwise_cosinor_model_test(data = fit$model, time = time, pairwise_t = pairwise_t, verbose = verbose)
}


#initialize output container object
output = new("kronosOut", 
             input      = data,
             fit        = fit,
             to_plot    = vals,
             ind_fit    = bygroup, 
             pairwise_t = pairwise_t)

return(output)
}

#' Figure out what variable represents time. Called by main kronos function. 
#' @description Extracts time from the formula and from the time argument. Also handles inconsistencies. 
#' @param formula A formula. Use the time() function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param data input data
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' 
get_time_test <- function(formula, time, data, verbose = T){
  #extract time from formula
  modelterms <- terms(formula, specials = "time", keep.order = T)
  
  model_time <-  all.vars(formula)[unlist(attr(modelterms, "specials"))]
  
  #If both formula and argument specifiy time, check them out. 
  if(!is.null(time) & !identical(model_time, character(0))){
    stopifnot("The time argument was stated in the formula and differs from the one given as the time argument" = identical(time, model_time))
  }
  
  #Preferrably use time from the formula. 
  if(!identical(model_time, character(0))){
    time <- model_time
  }
  #Check whether time argument exists in dataset
  stopifnot("The time argument argument was not found. Use the `time = <here>` argument or assign within the formula: `y ~ x + time(<here>)" = length(time) == 1)
  stopifnot("The time argument does not match any column in the data set." = time %in% colnames(data))
  
  return(time)
}

#' Get sine and cosine components
#' @description Based on cosinor and limorhyde packages
#' @param data input data
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' 
get_cos_sine_test <- function(data, period, colnamePrefix = NULL){
  b = cbind(cos(data/period * 2 * pi), 
            sin(data/period * 2 * pi))
  colnames(b) = c("cos", "sin")
  colnames(b) = paste0(colnamePrefix, colnames(b))
  return(b)
}

#' Fit cosinor model
#' @description Fit cosinor model for totality of data
#' @param formula A formula. Use the time() function to designate which variable represents time. 
#' @param data input data
#' @param time A string. Should be the column name containing the time values.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param for_pw A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' 
fit_cosinor_model_test <- function(formula, data, time = NULL, verbose = T, for_pw = F){
  
  #Standardize formula
  formula = update.formula(formula, formula(paste("~ . + time(", time, ")")))
  
  formula = update.formula(formula, formula(paste("~ . -time(", time, ")")))
  
  #Add sine and cosine components
  formula = update.formula(formula, formula(paste0("~ . * (", time, "_cos + ", time, "_sin)" )))
  
  #fix intercept for pf_fits
  if(!for_pw){
  formula = update.formula(formula, formula(paste0("~ . -1")))
    
  }
  #in case of singular group

  if(verbose){print(paste0("Using the following model: ", c(formula)))}
  
  fit = lm(formula, data = data)
  
  return(fit)
}

#' Give tracing information for plotting purposes
#' @description Generate data needed to plot cosinor trace line. 
#' @param fit A model fit. Can be found in KronosOut@fit
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' 
kronos_predict_test = function(fit, period, time, verbose = T){
  new_data <-  fit$xlevels
  
  new_data$zt     = seq(0, period, 0.25) 
  
  new_data = expand.grid(new_data)
  
  new_data$zt_cos = cos(new_data$zt * 2 * pi / period)
  new_data$zt_sin = sin(new_data$zt * 2 * pi / period)
  
  names(new_data)[names(new_data) == "zt"]     <- time
  names(new_data)[names(new_data) == "zt_cos"] <- paste0(time, "_cos")
  names(new_data)[names(new_data) == "zt_sin"] <- paste0(time, "_sin")
  return(cbind(new_data, y_hat = predict(fit, newdata = new_data)))
}



#' Fit cosinor model
#' @description Fit cosinor model for one aspect of data
fit_groupwise_model_test <- function(data, group, time = time, period = period, verbose = T){
  
  data = data[data[,"unique_name"] == group,]
  
  form = as.formula(paste0(colnames(data)[1], " ~ (", time, "_cos + ", time, "_sin)" ))
  
  if(verbose){print(paste0("Using the following model: ", c(form) ))}
  
  fit = lm(form, data = data)
  
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
  
  
  avg       = mean(vals$pred_value)
  #highpoint = which(vals$pred_value == max(vals$pred_value))
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
  
  colnames(out_df)[1] <- "unique_group"
  
  return(out_df)
}

#' Fit pairwise cosinor models as some sort of TukeyHSD.
#' @description Fit cosinor model for subset of data.
#'
pairwise_cosinor_model_test <- function(data = fit$model, formula = formula, time = time, pairwise_t = pairwise_t, verbose = verbose){
  stopifnot("You need more that two groups in order to sensibly do pairwise comparisons. " = (unique(data[,"unique_name"])) > 2)
  combos <- combn(c(unique(data[,"unique_name"])), m = 2)
  
  for(c in 1:ncol(combos)){
    group_pair = combos[,c]
    x_pair     = data[data[,"unique_name"] %in% group_pair,]
    
    pairwise_t[[c]] <- anova(fit_cosinor_model_test(formula = update.formula(formula, ~ unique_name), data = x_pair, time = time, verbose = T, for_pw = T))
    names(pairwise_t)[c] <- paste(group_pair, collapse = " vs ")
  }
  return(pairwise_t)
}

#' Calculate standard error
#' @description calculate standard error for plotting purposes
#' @param x input vector
#' @return standard error of the input vector
#' @export
std_test <- function(x, na.rm = TRUE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x) / length(x))
}
