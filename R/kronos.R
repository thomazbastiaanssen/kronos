#' Fit a cosinor model and extract relevant parameters
#' @description This is the main workhorse function in the 'kronos' package. It manages the individual functionalities of 'kronos', including rhythmicity analysis and differential rhythmicity.
#' @param data input data
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} variable  
#' @param pairwise A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#' @return A \code{kronosOut} S4 object containing coefficients and all operations.
#' @examples
#' #Load prepared data stored in Kronos library
#' data("kronos_demo")
#' output <- kronos(formula = Variable_1 ~ time(Timepoint), 
#' data = onevariable, period = 24, verbose = TRUE, pairwise = FALSE)
#' 
#' #Extracting data from the output object:
#' getKronos_fit(output)
#' getKronos_trace(output)
#' getKronos_groupwise(output)
#' 
#' #Plotting:
#' gg_kronos_circle(output)
#' gg_kronos_sinusoid(output)
#' 
#' #For high-dimensional data, use fw_kronos:
#' out_list = fw_kronos(x = bigdata[1:50,], formula = ~ Group + time(Timepoint), 
#' metadata = bigmeta, period = 24, verbose = FALSE, pairwise = TRUE) 
#'
#' #Extracting data from the output object: 
#' kronosListToTable(out_list)
#' 
#' 
#' #Plotting:
#' gg_kronos_acrogram(out_list)
#' 
#' @importFrom methods new
#' @importFrom stats anova as.formula coef lm pf predict terms update.formula get_all_vars
#' @importFrom utils combn
#' @export
kronos <- function(formula, data, time = NULL, period = 24, verbose = TRUE, pairwise = TRUE){

if(!is.null(time)){
  formula = as.formula(paste0(Reduce(paste, deparse(formula)), " + time(", time, ")"))
}    
#Clean and standardize input data
data <- get_all_vars(formula = formula, data = data)
  
#Identify time component
var_out  <- get_vars(formula = formula, time = time, data = data, verbose = verbose)  
time     <- var_out$time
response <- var_out$response
factors  <- var_out$factors

#convert all non-time and non-response variables to categorical data. 
data[,{{factors}}] <- sapply(data[,{{factors}}], "as.character")

if(is.null(factors)){data$unique_group <- TRUE}

if(!is.null(factors) & length(factors) == 1){
  data$unique_group = data[,{{factors}}]
}

if(!is.null(factors) & length(factors) > 1){
data$unique_group = do.call(what = "paste", c(data[,{{factors}}], sep = "_"))
}
#Set up sine and cosine component
data <- cbind(data, get_cos_sine(data = data[,time], colnamePrefix = paste0(time, "_"), period = period))

#Update formula to reflect sine and cosine component
formula <- build_kronos_formula(formula = formula, time = time, verbose = verbose)

#Fit the general model
fit  <- fit_cosinor_model(formula = formula, data = data, time = time, verbose = verbose, for_pw = FALSE)

#Predict values for plotting
vals <- kronos_predict(fit = fit, period = period, time = time, factors = factors, verbose = verbose)

#In case of no non-time factors
if(length(fit$xlevels) == 0){
  fit$model$unique_group = TRUE
  groupwise = vector(mode = "list", length = 1)
  groupwise[[1]] <- fit_groupwise_model(data = fit$model, group = TRUE, 
                                             time = time, period = period, verbose = verbose)
}

#In case of multiple groups
if(length(fit$xlevels) == 1){
#fit$model$unique_group <- paste(fit$model[,names(fit$xlevels)])

fit$model$unique_group <- do.call(what = "paste", list(data[,{{factors}}], sep = "_"))[!is.na(data[,response])]

groupwise = vector(mode = "list", length = length(unique(fit$model$unique_group)))

for(g in 1:length(groupwise)){
  groupwise[[g]] = fit_groupwise_model(data = fit$model, 
                                       group = unique(fit$model$unique_group)[g], 
                                       time = time, period = period, verbose = verbose)
  }
}

#In case of multiple groups
if(length(fit$xlevels) > 1){
  #fit$model$unique_group <- paste(fit$model[,names(fit$xlevels)])
  
  fit$model$unique_group <- do.call(what = "paste", c(data[,{{factors}}], sep = "_"))[!is.na(data[,response])]
  
  groupwise = vector(mode = "list", length = length(unique(fit$model$unique_group)))
  
  for(g in 1:length(groupwise)){
    groupwise[[g]] = fit_groupwise_model(data = fit$model, 
                                         group = unique(fit$model$unique_group)[g], 
                                         time = time, period = period, verbose = verbose)
  }
}


bygroup = do.call(rbind, groupwise)
row.names(bygroup) <- NULL

pairwise_models <- list("pairwise was not activated.")
if(pairwise){  
  pairwise_models <- list()
  if(verbose){print("Fitting pairwise models")}
  pairwise_models = pairwise_cosinor_model(data = fit$model, formula = formula, time = time, verbose = verbose)
}

pairwise_p_vals <- "pairwise was not activated."
if(pairwise){  
  pairwise_p_vals <- kronos_anova(fit = fit, time = time)
}

#initialize output container object
output = new("kronosOut", 
             input           = data,
             fit             = fit,
             to_plot         = vals,
             ind_fit         = bygroup, 
             pairwise_models = pairwise_models, 
             pairwise_p_vals = pairwise_p_vals, 
             plot_info       = list(time = time, response = response, period = period))

return(output)
}

#' Figure out what variable represents time. Called by main 'kronos' function. 
#' @description Extracts time from the formula and from the time argument. Also handles inconsistencies. 
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param data input data
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#' 
get_vars <- function(formula, time, data, verbose = verbose){
  #extract time from formula
  modelterms <- terms(formula, specials = "time", keep.order = TRUE)
  
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
  
  model_response <-  all.vars(formula)[unlist(attr(modelterms, "response"))]
  model_factors  <-  all.vars(formula)[-c(unlist(attr(modelterms, "response")), 
                                          unlist(attr(modelterms, "specials")))]
  if(identical(character(0), model_factors)){
    model_factors <- NULL
  }
  return(list(time = time, response = model_response, factors = model_factors))
}

#' Get sine and cosine components
#' @description Based on 'cosinor' and 'limorhyde' packages
#' @param data input data
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param colnamePrefix A character string. Typically the name of the "Time" variable. 
#' 
get_cos_sine <- function(data, period, colnamePrefix = NULL){
  b = cbind(cos(data/period * 2 * pi), 
            sin(data/period * 2 * pi))
  colnames(b) = c("cos", "sin")
  colnames(b) = paste0(colnamePrefix, colnames(b))
  return(b)
}

#' Fit cosinor model
#' @description Fit cosinor model for totality of data
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. 
#' @param data Input data
#' @param time A string. Should be the column name containing the time values.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @param for_pw A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' 
fit_cosinor_model <- function(formula, data, time = NULL, verbose = verbose, for_pw = FALSE){
  #fix intercept for pairwise fits
  
  if(!for_pw & length(all.vars(formula)) > 3){
  formula = update.formula(formula, formula(paste0("~ . -1")))
  }
  
  #in case of singular group
  if(verbose){print(paste0("Using the following model: ", Reduce(paste0, deparse(formula))))}
  
  fit = lm(formula, data = data)
  
  return(fit)
}

#' Give tracing information for plotting purposes
#' @description Generate data needed to plot cosinor trace line. 
#' @param fit A model fit
#' @param time A string. Should be the column name containing the time values.  
#' @param factors A vector. The names of the independent variables. 
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#' 
kronos_predict = function(fit, period, time, factors, verbose = verbose){
  new_data <-  fit$xlevels
  
  new_data$zt     = seq(0, period, 0.25) 
  
  new_data = expand.grid(new_data)
  
  new_data$zt_cos = cos(new_data$zt * 2 * pi / period)
  new_data$zt_sin = sin(new_data$zt * 2 * pi / period)
  
  names(new_data)[names(new_data) == "zt"]     <- time
  names(new_data)[names(new_data) == "zt_cos"] <- paste0(time, "_cos")
  names(new_data)[names(new_data) == "zt_sin"] <- paste0(time, "_sin")
  new_data = cbind(new_data, y_hat = predict(fit, newdata = new_data))
  
  if(is.null(factors)){new_data$unique_group <- TRUE}
  
  if(!is.null(factors) & length(factors) == 1){
    new_data$unique_group = new_data[,{{factors}}]
  }
  
  if(!is.null(factors) & length(factors) > 1){
    new_data$unique_group = do.call(what = "paste", c(new_data[,{{factors}}], sep = "_"))
  }
  
  return(new_data)
}



#' Fit cosinor model
#' @description Fit cosinor model for one aspect of data. Called by main 'kronos' function. 
#' @param data input data
#' @param group A character string. Signifies which group will be assessed. 
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#' 
fit_groupwise_model <- function(data, group, time, period, verbose){
  
  data = data[data[,"unique_group"] == group,]
  
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
#' @description Fit cosinor model for subset of data. Called by main 'kronos' function. 
#' @param data input data
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#'
pairwise_cosinor_model <- function(data, formula, time, verbose){
  stopifnot("You need more that two groups in order to sensibly do pairwise comparisons. " = length(unique(data[,"unique_group"])) > 1)
  combos <- combn(c(unique(data[,"unique_group"])), m = 2)
  
  pairwise_models = list()
  
  for(c in 1:ncol(combos)){
    group_pair = combos[,c]
    x_pair     = data[data[,"unique_group"] %in% group_pair,]
    formula_pair = update.formula(formula, ~ unique_group)
    formula_pair = build_kronos_formula(formula = formula_pair, time = time)
    pairwise_models[[c]] <- anova(fit_cosinor_model(formula = formula_pair, data = x_pair, time = time, verbose = verbose, for_pw = TRUE))
    names(pairwise_models)[c] <- paste(group_pair, collapse = " vs ")
  }
  return(pairwise_models)
}

#' Update 'kronos' formula in light of sine and cosine components
#' @description Update 'kronos' formula.
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. 
#' @param time A string. Should be the column name containing the time values.  
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large data sets.
#'
build_kronos_formula <- function(formula, time, verbose){
  #Trim off the time component
  temp_form = gsub("time\\(.*\\)", Reduce(paste0, deparse(formula)), replacement = "")
  temp_form = gsub("[ ]*\\*[ ]*$|[ ]*\\+[ ]*$", replacement = "", temp_form)
  
  #trim off trailing spaces
  temp_form = gsub("1$",    replacement = "", temp_form)
  temp_form = gsub("[ ]*$", replacement = "", temp_form)
  
  #if there were no other components to the formula:
  if(grepl("\\~$", x = temp_form)){
    return(formula(paste0(temp_form, " ", time, "_cos + ", time, "_sin")))
  }else{
  return(formula(paste0(temp_form, " * (", time, "_cos + ", time, "_sin)")))
    }
}

#' Extract p-value from full fit
#' @description Compute p-values from full fit. 
#' @param fit A lm model fit.  
#' @param time A string. Should be the column name containing the time values.  
#'
kronos_anova <- function(fit, time){
  anova.fit    <- anova(fit)
  time_int     <- grep(row.names(anova.fit), pattern = paste0(":", time, "_sin", "|", ":", time, "_cos"))
  pvals        <- anova.fit$`Pr(>F)`[time_int]
  pvals        <- pmin(pmin(pvals[1:length(pvals)%%2 ==0], pvals[1:length(pvals)%%2 ==1])*2, 1)
  names(pvals) <- unique(gsub(row.names(anova.fit), pattern = paste0(":", time, "_sin", "|", ":", time, "_cos"), replacement =  "")[time_int])
  return(pvals)
}


#' Fit a cosinor model and extract relevant parameters on a feature table.
#' @description This wrapper applies kronos(), the main workhorse function in the 'kronos' package. It manages the individual functionalities of 'kronos', including rhythmicity analysis and differential rhythmicity.
#' @param x Input data. A table with rows being features and columns being samples
#' @param formula A formula. Use the \code{time} function to designate which variable represents time. Leave the left-hand side of the formula empty as it will be sequentially replaced by every feature in the table.
#' @param metadata A metadata table, with rows being samples and columns being metadata entries
#' @param time A string. Should be the column name containing the time values.  
#' @param period A numeric. The length of a period, in the same format as the \code{time} parameter.  
#' @param pairwise A boolean. Toggles whether to perform pairwise ANOVAs as a TukeyHSD-like post-hoc. 
#' @param verbose A boolean. Toggles whether to print diagnostic information while running. Useful for debugging errors on large datasets.
#' @return A list of \code{kronosOut} S4 objects containing coefficients and all operations.
#' @importFrom methods new
#' @importFrom stats anova as.formula coef lm pf predict terms update.formula get_all_vars
#' @importFrom utils combn
#' @examples
#' #Load prepared data stored in 'kronos' library
#' data("kronos_demo")
#' output <- kronos(formula = Variable_1 ~ time(Timepoint), 
#' data = onevariable, period = 24, verbose = TRUE, pairwise = FALSE)
#' 
#' #Extracting data from the output object:
#' getKronos_fit(output)
#' getKronos_trace(output)
#' getKronos_groupwise(output)
#' 
#' #Plotting:
#' gg_kronos_circle(output)
#' gg_kronos_sinusoid(output)
#' 
#' #For high-dimensional data, use fw_kronos:
#' out_list = fw_kronos(x = bigdata[1:50,], formula = ~ Group + time(Timepoint), 
#' metadata = bigmeta, period = 24, verbose = FALSE, pairwise = TRUE) 
#'
#' #Extracting data from the output object: 
#' kronosListToTable(out_list)
#' 
#' 
#' #Plotting:
#' gg_kronos_acrogram(out_list)
#' @export
fw_kronos <- function(x, formula, metadata, time = NULL, period = 24, verbose = FALSE, pairwise = FALSE){
  formula = update.formula(x_feature ~ ., formula)
  apply(X = x, MARGIN = 1, FUN = function(y){
    kronos(formula = formula, 
           data = cbind("x_feature" = y, metadata),
           time = time, period = period,  
           verbose = verbose, 
           pairwise = pairwise)})
  
}
  