#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @param kronosOut a \code{kronosOut} output object. 
#' @param target the specific entry of the \code{kronosOut} object to be extracted. 
#'
#' @return A \code{data.frame} of results.
#' @importFrom methods slot slotNames
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
#' 
#' @export
getKronos <- function(kronosOut, target){
  stopifnot("The target entry could not be found in the kronosOut object." = target %in% slotNames(kronosOut))
  
  return(slot(object = kronosOut, name = target))
  }

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The data used as input for the model.
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
#' 
#' @export
getKronos_input <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "input"))
}

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The model fit used.
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
getKronos_fit <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "fit"))
  }

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The traces per group for plotting.
#' @examples
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
#' @export
getKronos_trace <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "to_plot"))
  }

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return Rhythmicity parameters per group.
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
#' 
#' @export
getKronos_groupwise <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "ind_fit"))
  }

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return Pairwise comparisons between groups.
#' @examples
#' 
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
#' 
#' @export
getKronos_pairwise <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "pairwise_models"))
  }

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return ANOVA-like adjusted p-values for how each factor interacts with time.
#' @examples
#' 
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
#' 
#' @export
getKronos_pairwise_p <- function(kronosOut){
  pairwise_p <- data.frame(getKronos(kronosOut = kronosOut, target = "pairwise_p_vals"))
  
  #Fix colname
  colnames(pairwise_p) <- "adj.p.val"
  return(pairwise_p)
}

#' Get Results from \code{KronosOut} Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The names and values of additional circadian model parameters, mostly for plotting purposes.
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
#' 
#' @export
getKronos_params <- function(kronosOut){
  
  return(data.frame(getKronos(kronosOut = kronosOut, target = "plot_info")))
  }

#' Get Results from list of \code{KronosOut} Objects
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a list of \code{kronosOut} objects.
#' @param kronos_list a list of preferrably named \code{kronosOut} objects.
#' @param padjust a boolean. Toggles FDR using Benjamini Hochbergs procedure.
#'
#' @return ANOVA-like adjusted p-values for how each factor interacts with time.
#' 
delistKronos_pairwise <- function(kronos_list, padjust = TRUE){
  res = do.call(rbind,(lapply(kronos_list, FUN = function(x){t(getKronos_pairwise_p(x))})))

  row.names(res) <- NULL
  if(!is.null(names(kronos_list))){
    row.names(res) <- names(kronos_list)
  }
  colnames(res) = paste0(colnames(res), ":Time_p.val")
  
  if(padjust){
    resadj = apply(res, 
                   MARGIN = 2, 
                   p.adjust, method = "BH", 
                   simplify = TRUE)
    colnames(resadj) = gsub(x = colnames(resadj),pattern = "p\\.val", replacement = "q\\.val")
    
    res = cbind(res, resadj)
  }
  
  return(res)
}

#' Get Results from list of KronosOut Objects
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a list of \code{kronosOut} objects.
#' @param kronos_list a list of preferrably named kronosOut objects.
#' @param padjust a boolean. Toggles FDR using Benjamini Hochbergs procedure.
#'
#' @return A table with circadian output stats per group per feature.
#' @importFrom stats p.adjust
#' 
delistKronos_groupwise <- function(kronos_list, padjust = TRUE){
  res = do.call(rbind,lapply(kronos_list, FUN = function(x){
    #Get table per feature
    out_table = getKronos_groupwise(x)
    #Vectorize results
    out_vec = unlist(out_table[,-1])
    #Fix colnames
    names(out_vec) = (paste(out_table[,1],rep(colnames(out_table[,-1]),each = nrow(out_table) ),sep = "_"))
    return(out_vec)}))
  
  row.names(res) <- NULL
  if(!is.null(names(kronos_list))){
    row.names(res) <- names(kronos_list)
  }
  
  if(padjust){
    resadj = apply(res[,grepl(pattern = "p\\.val", x = colnames(res))], 
                   MARGIN = 2, 
                   p.adjust, method = "BH", 
                   simplify = TRUE)
    colnames(resadj) = gsub(x = colnames(resadj),pattern = "p\\.val", replacement = "q\\.val")
    
    res = cbind(res, resadj)
  }
  return(res)
}

#' Wrangle results from list of \code{KronosOut} Objects to publication ready table. 
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a list of \code{kronosOut} objects.
#' @param kronos_list a list of preferrably named \code{kronosOut} objects.
#' @param padjust a boolean. Toggles FDR using Benjamini and Hochbergs procedure.
#'
#' @return A table with circadian output stats per group per feature.
#' @importFrom stats p.adjust
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
#' #Plotting:
#' gg_kronos_acrogram(out_list)
#' 
#' @export
kronosListToTable <- function(kronos_list, padjust = TRUE){
  cbind(delistKronos_groupwise(kronos_list = kronos_list, padjust = padjust),
        delistKronos_pairwise(kronos_list  = kronos_list, padjust = padjust))
} 
