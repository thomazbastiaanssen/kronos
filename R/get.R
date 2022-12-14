#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @param kronosOut a \code{kronosOut} output object. 
#' @param target the specific entry of the \code{kronosOut} object to be extracted. 
#'
#' @return A \code{data.frame} of results.
#' @importFrom methods slot slotNames
#' @export
getKronos <- function(kronosOut, target){
  stopifnot("The target entry could not be found in the kronosOut object." = target %in% slotNames(kronosOut))
  
  return(slot(object = kronosOut, name = target))
  }

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The data used as input for the model.
#'
#' @export
getKronos_input <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "input"))
}

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The model fit used.
#'
#' @export
getKronos_fit <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "fit"))
  }

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The traces per group for plotting.
#'
#' @export
getKronos_trace <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "to_plot"))
  }

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return Rhythmicity parameters per group.
#'
#' @export
getKronos_groupwise <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "ind_fit"))
  }

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return Pairwise comparisons between groups.
#'
#' @export
getKronos_pairwise <- function(kronosOut){
  
  return(getKronos(kronosOut = kronosOut, target = "pairwise_models"))
  }

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return ANOVA-like adjusted p-values for how each factor interacts with time.
#'
#' @export
getKronos_pairwise_p <- function(kronosOut){
  
  return(data.frame(getKronos(kronosOut = kronosOut, target = "pairwise_p_vals")))
}

#' Get Results from KronosOut Object
#'
#' These functions provides a unified wrapper to retrieve results
#'  from a \code{kronosOut} object.
#' @inheritParams getKronos
#'
#' @return The names and values of additional circadian model parameters, mostly for plotting purposes.
#'
#' @export
getKronos_params <- function(kronosOut){
  
  return(data.frame(getKronos(kronosOut = kronosOut, target = "plot_info")))
  }

