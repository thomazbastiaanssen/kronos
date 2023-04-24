#' An S4 class to contain all output parameters and plotting information from the 'kronos' main pipeline
#'
#' @slot input A data.frame with the data that was fed to the main workhorse function as 'x'
#' @slot fit An lm fit for the entire model for the purpose of assessing differential rhytmicity. 
#' @slot to_plot A data.frame with the traces required to plot individual sinusoid curves
#' @slot ind_fit A data.frame with the parameters from individual rhythmic model fits. 
#' @slot pairwise_t A data.frame with the p.values for pairwise comparisons, if applicable. 
#' @description kronosOut is the main output container of the main 'kronos' functions. 
#'
setClass("kronosOut",
         slots = c(
           input           = "data.frame",
           fit             = "lm",
           to_plot         = "data.frame", 
           ind_fit         = "data.frame", 
           pairwise_models = "list", 
           pairwise_p_vals = "vector", 
           plot_info       = "list"
         )
)

#' Show method for \code{kronosOut} object.
#' @description method to print \code{kronosOut} object by calling \code{show}.
#' Since kronosOut objects are typically unwieldy, also gives some tips on how to handle it.
#' @param object An \code{kronosOut} object.
#' @importFrom utils str
#' @importFrom methods show
#' @return Does not return anything, for efficiency reasons. The obvious side effect is output to the terminal.
#' @export
#'
setMethod("show", "kronosOut", function(object){
  str(object)
  cat("\nThis is an kronosOut S4 object. They tend to be large so here's a summary instead.
      You can extract information from it using the getKronos() functions,
      or you can manually explore it by using the @ operator. ")
}
)
