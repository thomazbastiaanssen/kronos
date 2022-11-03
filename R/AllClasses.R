#' An S4 class to contain all output parameters and plotting information from the kronos main pipeline
#'
#' @slot input A data.frame with the data that was fed to the main workhorse function as 'x'
#' @slot fit An lm fit for the entire model for the purpose of assessing differential rhytmicity. 
#' @slot to_plot A data.frame with the traces required to plot individual sinusoid curves
#' @slot ind_fit A data.frame with the parameters from individual rhythmic model fits. 
#' @slot pairwise_t A data.frame with the p.values for pairwise comparisons, if applicable. 
#' @description kronosOut is the main output container of the kronos function. 
#'
setClass("kronosOut",
         slots = c(
           input      = "data.frame",
           fit        = "lm",
           to_plot    = "data.frame", 
           ind_fit    = "data.frame", 
           pairwise_t = "list", 
           pairwise_p = "vector", 
           plot_info  = "list"
         )
)

