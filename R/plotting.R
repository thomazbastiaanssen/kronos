#' A plotting method for circadian plots using ggplot. 
#' @description Wrapper around ggplot to make figures with a sinusoid trace. 
#' @importFrom ggplot2 ggplot
#' @export
#' 
gg_kronos <- function(kronosOut, fill){
  
  d = merge(kronosOut@input, kronosOut@to_plot, by="row.names", all=TRUE)[,-1]
  x_obs  = paste0(kronosOut@plot_info$time, ".x")
  x_pred = paste0(kronosOut@plot_info$time, ".y")
  y_obs  = colnames(kronosOut@fit$model)[1]
  y_pred = "y_hat"
  fill_obs  = paste0(fill, ".x")
  fill_pred = paste0(fill, ".y")

  ggplot(d) +
    
    geom_rect(data=NULL, aes(xmin=12, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightgray") +
    
    geom_line(aes_string(x = x_pred, y = y_pred, colour = fill_pred)) +
    
    stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.75), width = 1, 
                 aes_string(x = x_obs, y = y_obs, group = fill_obs)) + 
    
    stat_summary(fun = mean, geom="point", position = position_dodge(0.75), size = 4, shape = 21, 
                 aes_string(x = x_obs, y = y_obs, fill = fill_obs)) +
    
    geom_point(aes_string(x = x_obs, y = y_obs, fill = fill_obs), shape = 21, position = position_dodge(0.75)) +
    
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24), expand = c(0, 0)) +
    
    theme_bw() +
    
    guides(fill = guide_legend(title=fill), colour = guide_none()) +
    
    xlab("Time (h)") +
    
    ylab("Gene expression")
}