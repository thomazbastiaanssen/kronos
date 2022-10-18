#' A plotting method for circadian plots using ggplot. 
#' @description Wrapper around ggplot to make figures with a sinusoid trace. 
#' @importFrom ggplot2 ggplot
#' @export
#' 
gg_kronos_sinusoid <- function(kronosOut, fill = NULL){
  
  d = merge(kronosOut@input, kronosOut@to_plot, by="row.names", all=TRUE)[,-1]
  x_obs  = paste0(kronosOut@plot_info$time, ".x")
  x_pred = paste0(kronosOut@plot_info$time, ".y")
  y_obs  = colnames(kronosOut@fit$model)[1]
  y_pred = "y_hat"
  fill_obs  = paste0(fill, ".x")
  fill_pred = paste0(fill, ".y")
  period    = kronosOut@plot_info$period
  
  if(is.null(fill)){
    d$.x <- ""
    d$.y <- ""
  }

  ggplot(d) +
    
    geom_rect(data=NULL, aes(xmin=12, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightgray") +
    
    geom_line(aes_string(x = x_pred, y = y_pred, colour = fill_pred)) +
    
    stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.75), width = 1, 
                 aes_string(x = x_obs, y = y_obs, group = fill_obs)) + 
    
    stat_summary(fun = mean, geom="point", position = position_dodge(0.75), size = 4, shape = 21, 
                 aes_string(x = x_obs, y = y_obs, fill = fill_obs)) +
    
    geom_point(aes_string(x = x_obs, y = y_obs, fill = fill_obs), shape = 21, position = position_dodge(0.75)) +
    
    scale_x_continuous(breaks = seq(0, period, period/4), expand = c(0, 0)) +
    
    theme_bw() +
    
    guides(fill = guide_legend(title=fill), colour = guide_none()) +
    
    xlab("Time (h)") +
    
    ylab("Gene expression")
}

#' A plotting method for circadian plots using ggplot. 
#' @description Wrapper around ggplot to make figures with a sinusoid trace. 
#' @importFrom ggplot2 ggplot
#' @export
#' 
gg_kronos_circle <- function(kronosOut){
  
  d = kronosOut@ind_fit
  d$amp_scale = d$amplitude / max(d$amplitude)
  period = kronosOut@plot_info$period
  
  ggplot(d) +
    aes(x = acro, y = amp_scale) +
    
    geom_segment(aes(x = acro, xend = acro, y = 0, yend = amp_scale, linetype = p.val < 0.05)) +
    
    geom_point(shape = 21, size = 3, aes(fill = unique_group)) +

    coord_polar(theta = "x") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
    scale_x_continuous(breaks = seq(0, period, period/6), limits = c(0, 24)) +
    scale_y_continuous(breaks = seq(0, 1, 1/3), limits = c(0, 1)) +
    
    theme_bw()  +
    
    guides(fill = guide_legend("Legend"), colour = guide_none(), linetype = guide_none()) +
    
    theme(#panel.grid.minor = element_blank(), 
          axis.text = element_text(size = 14, colour = "black"), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid.major.y = element_line(colour = "darkgray"), 
          panel.border = element_blank() ) +
    
    xlab("") +
    
    ylab("")
}
