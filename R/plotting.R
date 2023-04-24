#' A plotting method for circadian plots using 'ggplot2'. 
#' @description Wrapper around 'ggplot2' to make figures with a sinusoid trace. 
#' @param kronosOut an output object from the main \code{kronos} function. 
#' @param fill The name of the variable that should be used to mark different groups. In the case of a single group, leave empty. 
#' @return a 'ggplot2' compatible object. 
#' @import ggplot2
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
#' @export
#' 
gg_kronos_sinusoid <- function(kronosOut, fill = "unique_group"){
  
  requireNamespace("ggplot2")
  
  d = merge(kronosOut@input, kronosOut@to_plot, by="row.names", all=TRUE)[,-1]
  x_obs  = paste0(kronosOut@plot_info$time, ".x")
  x_pred = paste0(kronosOut@plot_info$time, ".y")
  y_obs  = colnames(kronosOut@fit$model)[1]
  y_pred = "y_hat"
  fill_obs  = "unique_group.x"
  fill_pred = "unique_group.y"
  period    = kronosOut@plot_info$period
  
  if(sum(grepl(pattern = "unique_group.x|unique_group.y", x = colnames(d))) == 0){
    d$unique_group.x <- ""
    d$unique_group.y <- ""
  }
  ggplot(d) +
    
    geom_rect(data=NULL, aes(xmin=12, xmax=Inf, ymin=-Inf, ymax=Inf), fill="lightgray") +
    
    # geom_text(x = -Inf, y = Inf, hjust = -1, vjust = 2, label = "🌞", size = 8) +
    # geom_text(x = Inf,  y = Inf, hjust = 2, vjust = 2,  label = "🌙️",  size  =8 ) +
    # 
    geom_line(aes_string(x = x_pred, y = y_pred, colour = fill_pred)) +
    
    stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.75), width = 1, 
                 aes_string(x = x_obs, y = y_obs, group = fill_obs)) + 
    
    stat_summary(fun = mean, geom="point", position = position_dodge(0.75), size = 4, shape = 21, 
                 aes_string(x = x_obs, y = y_obs, fill = fill_obs)) +
    
    geom_point(aes_string(x = x_obs, y = y_obs, fill = fill_obs), shape = 21, position = position_dodge(0.75)) +
    
    scale_x_continuous(breaks = seq(0, period, period/4), expand = c(0, 0)) +
    
    theme_bw() +
    guides(fill = guide_legend("Legend", override.aes = list(size = 3)), colour = guide_none()) +
    
    xlab("") +
    
    ylab("")

}

#' A plotting method for circadian plots using 'ggplot2'. 
#' @description Wrapper around 'ggplot2' to make circadian circleplots. 
#' @param kronosOut an output object from the main \code{kronos} function. 
#' @return a 'ggplot2' compatible object. 
#' @import ggplot2
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
#' @export
#' 
gg_kronos_circle <- function(kronosOut){
  
  requireNamespace("ggplot2")
  
  d = kronosOut@ind_fit
  d$amp_scale = d$amplitude / max(d$amplitude)
  
  #Declare variables for check
  period = kronosOut@plot_info$period
  acro = "acro"
  amp_scale = "amp_scale"
  p.val = "p.val"
  unique_group = "unique_group"
  p.sig = "p.sig"
  d$p.sig = d$p.val < 0.05

    ggplot(d) +
    aes_string(x = acro, y = amp_scale) +
    
    geom_segment(aes_string(x = acro, xend = acro, y = 0, yend = amp_scale, linetype = p.sig)) +
    
    geom_point(shape = 21, size = 3, aes_string(fill = unique_group)) +

    coord_polar(theta = "x") +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
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

#' A plotting method for acrophase circleplots using 'ggplot2'. 
#' @description Wrapper around 'ggplot2' to make circadian circleplots. 
#' @param kronosOutList A list of KronosOut output objects from the main \code{kronos} function. 
#' @return a 'ggplot2' compatible object. 
#' @import ggplot2
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
#' @export
#' 
gg_kronos_acrogram <- function(kronosOutList){

d      <- do.call(rbind, lapply(kronosOutList, function(x) x@ind_fit))

period <- kronosOutList[[1]]@plot_info$period

ggplot(d) +
  aes_string(x = "acro", fill = "unique_group") +
  
  geom_histogram(position = "identity", color = 'black', alpha = 3/4, bins = period * 2, boundary = 0) + 
  
  coord_polar() +
  
  scale_x_continuous(breaks = seq(0, period, period/6), limits = c(0, period)) +

  theme_bw()+
  
  guides(fill = guide_legend("Legend"), colour = guide_none(), linetype = guide_none()) +
  
  theme(
    axis.text = element_text(size = 14, colour = "black"), 
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid.major.y = element_line(colour = "darkgray"), 
    panel.border = element_blank() ) +
  
  xlab("") +
  
  ylab("")
}
