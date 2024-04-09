#' @title Sel
#' @name Selectivity
#'
#' @description The function provides the graph of the fortified sample curve, as well as
#' that of your white sample, for visual comparison purposes regarding selectivity
#' of analysis
#'
#' @param t Total data points
#' @param xm X axys multiplier
#' @param ym Y axys multiplier
#' @param y spiked analyte response
#' @param z blank sample response
#'
#' @return Graph for selectivity analysis
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#' It's difficult to explain how to use it here, check out my thesis to understand the
#' implementation
#'
#'@export
#Package chemometrics_UFVJM
#Algorithm for selectivity
Sel<- function(t,xm,ym,y,z){
  x <- c(1: t-1)
  xt<- xm*x
  yr<- y[1:t]
  zr<- z[1:t]
  yt<- ym*yr
  zt<- ym*zr

  require(ggplot2)
  dados<-data.frame(x,yr,zr)

  # Gráfico 1: xt vs yt
  plot_xt_yt <- ggplot(dados, aes(x = xt)) +
    geom_line(aes(y = yt, colour = "Fortified sample")) +
    scale_color_manual(values = c("Fortified sample" = "black"))+
    labs(x = "xt", y = "yt")
  # Gráfico 2: xt vs zt
  plot_xt_zt <- ggplot(dados, aes(x = xt)) +
    geom_line(aes(y = zt, color = "blank sample")) +
    scale_color_manual(values = c("blank sample" = "blue"))+
    labs(x = "xt", y = "zt")

  # Mostrar os dois gráficos na interface "plots"
  require("gridExtra")
  show_plots <- function(plots) {
    gridExtra::grid.arrange(grobs = plots, ncol = 1)
  }

  # Exibir os gráficos
  show_plots(list(plot_xt_yt, plot_xt_zt))
}
