#' @title Sel
#' @name Selectivity
#'
#' @description The function provides the graph of the fortified sample curve, as well as
#' that of your white sample, for visual comparison purposes regarding selectivity
#' of analysis
#'
#' @param tdp Total data points
#' @param t Time,in minutes, of the Chromatographic run
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
SEL<- function(tdp,t,ym,y,z){
  x <- c(1:tdp)
  xm<- t/tdp
  xt <- numeric(length(x))

  for (i in 1:(length(x))) {
    xt[i] <- (x[i]-1) * xm
  }

  yr<- y[1:tdp]
  zr<- z[1:tdp]
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
