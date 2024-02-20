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
  yt<- ym*y
  zt<- ym*z

  require(ggplot2)
  dados<-data.frame(x,yr,zr)
  ggplot(dados,aes(x=x,y=yr))+
    geom_line(aes(col="Fortified sample"))+
    geom_line(aes(y=zr,col="blank sample"))
}
