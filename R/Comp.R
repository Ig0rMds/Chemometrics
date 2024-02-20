#' @title COMP
#' @name Comparison
#'
#' @description The function provides the graph of the fortified sample curve, as well as
#' that of your white sample, for visual comparison purposes regarding selectivity
#' of analysis
#'
#' @param x is a data frame
#'
#' @return Test of comparations
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#' # Entrace of datas
#'DT = c(25, 17, 27, 21, 15,
#'           10, -2, 12, 4, 16,
#'          18, 8, 4, 14, 6,
#'          23, 29, 25, 35,33,
#'           11, 23, 5, 17, 9,
#'           8, -6, 6, 0, 2)
#'
#'# Data processing
#'(dats = data.frame(trat = factor(rep(c('A','B','C','D','E','Controle'),
#'                                      each=5)), DT))
#'attach(dats)
#'COMP (dats)
#'@export
#Package chemometrics_UFVJM
#Algorithm for multiples comparasions

COMP<- function(x){

require(easyanova)
ea1(x, design = 1)
}

