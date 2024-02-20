#' @title LDLQ
#' @name Limit of detection and quantification
#'
#' @description The function informs the analyte's Limit of Detection and its Limit of
#' Quantification
#'
#' @param x Value of the fortified area
#' @param y white area value
#'
#' @return Values of LD e LQ
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#'x = 49.480
#'y = 17.054
#'LDLQ(x,y)
#'
#'@export
#Function LDLQ simplified
LDLQ<- function(x,y){
  LD = (x/y) * 3
  LQ = (x/y) * 10
  cat("Detection limit:",LD)
  cat("\nLimit of quantification:",LQ)
}

