#' @title AP
#' @name Accuracy & Precision
#'
#' @description The function analyzes the accuracy and precision of the data provided
#'
#' @param x The value of the analyzed concentration
#' @param y A vector of the dependent variable Y
#' @param z The average of y
#'
#' @return Assessment regarding accuracy and precision
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#'Test with TCDF MELDD methodology. Data obtained from (SICUPIRA,2019)
#'x<- 0.02
#'y<- c(233153.0,257269.0,307816.0,256882.0,291418.0,264616.0,262520.0)
#'AP(x,y,z)
#'@export
#Packge Chemometrics_UFVJM
#Algorithm for  check accuracy and precision
AP<- function(x,y){
  dpr<- (sd(y)/median(y)*100)
  z= mean(y)
  rec<- (median(y)/z) *100
  cat("Result:",rec,"+-",dpr)
  cat("\n\nInterpreting the result:\nThe method is considered to be accurate
if the average recovery rates of the analyte are within
from the range of 70 to 120%. On the other hand, the method will
be considered accurate if the standard deviations
relatives are less than 20%.")
}

