#' @title ME
#' @name Matrix effect
#'
#' @description The function analyzes the matrix effect present from the sample comparison
#' obtained with the analysis done in Blank sample
#'
#' @param x A vector of the independent variable
#' @param y A vector of the dependent variable Y, with the sample response being fortified
#' @param z A vector of the dependent variable Z, being blank's response
#'
#' @return Assessment regarding the matrix effect in percentage
#'
#' @author Igor Samuel Mendes
#'
#' @examples
#' #test with TCDF by DLLME. Data obtained from (SICUPIRA,2019)
#'x<- c(5.3,5.3,5.3,13.3,13.3,13.3,21.3,21.3,21.3,29.3,29.3,29.3,37.3,37.3,37.3,45.3,45.3,45.3)
#'y<- c(423832,267896,370213,1545854,1093452,780549,2378730,2204902,2709421,3492212,3591964,3523445,4762476,4830668,4853425,5334327,5753813,5648307)
#'z<- c(733485,755427,760125,1827101,1860320,1853241,2904764,2913248,3025938,4418976,4492318,4516138,5722249,5661746,5535078,6741298,7079091,7005771)
#'ME(x,y,z)
#'
#'
#'@export
#Packge Chemometrics_UFVJM
#Algorithm for Matrix effect
###################################################################
ME<- function (x,y,z){
  library(car)
  regmat = lm(y~x)
  regsol = lm(z~x)

  b0mat<-coefficients(regmat)[2]
  b0sol<-coefficients(regsol)[2]

  Res <- (b0mat/b0sol)*100

  if (Res < 100) {
    cat("Result:",Res,"% Therefore, there is a reduction in the chromatographic response\n------------------------------------------------------------------------\n")
  }
  if (Res > 100) {
    cat("Result:",Res,"% Therefore, there is an increase in the chromatographic response\n------------------------------------------------------------------------\n")
  }
  if (Res == 100) {
    cat("Result:",Res,"% Therefore, there is no matrix effect\n------------------------------------------------------------------------\n")
  }

  require(ggplot2)
  dados<-data.frame(x,y,z)
  ggplot(dados,aes(x=x,y=y))+
    geom_line(aes(col="Fortified sample"))+
    geom_line(aes(y=z,col="blank sample"))
}
