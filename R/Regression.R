#' @title REG
#' @name Linear Regression
#'
#' @description It is a function that performs the linear regression of a model,
#' testing hypotheses that ensure that the linear model is adequate, as well as whether the
#' regression is significant.
#'
#' @param trat A vector of the independent variable
#' @param resp A vector of the dependent variable Y
#' @param degree Choice of regression degree, which can be linear (1), quadratic (2)
#' or cubic (3)
#'
#'
#' @details The function also incorporates polynomial regression, but the graphical part
#' of this regression model has not yet been done
#'
#' @return Model regression and analysis of regression parameters
#'
#' @author Igor Samuel Mendes
#'
#'
#'
#' @examples
#' #teste com TCDD pelo MELDD (Dados obtidos de SICUPIRA,2019)
#'a<- c(5.3,5.3,5.3,13.3,13.3,13.3,21.3,21.3,21.3,29.3,
#'      29.3,29.3,37.3,37.3,37.3,45.3,45.3,45.3)
#'b<-c(423832.00000,267896.00000,370213.00000,780549.00000,1545854.00000,1093452.00000,
#'     2378730.00000,2709421.00000,2204908.00000,3591964.00000,3492212.00000,
#'     3523445.00000,4762476.00000,4830668.00000,4853425.00000,5648307.00000,
#'     5334327.00000,5753813.00000)
#'REG(a,b,grau=1)
#' note: In these data, there is an outlier at point 4
#'
#'@export
#Package chemometrics_UFVJM
#Regression function
####################################################################
REG<-function (trat, resp, ylab = "Response", xlab = "Independent",
                     yname.poly = "y", xname.poly = "x", grau = NA, theme = theme_classic(),
                     point = "mean_sd", color = "gray80", posi = "top", textsize = 12,
                     se = FALSE, ylim = NA, family = "sans", pointsize = 4.5,
                     linesize = 0.8, width.bar = NA, n = NA, SSq = NA, DFres = NA)
{
  library(ggplot2)
  library(gridExtra)
  library(car)
  library(ggpubr)
  if (is.na(width.bar) == TRUE) {
    width.bar = 0.1 * mean(trat)
  }
  if (is.na(grau) == TRUE) {
    grau = 1
  }
  dados = data.frame(trat, resp)
  medias = c()
  dose = tapply(trat, trat, mean, na.rm = TRUE)
  mod = c()
  mod1 = c()
  mod2 = c()
  modm = c()
  mod1m = c()
  mod2m = c()
  text1 = c()
  text2 = c()
  text3 = c()
  mods = c()
  mod1s = c()
  mod2s = c()
  fparcial1 = c()
  fparcial2 = c()
  fparcial3 = c()
  media = tapply(resp, trat, mean, na.rm = TRUE)
  desvio = tapply(resp, trat, sd, na.rm = TRUE)
  erro = tapply(resp, trat, sd, na.rm = TRUE)/sqrt(table(trat))
  dose = tapply(trat, trat, mean, na.rm = TRUE)
  moda = lm(resp ~ trat)
  mod1a = lm(resp ~ trat + I(trat^2))
  mod2a = lm(resp ~ trat + I(trat^2) + I(trat^3))
  mods = summary(moda)$coefficients
  mod1s = summary(mod1a)$coefficients
  mod2s = summary(mod2a)$coefficients
  modm = lm(media ~ dose)
  mod1m = lm(media ~ dose + I(dose^2))
  mod2m = lm(media ~ dose + I(dose^2) + I(dose^3))
  modf1 = lm(resp ~ trat)
  modf2 = lm(resp ~ trat + I(trat^2))
  modf3 = lm(resp ~ trat + I(trat^2) + I(trat^3))
  modf1ql = anova(modf1)
  modf2ql = anova(modf2)
  modf3ql = anova(modf3)
  modf1q = aov(resp ~ as.factor(trat))
  res = anova(modf1q)
  fadj1 = anova(modf1, modf1q)[2, c(3, 4, 5, 6)]
  fadj2 = anova(modf2, modf1q)[2, c(3, 4, 5, 6)]
  fadj3 = anova(modf3, modf1q)[2, c(3, 4, 5, 6)]
  if (is.na(DFres) == TRUE) {
    DFres = res[2, 1]
  }
  if (is.na(SSq) == TRUE) {
    SSq = res[2, 2]
  }
  df1 = c(modf3ql[1, 1], fadj1[1, 1], DFres)
  df2 = c(modf3ql[1:2, 1], fadj2[1, 1], DFres)
  df3 = c(modf3ql[1:3, 1], fadj3[1, 1], DFres)
  sq1 = c(modf3ql[1, 2], fadj1[1, 2], SSq)
  sq2 = c(modf3ql[1:2, 2], fadj2[1, 2], SSq)
  sq3 = c(modf3ql[1:3, 2], fadj3[1, 2], SSq)
  qm1 = sq1/df1
  qm2 = sq2/df2
  qm3 = sq3/df3
  if (grau == "1") {
    fa1 = data.frame(cbind(df1, sq1, qm1))
    fa1$f1 = c(fa1$qm1[1:2]/fa1$qm1[3], NA)
    fa1$p = c(pf(fa1$f1[1:2], fa1$df1[1:2], fa1$df1[3], lower.tail = F),
              NA)
    colnames(fa1) = c("DF", "SSQ", "MSQ", "F", "p-value")
    rownames(fa1) = c("Linear", "Deviations", "Residual")
  }
  if (grau == "2") {
    fa2 = data.frame(cbind(df2, sq2, qm2))
    fa2$f2 = c(fa2$qm2[1:3]/fa2$qm2[4], NA)
    fa2$p = c(pf(fa2$f2[1:3], fa2$df2[1:3], fa2$df2[4], lower.tail = F),
              NA)
    colnames(fa2) = c("DF", "SSQ", "MSQ", "F", "p-value")
    rownames(fa2) = c("Linear", "Quadratic", "Deviation",
                      "Residual")
  }
  if (grau == "3") {
    fa3 = data.frame(cbind(df3, sq3, qm3))
    fa3$f3 = c(fa3$qm3[1:4]/fa3$qm3[5], NA)
    fa3$p = c(pf(fa3$f3[1:4], fa3$df3[1:4], fa3$df3[5], lower.tail = F),
              NA)
    colnames(fa3) = c("DF", "SSQ", "MSQ", "F", "p-value")
    rownames(fa3) = c("Linear", "Quadratic", "Cubic", "Deviation",
                      "Residual")
  }
  if (grau == "1") {
    r2 = round(summary(modm)$r.squared, 6)
  }
  if (grau == "2") {
    r2 = round(summary(mod1m)$r.squared, 6)
  }
  if (grau == "3") {
    r2 = round(summary(mod2m)$r.squared, 6)
  }
  if (grau == "1") {
    if (is.na(n) == FALSE) {
      coef1 = round(coef(moda)[1], n)
    }
    else {
      coef1 = coef(moda)[1]
    }
    if (is.na(n) == FALSE) {
      coef2 = round(coef(moda)[2], n)
    }
    else {
      coef2 = coef(moda)[2]
    }
    s1 = s <- sprintf("%s == %e %s %e*%s ~~~~~ italic(R^2) == %0.6f",
                      yname.poly, coef1, ifelse(coef2 >= 0, "+", "-"),
                      abs(coef2), xname.poly, r2)
  }
  if (grau == "2") {
    if (is.na(n) == FALSE) {
      coef1 = round(coef(mod1a)[1], n)
    }
    else {
      coef1 = coef(mod1a)[1]
    }
    if (is.na(n) == FALSE) {
      coef2 = round(coef(mod1a)[2], n)
    }
    else {
      coef2 = coef(mod1a)[2]
    }
    if (is.na(n) == FALSE) {
      coef3 = round(coef(mod1a)[3], n)
    }
    else {
      coef3 = coef(mod1a)[3]
    }
    s2 = s <- sprintf("%s == %e %s %e * %s %s %e * %s^2 ~~~~~ italic(R^2) ==  %0.6f",
                      yname.poly, coef1, ifelse(coef2 >= 0, "+", "-"),
                      abs(coef2), xname.poly, ifelse(coef3 >= 0, "+", "-"),
                      abs(coef3), xname.poly, r2)
  }
  if (grau == "3") {
    if (is.na(n) == FALSE) {
      coef1 = round(coef(mod2a)[1], n)
    }
    else {
      coef1 = coef(mod2a)[1]
    }
    if (is.na(n) == FALSE) {
      coef2 = round(coef(mod2a)[2], n)
    }
    else {
      coef2 = coef(mod2a)[2]
    }
    if (is.na(n) == FALSE) {
      coef3 = round(coef(mod2a)[3], n)
    }
    else {
      coef3 = coef(mod2a)[3]
    }
    if (is.na(n) == FALSE) {
      coef4 = round(coef(mod2a)[4], n)
    }
    else {
      coef4 = coef(mod2a)[4]
    }
    s3 = s <- sprintf("%s == %e %s %e * %s %s %e * %s^2 %s %0.e * %s^3 ~~~~~ italic(R^2) == %0.6f",
                      yname.poly, coef1, ifelse(coef2 >= 0, "+", "-"),
                      abs(coef2), xname.poly, ifelse(coef3 >= 0, "+", "-"),
                      abs(coef3), xname.poly, ifelse(coef4 >= 0, "+", "-"),
                      abs(coef4), xname.poly, r2)
  }

  ################################################### Graphic part ################################################################
  residuoJKF<-function(grau){
    if (grau == "1") {
      yp <- predict.lm(moda)      #residuo Jaknife
      Jack<-rstudent(moda)
    }
    if (grau == "2") {
      yp <- predict.lm(mod1a)      #residuo Jaknife
      Jack<-rstudent(mod1a)
    }
    if (grau == "3") {
      yp <- predict.lm(mod2a)      #residuo Jaknife
      Jack<-rstudent(mod2a)
    }
    data.frame(yp,Jack)
  }
  dadoJKF=residuoJKF(grau)
  data1=data.frame(trat,resp)
  data1=data.frame(trat=dose,#as.numeric(as.character(names(media))),
                   resp=media,
                   desvio, erro)
  data2=data.frame(yp=dadoJKF$yp,yjack=dadoJKF$Jack)

  grafico=ggplot(data1,aes(x=trat,y=resp))
  if(point=="all"){grafico=grafico+
    geom_point(data=dados,
               aes(y=resp,x=trat),shape=21,
               fill=color,color="black")}
  if(point=="mean_sd"){grafico=grafico+
    geom_errorbar(aes(ymin=resp-desvio,ymax=resp+desvio),width=width.bar,size=linesize)}
  if(point=="mean_se"){grafico=grafico+
    geom_errorbar(aes(ymin=resp-erro,ymax=resp+erro),width=width.bar,size=linesize)}
  if(point=="mean"){grafico=grafico}
  grafico=grafico+geom_point(aes(fill=as.factor(rep(1,length(resp)))),na.rm=TRUE,
                             size=pointsize,shape=21,
                             color="black")+
    theme+ylab(ylab)+xlab(xlab)
  if(is.na(ylim[1])==TRUE){grafico=grafico}else{grafico=grafico+ylim(ylim)}

  if(grau=="0"){grafico=grafico+geom_line(y=mean(resp),size=linesize,lty=2)}
  if(grau=="1"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x,size=linesize,color="black")}
  if(grau=="2"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2),size=linesize,color="black")}
  if(grau=="3"){grafico=grafico+geom_smooth(method = "lm",se=se, na.rm=TRUE, formula = y~x+I(x^2)+I(x^3),size=linesize,color="black")}
  if(grau=="0"){grafico=grafico+
    scale_fill_manual(values=color,label=paste("y =",round(mean(resp),3)),name="")}
  if(grau=="1"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s1)),name="")}
  if(grau=="2"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s2)),name="")}
  if(grau=="3"){grafico=grafico+
    scale_fill_manual(values=color,label=c(parse(text=s3)),name="")}

  if(color=="gray"){if(grau=="1"){grafico=grafico+
    scale_fill_manual(values="black",label=c(parse(text=s1)),name="")}
    if(grau=="2"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s2)),name="")}
    if(grau=="3"){grafico=grafico+
      scale_fill_manual(values="black",label=c(parse(text=s3)),name="")}
  }

  grafico=grafico+
    theme(text = element_text(size=textsize,color="black",family=family),
          axis.text = element_text(size=textsize,color="black",family=family),
          axis.title = element_text(size=textsize,color="black",family=family),
          legend.position = posi,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)

  ############
  ############
  nt<-length(resp)
  grafico2=ggplot(data2,aes(x=yp,y=yjack))+
    geom_point()+#data2,aes(x=yp,y=yjack),shape=21,color="black")+
    geom_hline(yintercept =qt(0.975,nt-2), linetype = "dashed",
               color="red",size=0.8)+
    geom_hline(yintercept =0, linetype = "dashed",
               color="gray80",size=0.8)+
    geom_hline(yintercept =qt(0.025,nt-2), linetype = "dashed",
               color="red",size=0.8)+
    labs(x='Residue',y='Jacknife Residue')+
    theme_light()

  graficos<-ggarrange(grafico,grafico2,ncol = 2, nrow = 1)

  print(graficos)


  ######################################################### degree of regression #########################################################


  regressao<-function(grau){
    if (grau == 1) {
      cat("\n----------------------------------------------------\n")
      cat("Regression model")
      cat("\n----------------------------------------------------\n")
      print(mods)

      cat("\n----------------------------------------------------\n")
      cat("Regression biases")
      cat("\n----------------------------------------------------\n")
      print(as.matrix(fa1), na.print = " ")

      if (fa1[1, "p-value"] < 0.05) {
        cat("\n----------Considerations----------")
        cat("\nAs p-value (Linear) < alpha, the linear model describes the data, with 95% confidence\n")
    }
        if (fa1[2, "p-value"] > 0.05) {
        cat("\nAs p-value (Deviations) > alpha, there are no deviations from linearity, with 95% confidence\n")
        }
  }
    if (grau == 2) {
      cat("\n----------------------------------------------------\n")
      cat("Regression Models")
      cat("\n----------------------------------------------------\n")
      print(mod1s)
      cat("\n----------------------------------------------------\n")
      cat("Deviations from regression")
      cat("\n----------------------------------------------------\n")
      print(as.matrix(fa2), na.print = " ")
    }
    if (grau == 3) {
      cat("\n----------------------------------------------------\n")
      cat("Regression Models")
      cat("\n----------------------------------------------------\n")
      print(mod2s)
      cat("\n----------------------------------------------------\n")
      cat("Deviations from regression")
      cat("\n----------------------------------------------------\n")
      print(as.matrix(fa3), na.print = " ")
    }
  }
  modelo<-regressao(grau)

  ########################################################### Normality test of model residuals ##############################################################
  #grau 1
  if (grau == "1") {
    pvalor.shapiro <- shapiro.test(modm$residuals)$p.value
    normal_value <- shapiro.test(modm$residuals)$statistic
    cat("\n------------------------------------------------------------------------\nShapiro-Wilk test for normality of residuals \n")
    cat("p-value: ", pvalor.shapiro, "\n")
    cat("W: ", normal_value, "\n")
    if (pvalor.shapiro <= 0.05) {
      cat("At 95% confidence, the residues cannot be considered normal.\n------------------------------------------------------------------------\n")
    }
    if (pvalor.shapiro > 0.05) {
      cat("According to the Shapiro-Wilk test at 95% confidence, the residuals can be considered normal.\n------------------------------------------------------------------------\n")
    }}
  #grau 2
  if (grau == "2") {
    pvalor.shapiro <- shapiro.test(mod1m$residuals)$p.value
    cat("\n------------------------------------------------------------------------\nTeste de normalidade dos residuos (Shapiro-Wilk)\n")
    cat("p-value: ", pvalor.shapiro, "\n")
    if (pvalor.shapiro <= 0.05) {
      cat("At 5% significance, the residues cannot be considered normal.\n------------------------------------------------------------------------\n")
    }
    if (pvalor.shapiro > 0.05) {
      cat("According to the Shapiro-Wilk test at 95% confidence, the residuals can be considered normal.\n------------------------------------------------------------------------\n\n")
    }}
  #grau 3
  if (grau == "3") {
    pvalor.shapiro <- shapiro.test(mod2m$residuals)$p.value
    cat("\n------------------------------------------------------------------------\nTeste de normalidade dos residuos (Shapiro-Wilk)\n")
    cat("p-value: ", pvalor.shapiro, "\n")
    if (pvalor.shapiro <= 0.05) {
      cat("According to the Shapiro-Wilk test at 95% confidence, the residuals can be considered normal.\n------------------------------------------------------------------------\n")
    }
    if (pvalor.shapiro > 0.05) {
      cat("According to the Shapiro-Wilk test at 95% confidence, the residuals can be considered normal.\n------------------------------------------------------------------------\n\n")
    }}
  ########################################################################### Homogeneity of variances  #############################################################
  pvalor.levene <- leveneTest(modf1q)$`Pr(>F)`[1]
  valor_F <- leveneTest(modf1q)$`F value`[1]
  cat("\n------------------------------------------------------------------------\nBrown-Forsythe Test for Homogeneity of Variances\n")
  cat("p-value: ", pvalor.levene, "\n")
  cat("F-value: ", valor_F, "\n")
  if (pvalor.levene <= 0.05) {
    cat("At 95% confidence, the variances cannot be considered homogeneous.\n------------------------------------------------------------------------\n")
    }
  if (pvalor.levene > 0.05) {
    cat("According to the Brown-Forsythe test, at 95% confidence, the variances can be considered homogeneous.\n------------------------------------------------------------------------\n")
  }

  ########################################################### Test of independence of model residuals ##############################################################
  #grau 1
  if (grau == "1") {
    pvalor.DW <-durbinWatsonTest(modm)$p
    valor_DW <- durbinWatsonTest(modm)$dw
    cat("\n------------------------------------------------------------------------\nDurbin-Watson Test for Independence of Residuals\n")
    cat("p-value: ", pvalor.DW, "\n")
    cat("dw-value: ", valor_DW, "\n")
    if (pvalor.DW <= 0.05) {
      cat("At 5% significance, the residues are dependent.\n------------------------------------------------------------------------\n")
    }
    if (pvalor.levene > 0.05) {
      cat("According to the Durbin-Watson test at 95% confidence, the residuals are independent.\n------------------------------------------------------------------------\n")
    }
  }
  #grau 2
  if (grau == "2") {
    pvalor.DW <-durbinWatsonTest(mod1m)$p
    cat("\n------------------------------------------------------------------------\nTeste de Durbin-Watson para Independencia dos residuos\n")
    cat("p-value: ", pvalor.DW, "\n")
    if (pvalor.DW <= 0.05) {
      cat("At 5% significance, the residues are dependent\n------------------------------------------------------------------------\n")
    }
    if (pvalor.levene > 0.05) {
      cat("According to the Durbin-Watson test at 95% confidence, the residuals are independent.\n------------------------------------------------------------------------\n\n")
    }  }
  #grau 3
  if (grau == "3") {
    pvalor.DW <-durbinWatsonTest(mod2m)$p
    cat("\n------------------------------------------------------------------------\nTeste de Durbin-Watson para Independencia dos residuos\n")
    cat("p-value: ", pvalor.DW, "\n")
    if (pvalor.DW <= 0.05) {
      cat("At 5% significance, the residues are dependent.\n------------------------------------------------------------------------\n")
    }
    if (pvalor.levene > 0.05) {
      cat("According to the Durbin-Watson test at 95% confidence, the residuals are independent.\n------------------------------------------------------------------------\n\n")
    }  }


  ################################################### Residuo ######################################################################
  pontosJKF<-function(grau){
    if (grau == "1") {
      yp <- predict.lm(moda)      #residuo Jaknife
      Jack<-rstudent(moda)
      n<-length(resp)
      QQ<-data.frame(resp,Jack,yp,ni=1:n)
      dd<-QQ[order(abs(QQ[,2]), decreasing=TRUE),]
      pts<-ifelse(abs(dd[,2])>qt(0.975,n-2),1,0)*dd$ni
    }
    if (grau == "2") {
      yp <- predict.lm(mod1a)      #residuo Jaknife
      Jack<-rstudent(mod1a)
      n<-length(resp)
      QQ<-data.frame(resp,Jack,yp,ni=1:n)
      dd<-QQ[order(abs(QQ[,2]), decreasing=TRUE),]
      pts<-ifelse(abs(dd[,2])>qt(0.975,n-2),1,0)*dd$ni
    }
    if (grau == "3") {
      yp <- predict.lm(mod2a)      #residuo Jaknife
      Jack<-rstudent(mod2a)
      n<-length(resp)
      QQ<-data.frame(resp,Jack,yp,ni=1:n)
      dd<-QQ[order(abs(QQ[,2]), decreasing=TRUE),]
      pts<-ifelse(abs(dd[,2])>qt(0.975,n-2),1,0)*dd$ni
    }
    # retirando os bugs de nao haver pontos discrepantes no codigo
    soma<-sum(pts)
    cat("Jacknife Test \n")
    #cat('\n')
    if (soma == 0) {
      cat("There are no outliers according to the Jaknife residue \n------------------------------------------------------------------------\n")
    }
    if (soma > 0) {
      pontos<-pts[min(which(pts!=0)):max(which(pts!=0))]
      cat("According to the Jacknife residue, the outliers are:", pontos,"\n------------------------------------------------------------------------\n")
    }

  }
  Jackres<-pontosJKF(grau)

  #Saida
  out<-list()
  out$modelo
  out$pvalor.shapiro
  out$pvalor.levene
  out$Jackres
  out$graficos
  invisible(out)
}
