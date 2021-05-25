## ---------------------------
##
## Script name: He-Hegarty2020DataAnalysis
##
## Purpose of script: Linear analysis and path model development
##
## Author: Chuanxiuyue (Carol) He
##
## Date Created: 2020-02-18
##
## Email: carol.hcxy@gmail.com
## 
## Content: 
##      1. Load Study 1 Data
##      2. Descriptive Stats (Table1)
##      3. Correlation Tables (Table2)
##      4. GGM - GMN paired t-test
##      5. Partial Correlation
##      6. Load Study 2 Data
##      7. Descriptive Stats (Table3)
##      8. GGM - GMN paired t-test
##      9. Compare Study 1 & Study 2
##      10. Correlation Tables (Table4)
##      11. Model 1: regression
##      12. Model 2.1
##      13. Model 2.2
## ---------------------------

## ---------------------------

## set working directory for Mac and PC

setwd("~/CarolHe/")      # Carol's working directory (mac)
setwd("C:/Users/carolhe/")    # if using PC

## ---------------------------

## load up the packages
## use install.packages() to install new packages
library(psych) ##mediation model & descriptive stats
library(tidyverse)
library(rstatix)## cohen's d

library(lm.beta)## standarized coefficients
library(ppcor)## semi-partial correlations
library(olsrr)## VIF

library(moments)
library(xlsx)
library(xtable)

## ---------------------------

#####----------Load Study 1 Data----------#####
df1 <- read.csv("HE&Hegarty2020-1.csv")
df <- df1
df[,c(1,3)]<- lapply(df[,c(1,3)], factor)
## view variable names
names(df)

#####----------Descriptive Stats (Table1)----------#####

## sex distribution
table(df$Sex)

## raw data distribution
psych::describe(df[,c(2,6,8:14)])

## OPT & MPT are skewed. Remedy it by using log transformation
df$LOPT <- log(df$OPT)
df$LMPT <- log(df$MPT)

## descriptive stats table
psych::describe(df[,c(6,8:16)])

#####----------Correlation Table (Table2)----------#####

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], use="complete.obs")
      p.mat[i, j] <- p.mat[j, i] <- round(tmp$p.value,4)
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

corr_S1<-round(cor(df[,c(9:16)]),2)
corr_p_S1 <- cor.mtest(df[,c(9:16)])
colnames(p.mat_35)=c('GGM','GMN','SBSOD','SA','ET','GPS','LPOT','LMPT')
rownames(p.mat_35)=c('GGM','GMN','SBSOD','SA','ET','GPS','LPOT','LMPT')
corr_S1
corr_p_S1

#####----------GGM - GMN paired t-test----------#####
t.test(df$GM,df$GNT,paired=T)

## Cohen's d
gather(df[,c('GM','GNT')],
       test,
       score,
       GM:GNT,
       factor_key = T)%>%
  rstatix::cohens_d(score~test,paired=T)

#####----------Partial Correlation----------#####
pcor.test(df$GM,df$SBSOD,df$GNT,method="pearson")
pcor.test(df$GNT,df$SBSOD,df$GM,method="pearson")
pcor.test(df$GM,df$SA,df$GNT,method="pearson")
pcor.test(df$GNT,df$SA,df$GM,method="pearson")


#####----------Load Study 2 Data----------#####
df <- read.csv("HE&Hegarty2020-2.csv")
df[,c(2,3)]<- lapply(df[,c(2,3)], factor)
## view variable names
names(df)

#####----------Descriptive Stats (Table3)----------#####

## sex distribution
table(df$Gender)

## OPT is skewed. Remedy it by using log transformation
df$LOPT <- log(df$OPT)

## descriptive stats table
psych::describe(df[,c(4:12)])


#####----------GGM - GMN paired t-test----------#####
t.test(df$GGM,df$GMN,paired=T)

## Cohen's d
gather(df[,c('GGM','GMN')],
       test,
       score,
       GGM:GMN,
       factor_key = T)%>%
  rstatix::cohens_d(score~test,paired=T)


#####----------Compare Study1 & Study2----------#####
t.test(df1$GNT,df$GMN,var.equal = T)

## Cohens'D
df1$Study <- 'Study1'
df$Study <- 'Study2'
test1 <- df1[,c('GNT','Study')]
colnames(test1) <- c('GMN','Study')
test <- rbind(test1,df[,c('GMN','Study')])
test%>%
  rstatix::cohens_d(GMN~Study,paired=F)

#####----------Correlation Table (Table4)----------#####

corr_S2<-round(cor(df[,c(5:10,12)],use="complete.obs"),2)
corr_p_S2 <- cor.mtest(df[,c(5:10,12)])
corr_S2
corr_p_S2

diff.corr <- function( r1, n1, r2, n2 ){
  
  Z1 <- 0.5 * log( (1+r1)/(1-r1) )
  Z2 <- 0.5 * log( (1+r2)/(1-r2) )
  
  diff   <- Z1 - Z2
  SEdiff <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
  diff.Z  <- diff/SEdiff
  
  p <- 2*pnorm( abs(diff.Z), lower=F)
  cat( "Two-tailed p-value", p , "\n" )
}

diff.corr( r1=0.6, n1=123, r2=0.42, n2=149)

#####----------Model 1----------#####
dfnm <- na.omit(df)
m1 <- lm(SBSOD~GMN+SA+ET+GPS,data=dfnm)
summary(m1)
confint(m1)
lm.beta(m1)
(spcor(cbind(dfnm$SBSOD,dfnm$GMN,dfnm$SA,dfnm$ET,dfnm$GPS))$estimate)^2
ols_vif_tol(m1)

#####----------Model 2.1----------#####
set.seed(12345)
m2.1.1 <- mediate(SBSOD~SA+(ET)+(GPS), data = df, std = T, n.iter = 5000)
summary(m2.1.1,digits = 4)

m2.1.2 <- mediate(SBSOD~GMN+(ET)+(GPS), data = df, std = T, n.iter = 5000)
summary(m2.1.2,digits = 4)

#####----------Model 2.2----------#####
m2.2 <- mediate(GPS+ET~SBSOD+(SA), data = df, std = T, n.iter = 5000)
summary(m2.2,digits = 4)

