

##===============================  Respiratory Gases for Tonle Sap Lake  ===============================##


##  Author(s):  B. Miller

##  Required dataframes:  Master_TSL-Corrected.csv, metab.csv

#rm(list=ls())

#install.packages("readr")
library(readr) 
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("mgcv")
library(mgcv)
#install.packages("AICcmodavg")
library(AICcmodavg)
#install.packages("rLakeAnalyzer")
library(rLakeAnalyzer) # Produced by GLEON
#install.packages("lubridate")
library(lubridate) # Includes functions for time
#install.packages("LakeMetabolizer")
library(LakeMetabolizer) # Produced by GLEON
#install.packages("plotrix")
library(plotrix)
install.packages("rcompanion")
library(rcompanion)

setwd("~/Desktop/Working")


###############################################
##  Table for Wilcoxon Tests & Effect Sizes  ##
###############################################

Table <- function(x) { # Calculate means and standard errors
  n <- length(na.omit(x))
  Mean <- mean(na.omit(x))
  SE <- sd(na.omit(x))/sqrt(length(na.omit(x))) 
  Results <- list(n, Mean, SE)
  return(Results)
}

d <- function(x,y) { # Calculate Cohen's d for effect size
  xM <- mean(na.omit(x))
  yM <- mean(na.omit(y))
  sdx <- sd(na.omit(x))
  sdy <- sd(na.omit(y))
  SD <- sqrt(((sdx^2)+(sdy^2))/2)
  d <- (xM-yM)/SD
  return(d)
}

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, D13C_CO2, D13C_CH4, ENVIRON, STAGE) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCH4 = mean(PCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CH4 = mean(D13C_CH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(eC = avg.D13C_CO2 - avg.D13C_CH4) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  ungroup() 
cGas <- cGas[!cGas$SITE == "KgPr", ]
cGas <- cGas[!cGas$SITE == "StSe", ]
cGas <- cGas[!cGas$SITE == "PtSd.O", ] 
cGas <- cGas[!cGas$SITE == "PtSd.E", ] 
cGas <- cGas[!cGas$SITE == "PtSd.F", ] 
cGas <- cGas[!cGas$SITE == "KgKl.O", ] 
cGas <- cGas[!cGas$SITE == "KgKl.E", ] 
cGas <- cGas[!cGas$STAGE == "Rising", ] 
cGas <- cGas[complete.cases(cGas), ]
#View(cGas)

# Determine whether data follow normal or non-normal distributions in order to choose appropriate statistical test

  # PCO2 by flood stage

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(High$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(High$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(High$avg.PCO2) # p-value = 0.05783 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Falling$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Falling$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Transform
shapiro.test(Falling$avg.PCO2) # p-value = 0.6506 # Distribution of data is not significantly different from normal distribution 

df <- group_by(cGas, STAGE) %>%
  dplyr::select(avg.PCO2)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ STAGE, data = df)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov) # p-value = 5e-07 # High avg.PCO2 is significantly different from Falling avg.PCO2

Table(High$avg.PCO2)
Table(Falling$avg.PCO2)
d(High$avg.PCO2, Falling$avg.PCO2)

  # PCO2 by depth

Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Surface$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Surface$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Surface$avg.PCO2) # p-value = 0.3565 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Bottom$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Bottom$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Bottom$avg.PCO2) # p-value = 0.02027 # Distribution of data is significantly different from normal distribution 

df <- group_by(cGas, CLASS_Z) %>%
  dplyr::select(avg.PCO2) %>%
  mutate(avg.PCO2_log <- log(abs(avg.PCO2)))
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov) # p-value = 0.2336275 # Surface avg.PCO2 is not significantly different from Bottom avg.PCO2

df <- group_by(Falling, CLASS_Z) %>%
  dplyr::select(avg.PCO2)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # p-value = 0.4500919

  # PCO2 by lake environment

High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")

Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.PCO2) # p-value = 0.156 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.PCO2) # p-value = 0.4802 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Floodplain$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Floodplain$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Floodplain$avg.PCO2) # p-value = 0.07363 # Distribution of data is not significantly different from normal distribution 

df <- group_by(High, ENVIRON) %>%
  dplyr::select(avg.PCO2)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) # p-value = 0.904
TukeyHSD(res.aov) 

Table(Open$avg.PCO2)

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

df <- group_by(Falling, ENVIRON) %>%
  dplyr::select(avg.PCO2)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov) 



summary(aov(High$avg.PCO2, Falling$avg.PCO2))

  # PCH4 data follow a non-normal distribution, even after transformation and accounting for a Bonferroni-corrected critical alpha-value of 0.025

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(cGas$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(cGas$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(cGas$avg.PCH4) # p-value = 0.1148 # Distribution of data are significantly different from normal distribution 

PCH4_sqrt <- sqrt(cGas$avg.PCH4)
qqnorm(PCH4_sqrt, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(PCH4_sqrt, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(PCH4_sqrt) # p-value = 5.502e-08 # Distribution of data are significantly different from normal distribution 

PCH4_cbrt <- sign(cGas$avg.PCH4) * abs(cGas$avg.PCH4)^(1/3)
qqnorm(PCH4_cbrt, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(PCH4_cbrt, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(PCH4_cbrt) # p-value = 0.002504 # Distribution of data are significantly different from normal distribution 

PCH4_log <- log(cGas$avg.PCH4) 
qqnorm(PCH4_log, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(PCH4_log, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(PCH4_log) # p-value = 0.002719 # Distribution of data are significantly different from normal distribution 

PCH4_tuk<- transformTukey(cGas$avg.PCH4) 
qqnorm(PCH4_tuk, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(PCH4_tuk, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(PCH4_tuk) # p-value = p-value = 0.01891 # Distribution of data are significantly different from normal distribution 

  # D13C_CO2 data follow a normal distribution after log-transformation

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(cGas$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(cGas$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(cGas$avg.D13C_CO2) # p-value = 0.001736 # Distribution of data are significantly different from normal distribution 

D13C_CO2_sqrt <- sqrt(abs(cGas$avg.D13C_CO2))
qqnorm(D13C_CO2_sqrt, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(D13C_CO2_sqrt, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(D13C_CO2_sqrt) # p-value = 0.01971 # Distribution of data are significantly different from normal distribution 

D13C_CO2_cbrt <- sign(cGas$avg.D13C_CO2) * abs(cGas$avg.D13C_CO2)^(1/3)
qqnorm(D13C_CO2_cbrt, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(D13C_CO2_cbrt, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(D13C_CO2_cbrt) # p-value = 0.03931 # Distribution of data are significantly different from normal distribution 

D13C_CO2_log <- log(abs(cGas$avg.D13C_CO2)) 
qqnorm(D13C_CO2_log, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(D13C_CO2_log, col = "gray40", lwd=3, lty=1) # Data are not skewed # Use ANOVA to compare groups
shapiro.test(D13C_CO2_log) # p-value = 0.1209 # Distribution of data are not significantly different from normal distribution 

  # D13C_CH4 data follow a normal distribution after sqrt-transformation

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(cGas$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(cGas$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed # Use ANOVA to compare groups
shapiro.test(cGas$avg.D13C_CH4) # p-value = 0.001736 # Distribution of data are not significantly different from normal distribution 

D13C_CH4_sqrt <- sqrt(abs(cGas$avg.D13C_CH4))
qqnorm(D13C_CH4_sqrt, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(D13C_CH4_sqrt, col = "gray40", lwd=3, lty=1) # Data are skewed # Use Kruskal-Wallis to compare groups
shapiro.test(D13C_CH4_sqrt) # p-value = 0.181 # Distribution of data are significantly different from normal distribution 



Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")
wilcox.test(Surface$avg.D13C_CH4, Bottom$avg.D13C_CH4)
d(Surface$avg.D13C_CH4, Bottom$avg.D13C_CH4)


High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")
query <- subset(cGas, STAGE=="High")
Open <- subset(query, ENVIRON=="Pelagic")
Edge <- subset(query, ENVIRON=="Edge")
Floodplain <- subset(query, ENVIRON=="Floodplain")

Table(query$avg.PCO2)
Table(query$avg.PCH4)
Table(query$avg.D13C_CO2)
Table(query$avg.D13C_CH4)

Table(Open$avg.PCO2)
Table(Open$avg.PCH4)
Table(Open$avg.D13C_CO2)
Table(Open$avg.D13C_CH4)

Table(Edge$avg.PCO2)
Table(Edge$avg.PCH4)
Table(Edge$avg.D13C_CO2)
Table(Edge$avg.D13C_CH4)

Table(Floodplain$avg.PCO2)
Table(Floodplain$avg.PCH4)
Table(Floodplain$avg.D13C_CO2)
Table(Floodplain$avg.D13C_CH4)

wilcox.test(Open$avg.PCH4, Edge$avg.PCH4)
d(Open$avg.PCH4, Edge$avg.PCH4)

setwd("~/Desktop/Diel TSL O2")

# Load the data

metab <- read.csv("metab.csv", header=T)

query <- subset(metab, STAGE=="Falling")
Open <- subset(query, DESCRIP=="Open")
Edge <- subset(query, DESCRIP=="Edge")
Floodplain <- subset(query, DESCRIP=="Forest")

High <- subset(metab, STAGE=="High")
Falling <- subset(metab, STAGE=="Falling")

Table(Floodplain$NEP)
wilcox.test(High$NEP, Falling$NEP)
d(Floodplain$PCH4, Open$PCH4)


############################################
##   Scatterplots for O2_mmol & CO2_mmol  ##
############################################

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, O2, KH_O2, gradConcCO2, ACTUAL_Z) %>%
  mutate(iCO2_mmolL = gradConcCO2*1000) %>% # Excess CO2 in mmol^L-1
  mutate(iO2_mmolL = O2 - ((KH_O2*(209500/1000000))*32)) # O2 deficit in mmol^L-1
#cGas <- cGas[cGas$CLASS_Z == "Surface", ]
cGas <- cGas[!cGas$SITE == "StSe", ] 
cGas <- cGas[!cGas$SITE == "KgPr", ] 
cGas <- cGas[!cGas$SITE == "PtSd.O", ]
cGas <- cGas[!cGas$SITE == "PtSd.E", ] 
cGas <- cGas[!cGas$SITE == "PtSd.F", ] 
cGas <- cGas[!cGas$SITE == "KgKl.O", ] 
cGas <- cGas[!cGas$SITE == "KgKl.E", ] 
cGas <- cGas[!cGas$STAGE == "Rising", ] 
cGas <- cGas[complete.cases(cGas), ]
cGas <- cGas %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, iO2_mmolL, iCO2_mmolL, ACTUAL_Z) %>%
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(CO2_mmolL = mean(iCO2_mmolL)) %>%
  mutate(O2_mmolL = mean(iO2_mmolL)) %>% 
  distinct(CO2_mmolL, .keep_all=TRUE) %>%
  ungroup() 
#View(cGas)

High <- cGas[cGas$STAGE == "High", ]
Falling <- cGas[cGas$STAGE == "Falling", ]

new_High <- data.frame(CO2_mmolL = seq(min(High$CO2_mmolL), max(High$CO2_mmolL), length.out = 100)) 
lm_y <- lm(O2_mmolL ~ CO2_mmolL, data=High)
summary(lm_y) #  -0.08998
pred_High <- predict(lm_y, newdata=new_High, se.fit=TRUE, type="response")
new_High$pred <- predict(lm_y, newdata=new_High)
#plot(High$CO2_mmolL, resid(lm_y))

new_Falling <- data.frame(CO2_mmolL = seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100)) 
lm_y <- lm(O2_mmolL ~ CO2_mmolL, data=Falling)
summary(lm_y) #  -0.07092
pred_Falling <- predict(lm_y, newdata=new_Falling, se.fit=TRUE, type="response")
new_Falling$pred <- predict(lm_y, newdata=new_Falling)
#plot(Falling$CO2_mmolL, resid(lm_y))

par(mfrow=c(3, 2))
par(oma=c(4, 4, 4, 4))
par(mar=c(4, 8, 4, 2)) # Bottom, left, top, right
#par(oma=c(2, 2, 2, 2))
#par(mar=c(9, 9, 9, 9)) # Bottom, left, top, right
range(cGas$CO2_mmolL) #  -0.5232426 40.3720710
range(cGas$O2_mmolL) # -8.5488738  0.7144784
par(xpd=FALSE)
plot(High$CO2_mmolL, High$O2_mmolL, xlim=c(-5, 41), ylim=c(-15, 5), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

abline(0, -1, col="black", lty=2, lwd=2)
abline(h=0, col="darkorange", lwd=2)
abline(v=0, col="darkorange", lwd=2)

#polygon(c(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100), rev(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100))), 
        #c(pred_Falling$fit - (1.96 * pred_Falling$se.fit), rev(pred_Falling$fit + (1.96 * pred_Falling$se.fit))), 
        #col=adjustcolor("blue3", alpha.f=0.4), border = NA) # Falling
#lines(new_Falling$CO2_mmolL, new_Falling$pred, lwd=2, lty=2, col="blue3") 

polygon(c(seq(min(High$CO2_mmolL), max(High$CO2_mmolL), length.out = 100), rev(seq(min(High$CO2_mmolL), max(High$CO2_mmolL), length.out = 100))), 
        c(pred_High$fit - (1.96 * pred_High$se.fit), rev(pred_High$fit + (1.96 * pred_High$se.fit))), 
        col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) # High
lines(new_High$CO2_mmolL, new_High$pred, lwd=2, lty=2, col="cornflowerblue") 

mtext(expression(O[2] ~ Deficit ~ (mmol ~ L^{-1})), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(Excess ~ CO[2] ~ (mmol ~ L^{-1})), side=1, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-350, 50), col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-5, 45), las=1, cex=2, cex.axis=2, cex.lab=2, font=1)  #las=1 makes horizontal labels

#points(jitter(Falling$CO2_mmolL, factor=1), Falling$O2_mmolL, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
points(jitter(High$CO2_mmolL, factor=1), High$O2_mmolL, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)

legend("topright", legend=expression("Aerobic" ~ "ER" ~ italic(m) ~ "=" ~ "-1.0"), 
       text.col=c("black"), text.font=3, cex=1.8, bty="n")
legend("bottomright", legend=expression("High-water" ~ italic(m) ~ "=" ~ "-0.1"), 
       text.col="cornflowerblue", text.font=3, cex=1.8, bty="n")

#arrows(-17, -10, -17, -30, lty=1, lwd=3, length=0.1, xpd=TRUE, col="gray40")
#mtext(expression("Undersaturated" ~ O[2]), side=2, at=-20, line=9.5, cex=2.5, font=1, xpd=TRUE, col="gray40")

#arrows(10, -53, 30, -53, lty=1, lwd=3, length=0.1, xpd=TRUE, col="gray40")
#mtext(expression("Oversaturated" ~ CO[2]), side=1, at=20, line=9.5, cex=2.5, font=1, xpd=TRUE, col="gray40")

par(xpd=TRUE)
#legend("topright", inset=-0.18, legend=c("Open", "Edge", "Floodplain"), 
       #text.col=c("black"), text.font=c(1, 1, 1),  pch=c(21, 22, 25), col="black", bg="white", pt.cex=c(2.5, 2.5, 2.5), cex=1.8, bty="n")
#legend("topleft", inset=-0.18, legend=c("High", "Falling"), 
       #text.col=c("cornflowerblue", "blue3"), text.font=c(1, 1),  col="black", bg="white", pt.cex=c(2.5, 2.5), cex=1.8, bty="n")
legend("topleft", inset=-0.21, legend=expression("a)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(4, 3, 4, 7)) # Bottom, left, top, right
#par(oma=c(2, 2, 2, 2))
#par(mar=c(9, 9, 9, 9)) # Bottom, left, top, right
range(cGas$CO2_mmolL) #  -0.5232426 40.3720710
range(cGas$O2_mmolL) # -8.5488738  0.7144784
par(xpd=FALSE)
plot(cGas$CO2_mmolL, cGas$O2_mmolL, xlim=c(-5, 41), ylim=c(-15, 5), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

abline(0, -1, col="black", lty=2, lwd=2)
abline(h=0, col="darkorange", lwd=2)
abline(v=0, col="darkorange", lwd=2)

polygon(c(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100), rev(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100))), 
        c(pred_Falling$fit - (1.96 * pred_Falling$se.fit), rev(pred_Falling$fit + (1.96 * pred_Falling$se.fit))), 
        col=adjustcolor("blue3", alpha.f=0.4), border = NA) # Falling
lines(new_Falling$CO2_mmolL, new_Falling$pred, lwd=2, lty=2, col="blue3") 

#polygon(c(seq(min(cGas$CO2_mmolL), max(cGas$CO2_mmolL), length.out = 100), rev(seq(min(cGas$CO2_mmolL), max(cGas$CO2_mmolL), length.out = 100))), 
        #c(pred_cGas$fit - (1.96 * pred_cGas$se.fit), rev(pred_cGas$fit + (1.96 * pred_cGas$se.fit))), 
        #col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) # High
#lines(new_cGas$CO2_mmolL, new_cGas$pred, lwd=2, lty=2, col="cornflowerblue") 

#mtext(expression(O[2] ~ Deficit ~ (mmol ~ L^{-1})), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(Excess ~ CO[2] ~ (mmol ~ L^{-1})), side=1, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-350, 50), col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-5, 45), las=1, cex=2, cex.axis=2, cex.lab=2, font=1)  #las=1 makes horizontal labels

points(jitter(Falling$CO2_mmolL, factor=1), Falling$O2_mmolL, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
#points(jitter(High$CO2_mmolL, factor=1), High$O2_mmolL, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)

legend("topright", legend=expression("Atm." ~ CO[2] ~ "&" ~ O[2]), 
       text.col="darkorange", text.font=3, cex=1.8, bty="n")
legend("bottomright", legend=expression("Falling-water" ~ italic(m) ~ "=" ~ "-0.1"), 
       text.col="blue3", text.font=3, cex=1.8, bty="n")

#arrows(-17, -10, -17, -30, lty=1, lwd=3, length=0.1, xpd=TRUE, col="gray40")
#mtext(expression("Undersaturated" ~ O[2]), side=2, at=-20, line=9.5, cex=2.5, font=1, xpd=TRUE, col="gray40")

#arrows(10, -53, 30, -53, lty=1, lwd=3, length=0.1, xpd=TRUE, col="gray40")
#mtext(expression("Oversaturated" ~ CO[2]), side=1, at=20, line=9.5, cex=2.5, font=1, xpd=TRUE, col="gray40")

par(xpd=TRUE)
legend("top", inset=-0.25, legend=c("Open-water", "Edge", "Floodplain"), 
       text.col=c("black"), horiz=TRUE, text.font=c(1, 1, 1),  pch=c(21, 22, 25), col="black", bg="white", pt.cex=c(2.5, 2.5, 2.5), cex=1.8, bty="n")
#legend("topright", inset=-0.26, legend=c("High", "Falling"), 
       #text.col=c("cornflowerblue", "blue3"), text.font=c(1, 1),  col="black", bg="white", pt.cex=c(2.5, 2.5), cex=2.75, bty="n")
legend("topleft", inset=-0.21, legend=expression("b)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


#################################################
##   Scatterplots for avg.PCH4 & avg.D13C-CO2  ##
#################################################

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>% 
  dplyr::select(DATE, CLASS_Z, STAGE, SITE, ENVIRON, PCH4, D13C_CH4, PCO2, gradConcCH4, D13C_CO2) %>%
  mutate(iCH4_mmolL = gradConcCH4*1000) # Excess CH4 in mmol^L-1
#cGas <- cGas[cGas$CLASS_Z == "Surface", ]
cGas <- cGas[!cGas$SITE == "StSe", ] # Omit Stueng Saen tributary site
cGas <- cGas[!cGas$SITE == "KgPr", ] # Omit Kampong Prak tributary site
#cGas$PCH4 <- log10(abs(cGas$PCH4)) 
#cGas <- cGas[!cGas$PCH4 == "-Inf", ]
cGas <- cGas[complete.cases(cGas), ]
#names(cGas)
#View(cGas)
cGas <- cGas %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, PCH4, D13C_CH4, PCO2, D13C_CO2, iCH4_mmolL) %>%
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCH4 = mean(PCH4)) %>%
  mutate(avg.D13C_CH4 = mean(D13C_CH4)) %>% 
  mutate(avg.PCO2 = mean(PCO2)) %>% 
  mutate(avg.D13C_CO2 = mean(D13C_CO2)) %>% 
  mutate(avg.CH4_mmolL = mean(iCH4_mmolL)) %>%
  distinct(avg.PCH4, .keep_all=TRUE) %>%
  ungroup() 
#View(cGas)

plot(cGas$avg.D13C_CO2, cGas$avg.D13C_CH4)
abline(lm(cGas$avg.D13C_CH4~cGas$avg.D13C_CO2))
summary(lm(cGas$avg.D13C_CH4~cGas$avg.D13C_CO2))

query <- cGas %>% 
  dplyr::select(avg.D13C_CH4, avg.D13C_CO2) %>%
  mutate(x = abs(avg.D13C_CH4+50)) %>% 
  mutate(y = abs(avg.D13C_CO2+35))
View(query)

plot(query$x, query$y)
abline(lm(query$y~query$x))
summary(lm(query$y~query$x))

poly_y <- lm(y ~ x + I(x^2), data=query)
summary(poly_y)

# High

Data <- subset(cGas, STAGE=="High") 
Table(cGas$avg.D13C_CH4)
#names(Data)
#View(Data)
sum(abs(cGas$avg.D13C_CO2)>35)
sum(abs(cGas$avg.D13C_CO2)<35)
sum(abs(cGas$avg.D13C_CH4)>50)
sum(abs(cGas$avg.D13C_CH4)<50)
sum((cGas$D13C_CO2 - cGas$D13C_CH4)<31)
sum((cGas$D13C_CO2 - cGas$D13C_CH4)>60)

new_Data <- data.frame(avg.PCH4 = seq(min(Data$avg.PCH4), max(Data$avg.PCH4), length.out = 100)) 
lm_y <- lm(avg.D13C_CO2 ~ avg.PCH4, data=Data)
poly_y <- lm(avg.D13C_CO2 ~ avg.PCH4 + I(avg.PCH4^2), data=Data)
gam_y <- gam(avg.D13C_CO2 ~ s(avg.PCH4, k=5), method="REML", data=Data)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) # Multiple R-squared:  0.3342 # p-value: 0.0002758
pred_Data <- predict(lm_y, newdata=new_Data, se.fit=TRUE, type="response")
exp(-0.5*delAICc1)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) # AICc Weights
exp(-0.5*delAICc2)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) 
exp(-0.5*delAICc3)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) 
new_Data$pred <- predict(lm_y, newdata=new_Data)
#plot(Data$avg.CH4_mmolL, resid(lm_y))

#par(mfrow=c(3, 2))
#par(oma=c(3, 3, 3, 3))
par(mar=c(4, 8, 4, 2)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$avg.PCH4) # 0.0001162435 1.1411083576
range(cGas$avg.D13C_CO2) # -65.762 -11.748
plot(Data$avg.PCH4, Data$avg.D13C_CO2, xlim=c(0, 60000), ylim=c(-85, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")
abline(h=-35, lty=1, lwd=2, col="gray40")
legend(35000, -24, legend=expression("Anaerobic"), text.col=c("gray40"), text.font=4, cex=1.8, bty="n")

polygon(c(seq(min(Data$avg.PCH4), max(Data$avg.PCH4), length.out = 100), rev(seq(min(Data$avg.PCH4), max(Data$avg.PCH4), length.out = 100))), 
        c(pred_Data$fit - (1.96 * pred_Data$se.fit), rev(pred_Data$fit + (1.96 * pred_Data$se.fit))), 
        col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) 
lines(new_Data$avg.PCH4, new_Data$pred, lwd=2, lty=2, col="cornflowerblue") 

points(jitter(Data$avg.PCH4, factor=2), Data$D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
#points(jitter(Edge$avg.PCH4, factor=2), Edge$D13C_CO2, bg=adjustcolor("blue", alpha.f=0.8), col="black", pch=22, cex=3)
#points(jitter(Floodplain$avg.PCH4, factor=2), Floodplain$D13C_CO2, bg=adjustcolor("blue", alpha.f=0.8), col="black", pch=25, cex=3)

mtext(expression(paste(delta^13, C-CO[2], "(\u2030)", sep="")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(italic(P)*CH[4] ~ (mu*atm)), side=1, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-55", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 60000), at=c(0, 15000, 30000, 45000, 60000), 
     labels=c("0", "15000", "30000", "45000", "60000"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1)

par(xpd=TRUE)
legend("topright", legend=c(expression(italic(R^2) ~ "=" ~ "0.33," ~ italic(p) ~ "<" ~ "0.001")), 
       text.col="cornflowerblue", text.font=3,  cex=1.8, bty="n")
legend("topleft", inset=-0.21, legend=expression("c)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

# Falling

Data <- subset(cGas, STAGE=="Falling") 
#names(Data)
#View(Data)
range(Data$avg.D13C_CH4)
sum(abs(Data$avg.D13C_CO2)>35)
sum(abs(Data$avg.D13C_CO2)<35)
sum(abs(Data$avg.D13C_CH4)>50)
sum(abs(Data$avg.D13C_CH4)<50)
sum((Data$avg.D13C_CO2 - Data$avg.D13C_CH4)<31)
sum((Data$avg.D13C_CO2 - Data$avg.D13C_CH4)>31)
sum((cGas$avg.D13C_CO2 - cGas$avg.D13C_CH4)>58)

new_Data <- data.frame(PCH4 = seq(min(Data$PCH4), max(Data$PCH4), length.out = 100)) 
lm_y <- lm(D13C_CO2 ~ PCH4, data=Data)
poly_y <- lm(D13C_CO2 ~ PCH4 + I(PCH4^2), data=Data)
gam_y <- gam(D13C_CO2 ~ s(PCH4, k=5), method="REML", data=Data)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) # Multiple R-squared:  0.02377 # p-value: 0.6324
pred_Data <- predict(lm_y, newdata=new_Data, se.fit=TRUE, type="response")
exp(-0.5*delAICc1)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) # AICc Weights
exp(-0.5*delAICc2)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) 
exp(-0.5*delAICc3)/sum(exp(-0.5*delAICc1), exp(-0.5*delAICc2), exp(-0.5*delAICc3)) 
new_Data$pred <- predict(lm_y, newdata=new_Data)
#plot(Data$PCH4, resid(lm_y))

#par(mfrow=c(3, 2))
#par(oma=c(3, 3, 3, 3))
par(mar=c(4, 3, 4, 7)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$PCH4) #  5.80 8036.83
range(cGas$D13C_CO2) # -65.762 -11.748
plot(Data$avg.PCH4, Data$avg.D13C_CO2, xlim=c(0, 60000), ylim=c(-80, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")
abline(h=-35, lty=1, lwd=2, col="gray40")
legend(35000, -24, legend=expression("Anaerobic"), text.col=c("gray40"), text.font=4, cex=1.8, bty="n")

points(jitter(Data$avg.PCH4, factor=2), Data$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
#points(jitter(Edge$avg.PCH4, factor=2), Edge$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=22, cex=3)

#polygon(c(seq(min(Pelagic$PCH4), max(Pelagic$PCH4), length.out = 100), rev(seq(min(Pelagic$PCH4), max(Pelagic$PCH4), length.out = 100))), 
        #c(pred_Pelagic$fit - (1.96 * pred_Pelagic$se.fit), rev(pred_Pelagic$fit + (1.96 * pred_Pelagic$se.fit))), 
        #col=adjustcolor("gray", alpha.f=0.4), border = NA) # Pelagic
#lines(new_Pelagic$PCH4, new_Pelagic$pred, lwd=4, lty=1, col="black") 

#polygon(c(seq(min(Edge$PCH4), max(Edge$PCH4), length.out = 100), rev(seq(min(Edge$PCH4), max(Edge$PCH4), length.out = 100))), 
        #c(pred_Edge$fit - (1.96 * pred_Edge$se.fit), rev(pred_Edge$fit + (1.96 * pred_Edge$se.fit))), 
        #col=adjustcolor("gray", alpha.f=0.4), border = NA) # Edge
#lines(new_Edge$PCH4, new_Edge$pred, lwd=4, lty=2, col="black") 

#mtext(expression(paste(delta^13, C-CO[2], "(\u2030)", sep="")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(italic(P)*CH[4] ~ (mu*atm)), side=1, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-80, -10), at=c(-10, -25, -40, -65, -80), 
     #labels=c("0", "-20", "-40", "-60", "-80"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 6000), at=c(0, 15000, 30000, 45000, 60000), 
     labels=c("0", "15000", "30000", "45000", "60000"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1)

par(xpd=TRUE)
#legend("topleft", legend=c("Open", "Edge", "Floodplain"), 
       #text.col=c("black"), text.font=c(1, 1, 1),  pch=c(21, 22, 25), col="black", bg="white", pt.cex=c(2.5, 2.5, 2.5), cex=1.75, bty="n")
#legend("topright", inset=-0.26, legend=c("High", "Falling"), 
       #text.col=c("cornflowerblue", "blue3"), text.font=c(1, 1),  col="black", bg="white", pt.cex=c(2.5, 2.5), cex=2.75, bty="n")
legend("topleft", inset=-0.21, legend=expression("d)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


#####################################################
##   Scatterplots for avg.D13C_CH4 & avg.D13C-CO2  ##
#####################################################

#newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>% 
  dplyr::select(DATE, CLASS_Z, STAGE, SITE, ENVIRON, PCH4, D13C_CH4, PCO2, D13C_CO2) 
#mutate(PCH4_x = (1/PCH4)*1000) # Corresponds to 1/ppb of CH4
#cGas <- cGas[cGas$CLASS_Z == "Surface", ]
cGas <- cGas[!cGas$SITE == "StSe", ] # Omit Stueng Saen tributary site
cGas <- cGas[!cGas$SITE == "KgPr", ] # Omit Kampong Prak tributary site
#cGas$PCH4 <- log10(abs(cGas$PCH4)) 
#cGas <- cGas[!cGas$PCH4 == "-Inf", ]
cGas <- cGas[complete.cases(cGas), ]
#names(cGas)
#View(cGas)
cGas <- cGas %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, PCH4, D13C_CH4, PCO2, D13C_CO2) %>%
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCH4 = mean(PCH4)) %>%
  mutate(avg.D13C_CH4 = mean(D13C_CH4)) %>% 
  mutate(avg.PCO2 = mean(PCO2)) %>% 
  mutate(avg.D13C_CO2 = mean(D13C_CO2)) %>% 
  distinct(avg.PCH4, .keep_all=TRUE) %>%
  ungroup() 
View(cGas)

# High

Data <- subset(cGas, STAGE=="High") 
range(Data$D13C_CO2)
range(Data$D13C_CH4)
#names(Data)
#View(Data)

#par(mfrow=c(3, 2))
#par(oma=c(3, 3, 3, 3))
par(mar=c(4, 8, 4, 2)) # Bottom, left, top, right par(mar=c(4, 2, 4, 8)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$D13C_CH4) # -54.325 -14.747
range(cGas$D13C_CO2) # -58.843 -29.194
plot(Data$avg.D13C_CH4, Data$avg.D13C_CO2, xlim=c(-85, -10), ylim=c(-85, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")

abline(5, 1, lty=2, lwd=1, col="black")
legend(-30, -15, legend=expression(epsilon[C] ~ "=" ~ "5"), text.col=c("black"), text.font=4, cex=1.8, bty="n")
abline(20, 1, lty=2, lwd=1, col="black")
#legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(30, 1, lty=2, lwd=1, col="black")
legend(-55, -15, legend=expression(epsilon[C] ~ "=" ~ "30"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(40, 1, lty=2, lwd=1, col="black")
legend(-70, -20, legend=expression(epsilon[C] ~ "=" ~ "40"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(55, 1, lty=2, lwd=1, col="black")
legend(-80, -15, legend=expression(epsilon[C] ~ "=" ~ "55"), text.col=c("black"), text.font=4, cex=1.75, bty="n")

polygon(x=c(-72, -50, -50, -72), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray", alpha.f=0.6), border=NA)
#text(-59, -16, expression(italic("Acetate \n Ferm.")), cex=1.55, font=3)

polygon(x=c(-100, -69.5, -69.5, -100), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray", alpha.f=0.6), border=NA)
text(-78, -16, expression(italic("Carbonate \n Red.")), cex=1.55, font=3)

polygon(x=c(-65, -10, -10, -65), y=c(-35, -35, -85, -85), col=adjustcolor("gray", alpha.f=0.6), border=NA)
text(-35, -57, expression(italic("Oxidation")), cex=1.55, font=3)

points(jitter(Data$avg.D13C_CH4, factor=4), Data$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
#points(jitter(Edge$avg.D13C_CH4, factor=4), Edge$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=22, cex=3)
#points(jitter(Floodplain$avg.D13C_CH4, factor=4), Floodplain$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=25, cex=3)
points(-48, -10, col="darkorange", pch=19, cex=1.5)
legend(-54, -7, legend=expression("Atm." ~ CO[2] ~ "&" ~ CH[4]), text.col=c("darkorange"), text.font=4, cex=1.8, bty="n")
text(-59, -16, expression(italic("Acetate \n Ferm.")), cex=1.55, font=3)
legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
mtext(expression(paste(delta^13, C-CO[2], "(\u2030)", sep="")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13, C-CH[4], "(\u2030)", sep="")), side=1, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-55", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-65", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
par(xpd=TRUE)
legend("topleft", inset=-0.21, legend=expression("e)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

# Falling

Data <- subset(cGas, STAGE=="Falling") 
#names(Data)
#View(Data)

#par(mfrow=c(3, 2))
#par(oma=c(3, 3, 3, 3))
par(mar=c(4, 3, 4, 7)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$D13C_CH4) # -54.325 -14.747
range(cGas$D13C_CO2) # -58.843 -29.194
plot(Data$avg.D13C_CH4, Data$avg.D13C_CO2, xlim=c(-85, -10), ylim=c(-85, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")

abline(5, 1, lty=2, lwd=1, col="black")
legend(-30, -15, legend=expression(epsilon[C] ~ "=" ~ "5"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(20, 1, lty=2, lwd=1, col="black")
legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(30, 1, lty=2, lwd=1, col="black")
legend(-55, -15, legend=expression(epsilon[C] ~ "=" ~ "30"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(40, 1, lty=2, lwd=1, col="black")
legend(-70, -20, legend=expression(epsilon[C] ~ "=" ~ "40"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(55, 1, lty=2, lwd=1, col="black")
legend(-80, -15, legend=expression(epsilon[C] ~ "=" ~ "55"), text.col=c("black"), text.font=4, cex=1.75, bty="n")

polygon(x=c(-72, -50, -50, -72), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray", alpha.f=0.6), border=NA)
text(-59, -16, expression(italic("Acetate \n Ferm.")), cex=1.55, font=3)

polygon(x=c(-100, -69.5, -69.5, -100), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray", alpha.f=0.6), border=NA)
text(-78, -16, expression(italic("Carbonate \n Red.")), cex=1.55, font=3)

polygon(x=c(-65, -10, -10, -65), y=c(-35, -35, -85, -85), col=adjustcolor("gray", alpha.f=0.6), border=NA)
text(-35, -57, expression(italic("Oxidation")), cex=1.55, font=3)

points(jitter(Data$avg.D13C_CH4, factor=2), Data$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
#points(jitter(Edge$avg.D13C_CH4, factor=2), Edge$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=22, cex=3)
#points(jitter(Floodplain$D13C_CH4, factor=2), Floodplain$D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=25, cex=3)
points(-48, -10, col="darkorange", pch=19, cex=1.5)
legend(-54, -7, legend=expression("Atm." ~ CO[2] ~ "&" ~ CH[4]), text.col=c("darkorange"), text.font=4, cex=1.8, bty="n")

#mtext(expression(paste(delta^13, C-CO[2], "(\u2030 PDB)", sep="")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13, C-CH[4], "(\u2030)", sep="")), side=1, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-60, -10), at=c(-10, -20, -30, -40, -50, -60), 
#labels=c("-10", "-20", "-30", "-40", "-50", "-60"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-65", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
par(xpd=TRUE)
legend("topleft", inset=-0.21, legend=expression("f)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

