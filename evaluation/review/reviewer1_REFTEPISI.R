setwd(file.path(path.expand(dirname("~")), "Downloads"))
library(tidyverse)
library(lme4)
library(car)
library(zoo)
library(ggplot2)
library(qqplotr)
library(patchwork)
library(interactions)
library(effects)
library(MuMIn) # for effect size
library(sjPlot)
library(dplyr)
library(performance)
library(lmtest)
library(astsa)

# From the MoCsEFC dataset:
minISI.MoCsEFC <- 2.246


dataset <- read.delim("REFTEP++_ISI__2023-07-04T122658.txt", sep=",")
dataset <- dataset[!is.infinite(dataset$ISI),]
trials.kept <- table(dataset$Subject)
dataset <- dataset[dataset$PreRangeAPB < 50 & dataset$PreRangeFDI < 50,]
trials.kept <- table(dataset$Subject) / trials.kept
summary(dataset)

dataset$PowerC3H <- log(dataset$PowerC3H)
dataset$PowerC4H <- log(dataset$PowerC4H)
dataset$CosPhaseC3H <- cos(dataset$PhaseC3H)
dataset$SinPhaseC3H <- sin(dataset$PhaseC3H)
dataset$CosPhaseC4H <- cos(dataset$PhaseC4H)
dataset$SinPhaseC4H <- sin(dataset$PhaseC4H)

#dataset$ISI <- dataset$ISI - mean(dataset$ISI)
dataset$ISI <- sqrt(sqrt(dataset$ISI - minISI.MoCsEFC))

dataset$Condition <- "mid"
dataset$SufficientMuSNR <- NA
dataset$SessionProgress <- NaN

good.mu <- c(T,T,T,T,T,F,F,F,T,T,T,F,T,T,T,T,T,T,T)
names(good.mu) <- c("REFTEP_018", "REFTEP_019", "REFTEP_020", "REFTEP_021", "REFTEP_022", "REFTEP_023", "REFTEP_024", "REFTEP_025", "REFTEP_026", "REFTEP_027", "REFTEP_028", "REFTEP_029", "REFTEP_030", "REFTEP_031(141)", "REFTEP_032", "REFTEP_032_2", "REFTEP_034", "REFTEP_035(198)", "REFTEP_036")


for (ID in unique(dataset$Subject)) {
  for (ses in 1:4) {
    mask <- dataset$Subject == ID & dataset$Block == ses
    x <- dataset$ResponseFDI[mask]
    # 1) Nonlinear Transform
    x <- sqrt(sqrt(x))
    
    # 2) Antimedian filter
    window.radius <- 150 
    if(length(x) > 2*window.radius+1) {
      padded <- c(rev(x[1:window.radius]), x, x[length(x)-(1:window.radius)+1])
      x <- x - zoo::rollmedian(padded, k=2*window.radius + 1)  
    } else {
      x <- x - median(x)
    }
    x <- detrend(x, order=1)
    
    # 3) z-score
    x <- (x - mean(x)) / sd(x)
    
    dataset$ResponseFDI[mask] <- x
    dataset$SessionProgress[mask] <- (1:length(x))/length(x)
    # Well of course, if i split the response!
    dataset$Condition[mask] <- ifelse(dataset$PLV[mask] < quantile(dataset$PLV[mask],0.25), "low", ifelse(dataset$PLV[mask] > quantile(dataset$PLV[mask],0.75), "high", "mid"))
    
    # Try out within-session normalization of Power:
    pL <- dataset$PowerC3H[mask]
    pR <- dataset$PowerC4H[mask]
    pL <- detrend(dataset$PowerC3H[mask], order=1)
    pR <- detrend(dataset$PowerC4H[mask], order=1)
    dataset$PowerC3H[mask] <- (pL - mean(pL)) / sd(pL)
    dataset$PowerC4H[mask] <- (pR - mean(pR)) / sd(pR)
    dataset$SufficientMuSNR[mask] <- good.mu[ID]
  }
}

#ggplot(data=dataset, aes(x=ResponseFDI, color=Subject)) + geom_density()
#ggplot(data=dataset, aes(x=PowerC3H, color=Subject)) + geom_density()
#ggplot(data=dataset, aes(x=PowerC4H, color=Subject)) + geom_density()
#ggplot(data=dataset, aes(x=PLV, color=Condition)) + geom_density()
#ggplot(data=dataset, aes(x=ResponseFDI, color=Condition)) + geom_density()

#ggplot(data=dataset, aes(x=SessionProgress, y=PowerC3H, color=Subject)) + geom_point()

sort(table(dataset$Subject))


ggplot(data=dataset, aes(x=ISI, color=Subject)) + geom_density()
m.C.REFTEP.ISI <- lmer(ResponseFDI ~ Condition*ISI*PowerC3H*PowerC4H + (1+Condition|Subject), data=dataset)
m.REFTEP.ISI <- lmer(ResponseFDI ~ PLV*ISI*PowerC3H*PowerC4H + (1+PLV|Subject), data=dataset) 
car::Anova(m.REFTEP.ISI)
m0.REFTEP.ISI <- lmer(ResponseFDI ~ (PLV*ISI*PowerC3H*PowerC4H - ISI:PowerC3H) + (1+PLV|Subject), data=dataset)
lrtest(m0.REFTEP.ISI, m.REFTEP.ISI)

m.i.REFTEP.ISI <- lmer(ResponseFDI ~ iPLV*ISI*PowerC3H*PowerC4H + (1+iPLV|Subject), data=dataset)


m.smaller.REFTEP.ISI <- lmer(ResponseFDI ~ PLV + ISI + PowerC3H + PowerC4H
                             + PLV:ISI + PLV:PowerC3H + PLV:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H + (1+PLV|Subject), data=dataset)

m0.smaller.REFTEP.ISI <- lmer(ResponseFDI ~ PLV + ISI + PowerC3H + PowerC4H
                             + PLV:PowerC3H + PLV:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H + (1+PLV|Subject), data=dataset)
m0b.smaller.REFTEP.ISI <- lmer(ResponseFDI ~ PLV + ISI + PowerC3H + PowerC4H
                              + PLV:ISI + PLV:PowerC3H + PLV:PowerC4H + PowerC3H:PowerC4H + (1+PLV|Subject), data=dataset)
summary(m.smaller.REFTEP.ISI)

lrtest(m0.smaller.REFTEP.ISI, m.smaller.REFTEP.ISI) # PLV:ISI is not supported by LR

# 50 for the ISI:Power
lrtest(m0b.smaller.REFTEP.ISI, m.smaller.REFTEP.ISI) # ISI:C3Power is not supported by LR









interact_plot(m.REFTEP.ISI, pred=ISI, modx=PowerC4H,plot.points=T)


m.ph.p.3 <- lm(ResponseFDI ~ PowerC3H*(CosPhaseC3H + SinPhaseC3H), data=dataset)
m.ph.p.2 <- lm(ResponseFDI ~ PowerC3H + CosPhaseC3H + SinPhaseC3H, data=dataset)
m.ph.p.1 <- lm(ResponseFDI ~ PowerC3H, data=dataset)

lrtest(m.ph.p.1, m.ph.p.2) # We find that phase is predictive of MEP,
lrtest(m.ph.p.2, m.ph.p.3) # But no phase-power interaction


dataset.only.good.SNR <- dataset[dataset$SufficientMuSNR,]
m.REFTEP.ISI.SNR <- lmer(ResponseFDI ~ PLV*ISI*PowerC3H*PowerC4H + (1+PLV|Subject), data=dataset.only.good.SNR)
m.smaller.REFTEP.ISI.SNR <- lmer(ResponseFDI ~ PLV + ISI + PowerC3H + PowerC4H
                             + PLV:ISI + PLV:PowerC3H + PLV:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H + (1+PLV|Subject), data=dataset.only.good.SNR)


# Weird results:
m.REFTEP.sessionProgress <- lmer(ResponseFDI ~ PLV*ISI*PowerC3H*PowerC4H*SessionProgress + (1+PLV|Subject), data=dataset)
car::Anova(m.REFTEP.sessionProgress)
