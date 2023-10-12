# Run check-volume-conduction.R and LME_in_R.R first
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
library(lmerTest)


anova_file <- file("p-values.log", "a")

dataset.raw <- read.delim("results___2023-06-13T171710.txt", sep=",")
# remove the excess trial of subj10
dataset.raw <- dataset.raw[setdiff(1:nrow(dataset.raw), 1800*(10-1)+900+621),]

dataset.raw <- dataset.raw[dataset.raw$Subject != 8,] # due to the patched-together nature, 
# correspondence between the two dataframes is not clear.



dataset.realtime.aug <- read.delim("realtime-cPLV.csv", sep=",")
dataset.realtime.aug <- dataset.realtime.aug[dataset.realtime.aug$Subject != 8,]


dataset <- data.frame(dataset.raw)


dataset$MuPowerR <- log(dataset.raw$MuPowerR)
dataset$MuPowerL <- log(dataset.raw$MuPowerL)

# WARNING: THIS SHIFT WAS NOT DONE WHEN ORIGINALLY SUBMITTING! -> Correct corresponding p-Values!
dataset$Waittime <- dataset.raw$Waittime - min(dataset.raw$Waittime) # remove min ISI
dataset$Waittime <- sqrt(sqrt(dataset$Waittime))
dataset$Waittime[dataset.raw$Waittime == 10] <- NaN

# Make Condition clearly categorical
dataset["Condition"][dataset["Condition"] == 0] <- "low"
dataset["Condition"][dataset["Condition"] == 1] <- "high"

dataset <- dataset %>% rename(ISI = Waittime)

# Checking that the conditions line up:
# dataset[which(dataset.realtime.aug$Condition != dataset$Condition),]
# Only subject 8 is not lining up -- but we'll exclude that one later!

dataset$RealPart <- dataset.realtime.aug$RealPart
dataset$ImagPart <- dataset.realtime.aug$ImagPart
dataset$PLV <- dataset.realtime.aug$PLV
dataset$iPLV <- dataset.realtime.aug$iPLV

# Throw out trials with too much preinnervation (BOTH muscles must LACK preinnervation)
dataset <- dataset[dataset$PreRangeAPB < 50 & dataset$PreRangeFDI < 50,]
# also do this BEFORE median-window, to remove the outliers (?)

for (ID in unique(dataset$Subject)) {
  for (ses in 1:2) {
    x <- dataset$ResponseFDI[dataset$Subject == ID & dataset$Session == ses]
    # 1) Nonlinear Transform
    x <- sqrt(sqrt(x))
    
    # 2) Antimedian filter
    window.radius <- 150
    #padded <- c(rep(x[1], window.radius), x, rep(x[length(x)], window.radius))
    padded <- c(rev(x[1:window.radius]), x, x[length(x)-(1:window.radius)+1])
    x <- x - zoo::rollmedian(padded, k=2*window.radius + 1)
    
    # 3) z-score
    x <- (x - mean(x)) / sd(x) # Pedro says yes
    # Really z-Score? Or instead do (1|Subject/Session)
    
    dataset$ResponseFDI[dataset$Subject == ID & dataset$Session == ses] <- x
  }
  
  # Make Subject explicitly categorical
  dataset$Subject[dataset$Subject == ID] <- sprintf("%02d", ID) #format(ID)
}





# Check percentage of trials removed due to pre-innervation:
trials.kept <- table(dataset$Subject) / table(dataset.raw$Subject)

# Make Session clearly categorical
dataset["Session"][dataset["Session"] == 1] <- "1"
dataset["Session"][dataset["Session"] == 2] <- "2"



# PLV across conditions: We observe that in spite of the moving individual 
# thresholds, the conditions are well-distinct (although overlapping)
#distr.PLV.cross.conditions <- ggplot(data=dataset, aes(x=PLV, group=Condition, color=Condition, fill=Condition)) + geom_density(alpha=0.5) + theme_bw() 




# Models

lmer.plv <- lmer(ResponseFDI ~ 1 + PLV*MuPowerR*MuPowerL*ISI + (1+Condition|Subject), 
                  data = dataset, REML=FALSE)
car::Anova(lmer.plv)

lmer.iplv <- lmer(ResponseFDI ~ 1 + iPLV*MuPowerR*MuPowerL*ISI + (1+Condition|Subject), 
                 data = dataset, REML=FALSE)
car::Anova(lmer.iplv)

# Checking interactions: Really, the same as for Condition.
# interact_plot(lmer.plv, pred = MuPowerL, modx = PLV, interval = TRUE, int.width = 0.8, plot.points=TRUE, colors = "Set1")
# interact_plot(lmer.plv, pred = ISI, modx = PLV, interval = TRUE, int.width = 0.8, plot.points=TRUE, colors = "Set1")
# interact_plot(lmer.iplv, pred = MuPowerL, modx = iPLV, interval = TRUE, int.width = 0.8, plot.points=TRUE, colors = "Set1")
# interact_plot(lmer.iplv, pred = ISI, modx = iPLV, interval = TRUE, int.width = 0.8, plot.points=TRUE, colors = "Set1")


sink(anova_file, append=T)
cat("\n\n\n\n\n")
cat(rep('#', 50), sep="#", end="\n")
cat('\n Posthoc Realtime \n\n')
lmer.plv.red <- lmer(ResponseFDI ~ 1 + PLV+MuPowerR+MuPowerL
                     + PLV:MuPowerR + PLV:MuPowerL + MuPowerR:MuPowerL
                     + PLV:ISI #+ PLV:MuPowerR:MuPowerL # not helpful
                     + (1+Condition|Subject), 
                 data = dataset, REML=FALSE)
car::Anova(lmer.plv.red)

lmer.iplv.red <- lmer(ResponseFDI ~ 1 + iPLV+MuPowerR+MuPowerL
                      + iPLV:MuPowerR + iPLV:MuPowerL + MuPowerR:MuPowerL
                      + iPLV:ISI
                      + (1+Condition|Subject), 
                      data = dataset, REML=FALSE)
car::Anova(lmer.iplv.red)
sink()
close(anova_file)










# RUN check_phase.R before this:
lmer.best.PLV   <- lmer(ResponseFDI ~ 1 + PLV + MuPowerR + MuPowerL
                        + PLV:MuPowerL + MuPowerR:MuPowerL + PLV:ISI 
                        + (1+PLV|Subject),
                        data=dataset, REML=FALSE)
lmer.best.PLV.w.phase <- lmer(ResponseFDI ~ 1 + PLV + MuPowerR + MuPowerL
                              + PLV:MuPowerL + MuPowerR:MuPowerL + PLV:ISI + CosPhaseC3 + SinPhaseC3
                              + (1+PLV|Subject),
                              data=dataset, REML=FALSE)
lmer.best.PLV.w.phase.power <- lmer(ResponseFDI ~ 1 + PLV + MuPowerR + MuPowerL
                                    + PLV:MuPowerL + MuPowerR:MuPowerL + PLV:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL
                                    + (1+PLV|Subject),
                                    data=dataset, REML=FALSE)
lmer.best.PLV.w.all.interactions <- lmer(ResponseFDI ~ 1 + PLV + MuPowerR + MuPowerL
                                         + PLV:MuPowerL + MuPowerR:MuPowerL + PLV:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL + (CosPhaseC3 + SinPhaseC3):Condition
                                         + (1+PLV|Subject),
                                         data=dataset, REML=FALSE)
cat(rep('_', 50), sep="_", end="\n")
anova(lmer.best.PLV.w.phase, lmer.best.PLV)
anova(lmer.best.PLV.w.phase.power, lmer.best.PLV.w.phase)
anova(lmer.best.PLV.w.all.interactions, lmer.best.PLV.w.phase.power)

lmer.best.iPLV   <- lmer(ResponseFDI ~ 1 + iPLV + MuPowerR + MuPowerL
                        + iPLV:MuPowerL + MuPowerR:MuPowerL + iPLV:ISI 
                        + (1+iPLV|Subject),
                        data=dataset, REML=FALSE)
lmer.best.iPLV.w.phase <- lmer(ResponseFDI ~ 1 + iPLV + MuPowerR + MuPowerL
                              + iPLV:MuPowerL + MuPowerR:MuPowerL + iPLV:ISI + CosPhaseC3 + SinPhaseC3
                              + (1+iPLV|Subject),
                              data=dataset, REML=FALSE)
lmer.best.iPLV.w.phase.power <- lmer(ResponseFDI ~ 1 + iPLV + MuPowerR + MuPowerL
                                    + iPLV:MuPowerL + MuPowerR:MuPowerL + iPLV:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL
                                    + (1+iPLV|Subject),
                                    data=dataset, REML=FALSE)
lmer.best.iPLV.w.all.interactions <- lmer(ResponseFDI ~ 1 + iPLV + MuPowerR + MuPowerL
                                         + iPLV:MuPowerL + MuPowerR:MuPowerL + iPLV:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL + (CosPhaseC3 + SinPhaseC3):Condition
                                         + (1+iPLV|Subject),
                                         data=dataset, REML=FALSE)
cat(rep('_', 50), sep="_", end="\n")
anova(lmer.best.iPLV.w.phase, lmer.best.iPLV)
anova(lmer.best.iPLV.w.phase.power, lmer.best.iPLV.w.phase)
anova(lmer.best.iPLV.w.all.interactions, lmer.best.iPLV.w.phase.power)
