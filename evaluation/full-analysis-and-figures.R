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



WRITE.PLOTS <- T
anova_file <- file("p-values.log", "w")


# 2023-03-15T165605

# New phases and power:
dataset.raw <- read.delim("results___2023-06-13T171710.txt", sep=",")


#dataset.realtime.aug <- read.delim("realtime-cPLV.csv", sep=",")

dataset <- data.frame(dataset.raw)
dataset$MuPowerR <- log(dataset.raw$MuPowerR)
dataset$MuPowerL <- log(dataset.raw$MuPowerL)

#cutoff.ISI <- 2.35 # 2.35 for first group, 2.46 for second group
dataset$Waittime <- dataset.raw$Waittime - min(dataset.raw$Waittime) # remove min ISI
dataset$Waittime[dataset.raw$Waittime == 10] <- NaN
#dataset <- dataset[dataset.raw$Waittime >= cutoff.ISI,] # NEW! May omit short ISIs
dataset$Waittime <- sqrt(sqrt(dataset$Waittime))


# Consideration: We could omit the very short ISIs as "artifact-affected".
# This reduces our sample size by about 1/4, and we still in fact find the
# interaction Condition:ISI, meaning this is not *purely* due to the difference 
# between short and long ISIs


# Why do we remove the min ISI?
# -> The raw ISIs are >2 (~ > 1.189 respectively after fourth root)
# -> The main effect of e.g. Condition will however be given at ISI=0 -- leading
#    to extrapolation (inflation of the estimated effect)
# If we instead do (ISI_i - minISI)^0.25, we can still go back from the derived value quite easily
# and the estimated main effect becomes the effect at minISI


# Make Condition clearly categorical
dataset["Condition"][dataset["Condition"] == 0] <- "low"
dataset["Condition"][dataset["Condition"] == 1] <- "high"

dataset <- dataset %>% rename(ISI = Waittime)

# Throw out trials with too much preinnervation (BOTH muscles must LACK preinnervation)
dataset <- dataset[dataset$PreRangeAPB < 50 & dataset$PreRangeFDI < 50,]
# also do this BEFORE median-window, to remove the outliers

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





##### L M E Models:
sink(anova_file,append=F)
cat('Naive model:\n')
naive.lmer <- lmer(ResponseFDI ~ 1+ Condition + (1+Condition|Subject), 
                   data = dataset, REML=FALSE)

car::Anova(naive.lmer)
naive.effect.size <- r.squaredGLMM(naive.lmer)[1] # R2m 


cat('\n______________________________\nMixed (Power) model:\n')
mixed.lmer <- lmer(ResponseFDI ~ 1+ Condition*MuPowerR*MuPowerL + (1+Condition|Subject), # 
                   data = dataset, REML=FALSE)
car::Anova(mixed.lmer)
mixed.effect.size <- r.squaredGLMM(mixed.lmer)[1] # R2m 

cat('\n______________________________\nMaximal (Power,ISI) model:\n')
lmer.time <- lmer(ResponseFDI ~ 1 + Condition*MuPowerR*MuPowerL*ISI + (1+Condition|Subject), 
                  data = dataset, REML=FALSE)
car::Anova(lmer.time)

cat('\n______________________________\nSummary model:\n')
lmer.best   <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL + ISI
                    + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI 
                    + Condition:MuPowerR:MuPowerL
                    + (1+Condition|Subject), #+ (1+Condition|Subject) + (1|Session), # The results are a bit stronger with this, but very likely overfitting
                    data=dataset, REML=FALSE)
car::Anova(lmer.best)
summary(lmer.best)

cat('\n______________________________\nAkaike Information Criteria:\n')
AIC(naive.lmer, mixed.lmer, lmer.time, lmer.best)
sink()
close(anova_file)
# Continuous plv/iplv results are in posthoc-realtime-prediction


gdata <- dataset %>% group_by(Subject, Condition) %>% summarise(SD = sd(ResponseFDI), ResponseFDI = mean(ResponseFDI), n=n())

df.fx <- as.data.frame(effect('Condition', naive.lmer))
df.fx$Group = "all"

(gt <- ggplot(data=df.fx, aes(x=Condition, y=fit, group=Group)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .4, colour = NA, fill="dodgerblue") + 
    geom_line(lwd=2, color="dodgerblue4") +
    geom_point(size=4, color="dodgerblue4") + 
    geom_line(data=gdata, aes(x = Condition, y = ResponseFDI, group=Subject), alpha=0.4) + 
    geom_point(data=gdata, aes(x = Condition, y = ResponseFDI, group=Subject), alpha=0.4, size=1.5) + 
    theme_bw())
if(WRITE.PLOTS) {
  ggsave(
    plot = gt,
    filename = "mainResult.pdf",
    bg = "transparent", width=1600, height=1300,units="px"
  )
}



##### Power-overview

p.L <- ggplot(data=dataset, aes(x=MuPowerL, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("log \u03bc power (left)") +
  theme_bw()


p.R <- ggplot(data=dataset, aes(x=MuPowerR, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("log \u03bc power (right)") +
  theme_bw()


plot.r.pow.resp <- ggplot(dataset, aes(x=MuPowerR, y=ResponseFDI) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", option="A") +
  xlab("log \u03bc power (right)") +
  ylab("ResponseFDI (z-score)") +
  theme_bw()


plot.l.pow.resp <- ggplot(dataset, aes(x=MuPowerL, y=ResponseFDI) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", option="A") +
  xlab("log \u03bc power (left)") +
  ylab("ResponseFDI (z-score)") +
  theme_bw()

p.R <- p.R + theme(legend.position="none")
if(WRITE.PLOTS) {
  ggsave(
    plot = plot.l.pow.resp + plot.r.pow.resp + p.L + p.R,
    filename = "power-data.pdf",
    bg = "transparent", width=2100, height=1500,units="px"
  )
}





##### Mixed and final LMER-plots

sjRes.mixed <- sjPlot::plot_model(mixed.lmer,show.p = TRUE,show.value=TRUE,
                             vline.color="gray90",value.offset=0.3,
                             dot.size=1,
                             colors="black") + 
    ylim(-0.3,0.3) +
    labs(title=element_blank(), y="Estimated fixed effect on ResponseFDI") +
    scale_x_discrete(labels=rev(rownames(as.data.frame(car::Anova(mixed.lmer))))) +
    theme_bw() + 
    theme(axis.text=element_text(size=10), aspect.ratio = 1)
# correct p-values to Wald chi2 type II
p <- as.data.frame(car::Anova(mixed.lmer))$`Pr(>Chisq)`
star.labels <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fe <- round(fixef(mixed.lmer), digits=3)
sjRes.mixed$data$p.stars <- star.labels
sjRes.mixed$data$p.value <- p
sjRes.mixed$data$p.label <- sprintf("%g %s", fe[2:length(fe)], star.labels)

ipl.mixed <- interact_plot(mixed.lmer, pred = MuPowerL, modx = Condition, 
                           interval = TRUE, int.width = 0.8, 
                           plot.points=TRUE, colors = "Set1") + 
    labs(title=element_blank(), x="log \u03bc power (left)", y="ResponseFDI (z-score)") +
    theme_bw() + theme(aspect.ratio = 1)
ipl.mixed$layers[[3]]$aes_params$alpha = 0.02
#ipl$theme$panel.background$colour = "black"





sjRes.best <- sjPlot::plot_model(lmer.best, show.p = TRUE,show.value=TRUE, 
                                 vline.color = "gray90",value.offset=0.3,
                                 p.val="wald",
                                 dot.size=1, colors="black")  + 
    labs(title=element_blank(), y="Estimated fixed effect on ResponseFDI") +
    scale_x_discrete(labels=rev(rownames(as.data.frame(car::Anova(lmer.best))))) +
    theme_bw() + 
    theme(axis.text=element_text(size=10), aspect.ratio = 1)
# correct p-values to Wald chi2 type II
p <- as.data.frame(car::Anova(lmer.best))$`Pr(>Chisq)`
star.labels <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fe <- round(fixef(lmer.best), digits=3)
sjRes.best$data$p.stars <- c(star.labels)#, star.labels[length(star.labels)]) # Conditions high and low should have the same p-value -- it is not an option to "just add one"
sjRes.best$data$p.label <- sprintf("%g %s", fe[2:length(fe)], c(star.labels))#, star.labels[length(star.labels)]))
sjRes.best$data$p.value <- p #c(p, p[length(p)])

ipl.best <- interact_plot(lmer.best, pred = ISI, modx = Condition, 
                          interval = TRUE, int.width = 0.8, 
                          plot.points=TRUE, colors = "Set1", point.alpha=0.02) + 
    labs(x = expression(sqrt(sqrt(ISI - min))), y="ResponseFDI (z-score)") +
    theme_bw() + theme(aspect.ratio = 1)
ipl.best$layers[[3]]$aes_params$alpha = 0.02

#(ipls <- interact_plot(lmer.best, pred=MuPowerL, modx=Condition, mod2=MuPowerR, colors="Set1", interval=T) + theme(aspect.ratio=1))





# POSTER:
library(ggthemes)
ipls <- interact_plot(lmer.best, pred=MuPowerL, modx=Condition, mod2=MuPowerR, colors="Set1", interval=T, vary.lty=F)
ipls$data$MuPowerR <- factor(ipls$data$MuPowerR, levels=sort(unique(ipls$data$MuPowerR), decreasing=T))
#levels(ipls$data$MuPowerR) <- c("High MuPowerR", "Mean MuPowerR", "Low MuPowerR")
(ipls <- ipls + facet_wrap(~ MuPowerR, ncol=1, dir="v", labeller=as_labeller(c(`-1.17274552030741`='Mean MuPowerR - 1 SD', `0.385373869056703`='Mean MuPowerR', `1.94349325842081`='Mean MuPowerR + 1 SD')))
  + theme_base() + theme(aspect.ratio=1))
# These are those magic values: unique(ipls$data$MuPowerR), but they are for some reason printed with 14 digits in interact_plot

ipl.ISI.poster <- interact_plot(lmer.best, pred = ISI, modx = Condition, 
                          interval = TRUE, int.width = 0.8, 
                          plot.points=T, colors = "Set1", point.alpha=0.02, vary.lty=F) + 
    theme_base() + labs(x = expression(sqrt(sqrt(ISI - min))), y="ResponseFDI (z-score)") + theme(aspect.ratio = 1)
ipl.ISI.poster


if(WRITE.PLOTS) {
  ggsave(
    plot = sjRes.mixed+ipl.mixed,
    filename = "mixedResult.pdf",
    bg = "transparent", width=3000, height=1100,units="px"
  )
  ggsave(
    plot = ipl.mixed,
    filename = "interaction-power-FC.pdf",
    bg = "transparent", width=1200, height=1100,units="px"
  )
  ggsave(
    plot = sjRes.best+ipl.best,
    filename = "bestModelOverview.pdf", # -cropped-ISI
    bg = "transparent", width=3000, height=1100,units="px"
  )
  ggsave(
    plot = ipl.best,
    filename = "interaction-FC-ISI.pdf",
    bg = "transparent", width=1200, height=1100,units="px"
  )
  
  ggsave(
    plot = (sjRes.mixed | ipl.mixed) / (sjRes.best | ipl.best),
    filename = "nonnaiveResults.pdf",
    bg = "transparent", width=3000, height=2200,units="px"
  )
  
  ggsave(
    plot = ipls,
    filename = "top_level_interaction.pdf",
    bg = "transparent", width=1400, height=2400,units="px"
  )
  ggsave(
    plot = ipl.ISI.poster,
    filename = "poster-isi.pdf",
    bg = "transparent", width=1400, height=1090,units="px"
  )
  ggsave(
    plot = ipl.ISI.poster,
    filename = "poster-isi.png",
    bg = "transparent", width=1400, height=1090,units="px"
  )
}






##### 
# Individual results:
condition.p.values <- c();
t.p.values <- c();
t.p.1 <- c();
t.p.2 <- c();

for(i in unique(dataset$Subject)) {
  dataset.ind <- data.frame(dataset)
  dataset.ind <- dataset.ind[dataset.ind$Subject == i,]
  lmer.s1 <- lmer(ResponseFDI ~ 1 + Condition + (1|Session), data=dataset.ind)
  #lmer.s1 <- lmer(ResponseFDI ~ 1 + Condition*MuPowerR*MuPowerL + (1|Session), data=dataset.ind)
  cat("\n\n\nIndividual result for participant ", i, "\n\n")
  res.anova <- car::Anova(lmer.s1)
  print(res.anova)
  condition.p.values <- append(condition.p.values, res.anova$`Pr(>Chisq)`[1])
  
  
  dataset.ind.1 <- data.frame(dataset.ind)
  dataset.ind.1 <- dataset.ind.1[dataset.ind.1$Session == 1,]
  t.res <- wilcox.test(dataset.ind.1$ResponseFDI ~ dataset.ind.1$Condition, alternative="greater", paired=FALSE)
  t.p.1 <- append(t.p.1, t.res$p.value)
  cat("\n            t-test (1st sessions):   ",t.res$p.value)
  
  dataset.ind.2 <- data.frame(dataset.ind)
  dataset.ind.2 <- dataset.ind.2[dataset.ind.2$Session == 2,]
  t.res <- wilcox.test(dataset.ind.2$ResponseFDI ~ dataset.ind.2$Condition, alternative="greater", paired=FALSE)
  t.p.2 <- append(t.p.2, t.res$p.value)
  cat("\n            t-test (2nd sessions):   ",t.res$p.value, "\n\n")
  
  #x1 <- dataset.ind$ResponseFDI[dataset.ind$Session == 1]
  #dataset.ind$ResponseFDI[dataset.ind$Session == 1] <- (x1 - mean(x1)) / sd(x1)
  #x2 <- dataset.ind$ResponseFDI[dataset.ind$Session == 2]
  #dataset.ind$ResponseFDI[dataset.ind$Session == 2] <- (x2 - mean(x2)) / sd(x2)
  
  t.res <- wilcox.test(dataset.ind$ResponseFDI ~ dataset.ind$Condition, alternative="greater", paired=FALSE)
  cat("\nt-test (one-sided, both sessions):   ",t.res$p.value)
  t.p.values <- append(t.p.values, t.res$p.value)
  
}


# Format given dataframe as Latex table
format_as_table <- function(df, digits=4) {
  nCol <- length(colnames(df))
  
  header <- paste("\\begin{tabular}{", paste(rep("l", nCol), collapse=""), "}", sep ="")
  foot <- "\\end{tabular}"
  
  body <- paste("\t", paste(colnames(df), collapse=" & "), "\\\\\\hline\n\t")
  for(row in 1:nrow(df)) {
    body <- paste(body, paste(format(round(df[row,], digits), scientific=FALSE), collapse=" & "), sep="\\\\\n\t")
  }
  
  return(paste(header, body, foot, sep="\n"))
}

cat(format_as_table(data.frame(S1=t.p.1, S2=t.p.2, both=t.p.values)))





####
# Checking what happens when we exclude the two participants with the strongest 
# bias to 0deg phase-shift:
anova_file <- file("p-values.log", "a")
sink(anova_file,append=T)
cat('\n Best with 03 and 14 removed (volume conduction)\n')
dataset.reduced <- data.frame(dataset)
dataset.reduced <- dataset.reduced %>% filter((Subject != "03") & (Subject != "14"))

lmer.best.reduced <- lmer(formula(lmer.best), data=dataset.reduced, REML=FALSE)
car::Anova(lmer.best.reduced)

dataset.reduced.2 <- data.frame(dataset)
dataset.reduced.2 <- dataset.reduced %>% filter((Subject != "03") & (Subject != "14")  & (Subject != "11"))

cat('\n Best with 03, 14 and 11 removed (volume conduction)\n')
lmer.best.reduced.2 <- lmer(formula(lmer.best), data=dataset.reduced.2, REML=FALSE)
car::Anova(lmer.best.reduced.2)

sink()
close(anova_file)




###
# Limit waittime: 
# maximum allowed ISI = 4 s
max.allowed = 4
max.allowed = (max.allowed - min(dataset.raw$Waittime))^(1/4)
dataset.shortISI <- dataset[!is.nan(dataset$ISI) & dataset$ISI <= max.allowed,]

lmer.best.shortISI   <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL + ISI
                    + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI 
                    + Condition:MuPowerR:MuPowerL
                    + (1+Condition|Subject), #+ (1+Condition|Subject) + (1|Session), # The results are a bit stronger with this, but very likely overfitting
                    data=dataset.shortISI, REML=FALSE)
summary(lmer.best.shortISI)
car::Anova(lmer.best.shortISI)

sjRes.best.shortISI <- sjPlot::plot_model(lmer.best.shortISI, show.p = TRUE,show.value=TRUE, 
                                 vline.color = "gray90",value.offset=0.3,
                                 p.val="wald",
                                 dot.size=1, colors="black")  + 
  labs(title=element_blank(), y="Estimated fixed effect on ResponseFDI") +
  scale_x_discrete(labels=rev(rownames(as.data.frame(car::Anova(lmer.best.shortISI))))) +
  theme_bw() + 
  theme(axis.text=element_text(size=10), aspect.ratio = 1)
# correct p-values to Wald chi2 type II
p <- as.data.frame(car::Anova(lmer.best.shortISI))$`Pr(>Chisq)`
star.labels <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fe <- round(fixef(lmer.best.shortISI), digits=3)
sjRes.best.shortISI$data$p.stars <- c(star.labels)
sjRes.best.shortISI$data$p.label <- sprintf("%g %s", fe[2:length(fe)], c(star.labels))
sjRes.best.shortISI$data$p.value <- p

ipl.best.shortISI <- interact_plot(lmer.best.shortISI, pred = ISI, modx = Condition, 
                          interval = TRUE, int.width = 0.8, 
                          plot.points=TRUE, colors = "Set1") + 
  labs(x = expression(sqrt(sqrt(ISI - min))), y="ResponseFDI (z-score)") +
  theme_bw() + theme(aspect.ratio = 1)
ipl.best.shortISI$layers[[3]]$aes_params$alpha = 0.02


if(WRITE.PLOTS) {
  ggsave(
    plot = sjRes.best.shortISI+ipl.best.shortISI,
    filename = "bestModelOverview-cropped-ISI.pdf", # -cropped-ISI
    bg = "transparent", width=3000, height=1100,units="px"
  )
}
