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
library(corrplot)



#anova_file <- file("p-values.log", "w")

### 1: LOAD AND FORMAT DATA
# New phases and power:
dataset.raw <- read.delim("full_results__2023-06-29T173156.txt", sep=",")

dataset <- data.frame(dataset.raw)



# Channels (H for Hjorth):
# C3H   C4H
# FC3H  FC4H
# F3H   F4H
# O1H   O2H
dataset$PowerC3H <- log(dataset.raw$PowerC3H)
dataset$PowerC4H <- log(dataset.raw$PowerC4H)
dataset$PowerFC3H <- log(dataset.raw$PowerFC3H)
dataset$PowerFC4H <- log(dataset.raw$PowerFC4H)
dataset$PowerF3H <- log(dataset.raw$PowerF3H)
dataset$PowerF4H <- log(dataset.raw$PowerF4H)
dataset$PowerO1H <- log(dataset.raw$PowerO1H)
dataset$PowerO2H <- log(dataset.raw$PowerO2H)


dataset$ISI <- sqrt(sqrt(dataset$ISI - min(dataset.raw$ISI)))
#dataset <- dataset[dataset.raw$ISI > 2.35,]
dataset <- dataset[!is.infinite(dataset$ISI),]

dataset <- dataset[dataset$PreRangeAPB < 50 & dataset$PreRangeFDI < 50,]



dataset$SessionProgress <- NaN

for (ID in unique(dataset$Subject)) {
  for (ses in 1:2) {
    mask <- dataset$Subject == ID & dataset$Session == ses
    x <- dataset$ResponseFDI[mask]
    # 1) Nonlinear Transform
    x <- sqrt(sqrt(x))
    
    # 2) Antimedian filter
    window.radius <- 150
    padded <- c(rev(x[1:window.radius]), x, x[length(x)-(1:window.radius)+1])
    x <- x - zoo::rollmedian(padded, k=2*window.radius + 1)
    
    # 3) z-score
    x <- (x - mean(x)) / sd(x)
    
    dataset$ResponseFDI[mask] <- x
    dataset$SessionProgress[mask] <- (1:length(x))/length(x)
    
    # Try out within-session normalization of Power:
    pL <- dataset$PowerC3H[mask]
    pR <- dataset$PowerC4H[mask]
    dataset$PowerC3H[mask] <- (pL - mean(pL)) / sd(pL)
    dataset$PowerC4H[mask] <- (pR - mean(pR)) / sd(pR)
    
    pLPM1 <- dataset$PowerFC3H[mask]
    dataset$PowerFC3H[mask] <- (pLPM1 - mean(pLPM1)) / sd(pLPM1)
    pLPM2 <- dataset$PowerF3H[mask]
    dataset$PowerF3H[mask] <- (pLPM2 - mean(pLPM2)) / sd(pLPM2)
    pRPM1 <- dataset$PowerFC4H[mask]
    dataset$PowerFC4H[mask] <- (pRPM1 - mean(pRPM1)) / sd(pRPM1)
    pRPM2 <- dataset$PowerF4H[mask]
    dataset$PowerF4H[mask] <- (pRPM2 - mean(pRPM2)) / sd(pRPM2)
  
    pLO <- dataset$PowerO1H[mask]
    dataset$PowerO1H[mask] <- (pLO - mean(pLO)) / sd(pLO)
    pRO <- dataset$PowerO2H[mask]
    dataset$PowerO2H[mask] <- (pRO - mean(pRO)) / sd(pRO)
  }
  # Make Subject explicitly categorical
  dataset$Subject[dataset$Subject == ID] <- sprintf("%02d", ID) #format(ID)
}


# Check percentage of trials removed due to pre-innervation:
trials.kept <- table(dataset$Subject) / table(dataset.raw$Subject)


# Make Session clearly categorical
dataset["Session"][dataset["Session"] == 1] <- "1"
dataset["Session"][dataset["Session"] == 2] <- "2"

# Remove trials with outlier power (those are 5 trials in total):
sum(!(dataset$PowerC3H < 5 & dataset$PowerC3H > -5 & dataset$PowerC4H < 5 & dataset$PowerC4H > -5))
dataset <- dataset[dataset$PowerC3H < 5 & dataset$PowerC3H > -5 & dataset$PowerC4H < 5 & dataset$PowerC4H > -5,]



### 2: Inspect data:

ggplot(data=dataset, aes(x=ISI, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=ResponseFDI, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerC3H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerC4H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerFC3H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerFC4H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerF3H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerF4H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerO1H, color=Subject)) + geom_density()
ggplot(data=dataset, aes(x=PowerO2H, color=Subject)) + geom_density()

CV <- cov(dataset[c("PowerC3H", "PowerFC3H", "PowerF3H", "PowerO1H", "PowerC4H", "PowerFC4H", "PowerF4H", "PowerO2H")])
corrplot(CV, method="color", type="upper", col=colorRampPalette(c("#FF00FF", "#2000B0", "#000000", "#006655", "#EFFF00"))(200))
# We can see that:
# 1) No pair of powers is strongly correlated
# 2) Physically adjacent Hjorth montages yield correlated power (C3H-FC3H, FC3H-F3H, etc)
# 3) Power is correlated across hemispheres on a per-area basis (C3H-C4H, FC3H-FC4H, F3H-F4H, O1H-O2H)


ggplot(data=dataset, aes(x=PLVC4H, color=Condition))  + geom_density()
ggplot(data=dataset, aes(x=PLVFC3H, color=Condition)) + geom_density()
ggplot(data=dataset, aes(x=PLVFC4H, color=Condition)) + geom_density()
ggplot(data=dataset, aes(x=PLVF3H, color=Condition))  + geom_density()
ggplot(data=dataset, aes(x=PLVF4H, color=Condition))  + geom_density()
ggplot(data=dataset, aes(x=PLVO1H, color=Condition))  + geom_density()
ggplot(data=dataset, aes(x=PLVO2H, color=Condition))  + geom_density()




##### Power-overview

p.L <- ggplot(data=dataset, aes(x=PowerC3H, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("log \u03bc power (left)") +
  theme_bw()


p.R <- ggplot(data=dataset, aes(x=PowerC4H, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("log \u03bc power (right)") +
  theme_bw()


plot.r.pow.resp <- ggplot(dataset, aes(x=PowerC4H, y=ResponseFDI) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", option="A") +
  xlab("log \u03bc power (right)") +
  ylab("ResponseFDI (z-score)") +
  theme_bw()


plot.l.pow.resp <- ggplot(dataset, aes(x=PowerC3H, y=ResponseFDI) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis", option="A") +
  xlab("log \u03bc power (left)") +
  ylab("ResponseFDI (z-score)") +
  theme_bw()

p.R <- p.R + theme(legend.position="none")

ggsave(
  plot = plot.l.pow.resp + plot.r.pow.resp + p.L + p.R,
  filename = "power-data.pdf",
  bg = "transparent", width=2100, height=1500,units="px"
)



cor.test(dataset$PowerC3H, dataset$PowerC4H, method="pearson")
cor.test(dataset$ResponseFDI, dataset$PowerC3H, method="pearson")
cor.test(dataset$ResponseFDI, dataset$PowerC4H, method="pearson")

cor(dataset[c("ResponseFDI", "PowerC3H", "PowerC4H")], method="pearson")



### 3: Summary Models
cat('Naive model:\n')
naive.lmer <- lmer(ResponseFDI ~ 1+ Condition + (1+Condition|Subject), 
                   data = dataset, REML=FALSE)

car::Anova(naive.lmer)
naive.effect.size <- r.squaredGLMM(naive.lmer)[1] # R2m 

gdata <- dataset %>% group_by(Subject, Condition) %>% summarise(SD = sd(ResponseFDI), ResponseFDI = mean(ResponseFDI), n=n())
df.fx <- as.data.frame(effect('Condition', naive.lmer))
df.fx$Group = "all"
(gt <- ggplot(data=df.fx, aes(x=Condition, y=fit, group=Group)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .4, colour = NA, fill="dodgerblue") + 
    geom_line(lwd=1.5, color="dodgerblue4") +
    geom_point(size=2, color="dodgerblue4") + 
    geom_line(data=gdata, aes(x = Condition, y = ResponseFDI, group=Subject), alpha=0.4) + 
    geom_point(data=gdata, aes(x = Condition, y = ResponseFDI, group=Subject), alpha=0.4, size=1) + 
    theme_bw() + labs(y='ResponseFDI (z-Score)'))
ggsave(
  plot = gt,
  filename = "mainResult.pdf",
  bg = "transparent", width=1000, height=900,units="px"
)

cat('\n______________________________\nMixed (Power) model:\n')
mixed.lmer <- lmer(ResponseFDI ~ 1+ Condition*PowerC3H*PowerC4H + (1+Condition|Subject), # 
                   data = dataset, REML=FALSE)
car::Anova(mixed.lmer)
mixed.effect.size <- r.squaredGLMM(mixed.lmer)[1] # R2m 

cat('\n______________________________\nMaximal (Power,ISI) model:\n')
lmer.time <- lmer(ResponseFDI ~ 1 + Condition*PowerC3H*PowerC4H*ISI + (1+Condition|Subject), 
                  data = dataset, REML=FALSE)
car::Anova(lmer.time)

cat('\n______________________________\nSummary model:\n')
lmer.best   <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + PowerC4H + ISI
                    + Condition:PowerC3H + PowerC3H:PowerC4H + Condition:ISI 
                    + Condition:PowerC3H:PowerC4H
                    + (1+Condition|Subject), #+ (1+Condition|Subject) + (1|Session), # The results are a bit stronger with this, but very likely overfitting
                    data=dataset, REML=FALSE)
car::Anova(lmer.best)
summary(lmer.best)

cat('\n______________________________\nAkaike Information Criteria:\n')
AIC(naive.lmer, mixed.lmer, lmer.time, lmer.best)



##### 
# Individual results:
#####
condition.p.values <- c();
t.p.values <- c();
t.p.1 <- c();
t.p.2 <- c();

for(i in unique(dataset$Subject)) {
  dataset.ind <- data.frame(dataset)
  dataset.ind <- dataset.ind[dataset.ind$Subject == i,]
  lmer.s1 <- lm(ResponseFDI ~ 1 + Condition, data=dataset.ind)
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

individual.results <- data.frame(S1=t.p.1, S2=t.p.2, both=t.p.values)
cat(format_as_table(individual.results[order(t.p.values, decreasing = T),]))







# Model chain from basic model (Response ~ Condition) to summary model:
#
m0 <- lmer(ResponseFDI ~ 1 + (1|Subject), data=dataset, REML=F)
m0b <- lmer(ResponseFDI ~ 1 + Condition + (1|Subject), data=dataset, REML=F)
m1 <- lmer(ResponseFDI ~ 1 + Condition + (1+Condition|Subject), data=dataset, REML=F)

lrtest(m0,m0b)
lrtest(m0b,m1)
lrtest(m0,m1)
# All of these justify the inclusion of Condition and (Condition|Subject)


# Main effects: PowerC3H, PowerC4H, but not ISI!
m2 <- lmer(ResponseFDI ~ 1 + Condition + ISI + (1+Condition|Subject), data=dataset, REML=F)
m3 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + (1+Condition|Subject), data=dataset, REML=F)
m4 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + PowerC4H + (1+Condition|Subject), data=dataset, REML=F) 
lrtest(m1,m2)
lrtest(m1,m3)
lrtest(m3,m4)

m5 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + PowerC4H + Condition:PowerC3H + (1+Condition|Subject), data=dataset, REML=F) 
m6 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + PowerC4H + Condition:PowerC3H + PowerC3H:PowerC4H + (1+Condition|Subject), data=dataset, REML=F) 
m7 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H + PowerC4H + Condition:PowerC3H + Condition:PowerC4H + PowerC3H:PowerC4H + (1+Condition|Subject), data=dataset, REML=F) 

lrtest(m4,m5)
lrtest(m5,m6)
lrtest(m6,m7)

m8 <- lmer(ResponseFDI ~ 1 + Condition + ISI + PowerC3H + PowerC4H + Condition:ISI + Condition:PowerC3H + Condition:PowerC4H + PowerC3H:PowerC4H + (1+Condition|Subject), data=dataset, REML=F) 
lrtest(m7,m8)

m9 <- lmer(ResponseFDI ~ 1 + Condition + ISI + PowerC3H + PowerC4H + Condition:ISI + Condition:PowerC3H + Condition:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H + (1+Condition|Subject), 
           data=dataset, REML=F) 
lrtest(m8,m9)

m10 <- lmer(ResponseFDI ~ 1 + Condition + ISI + PowerC3H + PowerC4H 
           + Condition:ISI + Condition:PowerC3H + Condition:PowerC4H 
           + ISI:PowerC3H + ISI:PowerC4H
           + PowerC3H:PowerC4H + (1+Condition|Subject), 
           data=dataset, REML=F) 
lrtest(m9,m10)
m10 <- lmer(ResponseFDI ~ 1 + Condition + ISI + PowerC3H + PowerC4H 
            + Condition:ISI + Condition:PowerC3H + Condition:PowerC4H 
            + ISI:PowerC3H + PowerC3H:PowerC4H + Condition:PowerC3H:PowerC4H + (1+Condition|Subject), 
            data=dataset, REML=F) 
lrtest(m9,m10)


AIC(m0, m0b, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
AIC(naive.lmer, mixed.lmer, lmer.best, lmer.time, m9)
# m9 is the best in AIC, and justified by LR and Wald-II-tests.


sjRes.m9 <- sjPlot::plot_model(m9, show.p = TRUE,show.value=TRUE, 
                                 vline.color = "gray90",value.offset=0.3,
                                 dot.size=1, colors="black")  + 
  labs(title=element_blank(), y="Estimated fixed effect on ResponseFDI") +
  scale_x_discrete(labels=rev(rownames(as.data.frame(car::Anova(m9))))) +
  ylim(-0.35, 0.35) + 
  theme_bw() + 
  theme(axis.text=element_text(size=10), aspect.ratio = 1)
# correct p-values to Wald chi2 type II
p <- as.data.frame(car::Anova(m9))$`Pr(>Chisq)`
star.labels <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fe <- round(fixef(m9), digits=3)
sjRes.m9$data$p.stars <- star.labels
sjRes.m9$data$p.value <- p
sjRes.m9$data$p.label <- sprintf("%g %s", fe[2:length(fe)], star.labels)
sjRes.m9

ipl.ISI <- interact_plot(m9, pred = ISI, modx = Condition, 
                          interval = TRUE, int.width = 0.8, 
                          plot.points=TRUE, colors = "Set1", point.alpha=0.02) + 
  labs(x = expression(sqrt(sqrt(ISI - min))), y="ResponseFDI (z-score)") +
  theme_bw() + theme(aspect.ratio = 1)
ipl.ISI$layers[[3]]$aes_params$alpha = 0.02

ipl.PowerC3H <- interact_plot(m9, pred = PowerC3H, modx = Condition, 
                         interval = TRUE, int.width = 0.8, 
                         plot.points=TRUE, colors = "Set1", point.alpha=0.02) + 
  labs(x = "Power left M1 (log, z-score)", y="ResponseFDI (z-score)") +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
ipl.PowerC3H$layers[[3]]$aes_params$alpha = 0.02

ipl.PowerC4H <- interact_plot(m9, pred = PowerC4H, modx = Condition, 
                              interval = TRUE, int.width = 0.8, 
                              plot.points=TRUE, colors = "Set1", point.alpha=0.02) + 
  labs(x = "Power right M1 (log, z-score)", y="ResponseFDI (z-score)") +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "none")
ipl.PowerC4H$layers[[3]]$aes_params$alpha = 0.02


ggsave(
  plot = (sjRes.m9 | ipl.ISI) / (ipl.PowerC3H | ipl.PowerC4H),
  filename = "nonnaiveResults.pdf",
  bg = "transparent", width=3000, height=2200,units="px"
)


### Inspect collinearity:
check_model(m9)

int.PLV.ISI <- dataset$PLVC4H * dataset$ISI
int.PLV.PowerC3H <- dataset$PLVC4H * dataset$PowerC3H
int.PLV.PowerC4H <- dataset$PLVC4H * dataset$PowerC4H
int.ISI.PowerC3H <- dataset$ISI * dataset$PowerC3H
int.PowerC3H.PowerC4H <- dataset$PowerC3H * dataset$PowerC4H
predictor.names <- c("PLVC4H", "ISI", "PowerC3H", "PowerC4H", "PLVC4H:ISI", "PLVC4H:PowerC3H", "PLVC4H:PowerC4H", "ISI:PowerC3H", "PowerC3H:PowerC4H")

test.M <- matrix(c(dataset$PLVC4H, dataset$ISI, dataset$PowerC3H, dataset$PowerC4H, int.PLV.ISI, int.PLV.PowerC3H, int.PLV.PowerC4H, int.ISI.PowerC3H, int.PowerC3H.PowerC4H), ncol=length(predictor.names))
colnames(test.M) <- predictor.names
all(test.M[,1] == dataset$PowerC3H)
all(test.M[,2] == int.ISI.PowerC3H)
all(test.M[,3] == int.PLV.ISI)

Cor.predictors <- cor(test.M)
corrplot(Cor.predictors, method="color", col=colorRampPalette(c("#FF00FF", "#2000B0", "#000000", "#006655", "#EFFF00"))(200))






###                   ###
#                       #
#   L O C A T I O N S   #
#                       #
###                   ###

### Best model (m9) with posthoc PLV lM1--rM1:
m.big.C4H <- lmer(ResponseFDI ~ 1 + PLVC4H*ISI*PowerC3H*PowerC4H + (1+PLVC4H|Subject), data=dataset, REML=FALSE)
m.plv.C4H <- lmer(ResponseFDI ~ 1 + PLVC4H + ISI + PowerC3H + PowerC4H
                                 + PLVC4H:ISI + PLVC4H:PowerC3H + PLVC4H:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H
                                 + (1+PLVC4H|Subject),
                                 data=dataset, REML=FALSE)
car::Anova(m.plv.C4H) # Already differs from real-time Conditions: FC no longer interacts with Power on either side.
lrtest(m.big.C4H, m.plv.C4H) # To check if the big model improves on the best-type

m0.plv.C4H <- lmer(ResponseFDI ~ 1 + PLVC4H + ISI + PowerC3H + PowerC4H
                                   + PLVC4H:PowerC3H + PLVC4H:PowerC4H + ISI:PowerC3H + PowerC3H:PowerC4H
                                   + (1+PLVC4H|Subject), data=dataset, REML=FALSE) # W/o PLV:ISI interactions
lrtest(m.plv.C4H, m0.plv.C4H)
m0.plv.C4H <- lmer(ResponseFDI ~ 1 + PLVC4H + ISI + PowerC3H + PowerC4H
                   + PLVC4H:ISI + PLVC4H:PowerC3H + PLVC4H:PowerC4H + PowerC3H:PowerC4H
                   + (1+PLVC4H|Subject),
                   data=dataset, REML=FALSE) # w/O ISI:PowerC3H interaction
lrtest(m.plv.C4H, m0.plv.C4H)

# Note: From AIC(m.plv.C4H, m9) and lrtest(m.plv.C4H, m9), we can observe that m.plv.C4H is a worse model than
#       m9 [higher AIC, and lower logLikelihood]
AIC(m9, m.plv.C4H)
lrtest(m.plv.C4H, m9)

#####
##### P R E M O T O R:
#####
# Ipsilateral to stimulation side
m.big.FC3H <- lmer(ResponseFDI ~ 1 + PLVFC3H*ISI*PowerC3H*PowerFC3H + (1+PLVFC3H|Subject), data=dataset, REML=FALSE)
m.plv.FC3H <- lmer(ResponseFDI ~ 1 + PLVFC3H + ISI + PowerC3H + PowerFC3H
                  + PLVFC3H:ISI + PLVFC3H:PowerC3H + PLVFC3H:PowerFC3H + ISI:PowerC3H + PowerC3H:PowerFC3H
                  + (1+PLVFC3H|Subject),
                  data=dataset, REML=FALSE)
car::Anova(m.plv.FC3H) # Powers still interact
lrtest(m.big.FC3H, m.plv.FC3H) # To check if the big model improves on the best-type

# Contralateral to stimulation side
m.big.FC4H <- lmer(ResponseFDI ~ 1 + PLVFC4H*ISI*PowerC3H*PowerFC4H + (1+PLVFC4H|Subject), data=dataset, REML=FALSE)
m.plv.FC4H <- lmer(ResponseFDI ~ 1 + PLVFC4H + ISI + PowerC3H + PowerFC4H
                   + PLVFC4H:ISI + PLVFC4H:PowerC3H + PLVFC4H:PowerFC4H + ISI:PowerC3H + PowerC3H:PowerFC4H
                   + (1+PLVFC4H|Subject),
                   data=dataset, REML=FALSE)
car::Anova(m.plv.FC4H) # Powers still interact
lrtest(m.big.FC4H, m.plv.FC4H) # To check if the big model improves on the best-type


# More anterior:
# Ipsilateral to stimulation side
m.big.F3H <- lmer(ResponseFDI ~ 1 + PLVF3H*ISI*PowerC3H*PowerF3H + (1+PLVF3H|Subject), data=dataset, REML=FALSE)
m.plv.F3H <- lmer(ResponseFDI ~ 1 + PLVF3H + ISI + PowerC3H + PowerF3H
                   + PLVF3H:ISI + PLVF3H:PowerC3H + PLVF3H:PowerF3H + ISI:PowerC3H + PowerC3H:PowerF3H
                   + (1+PLVF3H|Subject),
                   data=dataset, REML=FALSE)
car::Anova(m.plv.F3H) # Interestingly: Main effect of PLVF3H; and powers still interact
lrtest(m.big.F3H, m.plv.F3H) # To check if the big model improves on the best-type

# Contralateral to stimulation side
m.big.F4H <- lmer(ResponseFDI ~ 1 + PLVF4H*ISI*PowerC3H*PowerF4H + (1+PLVF4H|Subject), data=dataset, REML=FALSE)
m.plv.F4H <- lmer(ResponseFDI ~ 1 + PLVF4H + ISI + PowerC3H + PowerF4H
                   + PLVF4H:ISI + PLVF4H:PowerC3H + PLVF4H:PowerF4H + ISI:PowerC3H + PowerC3H:PowerF4H
                   + (1+PLVF4H|Subject),
                   data=dataset, REML=FALSE)
car::Anova(m.plv.F4H) # Powers still interact
lrtest(m.big.F4H, m.plv.F4H) # To check if the big model improves on the best-type


#####
##### O C C I P I T A L
#####
# Ipsilateral to stimulation side
m.big.O1H <- lmer(ResponseFDI ~ 1 + PLVO1H*ISI*PowerC3H*PowerO1H + (1+PLVO1H|Subject), data=dataset, REML=FALSE)
m.plv.O1H <- lmer(ResponseFDI ~ 1 + PLVO1H + ISI + PowerC3H + PowerO1H
                  + PLVO1H:ISI + PLVO1H:PowerC3H + PLVO1H:PowerO1H + ISI:PowerC3H + PowerC3H:PowerO1H
                  + (1+PLVO1H|Subject),
                  data=dataset, REML=FALSE)
car::Anova(m.plv.O1H) # Interestingly: Main effect of PLVO1H; but powers don't interact
lrtest(m.big.O1H, m.plv.O1H) # To check if the big model improves on the best-type

# Contralateral to stimulation side
m.big.O2H <- lmer(ResponseFDI ~ 1 + PLVO2H*ISI*PowerC3H*PowerO2H + (1+PLVO2H|Subject), data=dataset, REML=FALSE)
m.plv.O2H <- lmer(ResponseFDI ~ 1 + PLVO2H + ISI + PowerC3H + PowerO2H
                  + PLVO2H:ISI + PLVO2H:PowerC3H + PLVO2H:PowerO2H + ISI:PowerC3H + PowerC3H:PowerO2H
                  + (1+PLVO2H|Subject),
                  data=dataset, REML=FALSE)
car::Anova(m.plv.O2H) # Nothing. Not even powers
lrtest(m.big.O2H, m.plv.O2H) # To check if the big model improves on the best-type


# In summary: Mostly nothing, powers do frequently interact [unlike in REFTEP]
# - Ipislateral Premotor (F3H) shows main effect of PLV (and no interaction of PLV with ISI!)
# - Ipsilateral Occipital
# - Powers do NOT interact for occipital cortex (O1, O2)


# Oh, notice: Does this also serve the ISI-investigation?
# Taking different regions of interest, we see no role of ISI (actually we see a trend in F3H)
# Is this helpful?
# + in REFTEP, we also see no interaction of PLV and ISI (but we do see a main effect of ISI in REFTEP)



















### PHASE:

p.C3ph <- ggplot(data=dataset, aes(x=PhaseC3H, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("phase (left)") +
  theme_bw()

p.C4ph <- ggplot(data=dataset, aes(x=PhaseC4H, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("phase (left)") +
  theme_bw()

p.C3ph / p.C4ph

dataset$CosPhaseC3 = cos(dataset$PhaseC3H)
dataset$SinPhaseC3 = sin(dataset$PhaseC3H)
dataset$CosPhaseC4 = cos(dataset$PhaseC4H)
dataset$SinPhaseC4 = sin(dataset$PhaseC4H)
m.ph.1 <- lmer(ResponseFDI ~ 1 + (1|Subject), data=dataset, REML=FALSE)
m.ph.2 <- lmer(ResponseFDI ~ 1 + PowerC3H + (1|Subject), data=dataset, REML=FALSE)
m.ph.3 <- lmer(ResponseFDI ~ 1 + PowerC3H + CosPhaseC3 + SinPhaseC3 + (1|Subject), data=dataset, REML=FALSE)
m.ph.4 <- lmer(ResponseFDI ~ 1 + PowerC3H*(CosPhaseC3 + SinPhaseC3) + (1|Subject), data=dataset, REML=FALSE)
m.ph.5 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H*(CosPhaseC3 + SinPhaseC3) + (1+Condition|Subject), data=dataset, REML=FALSE)
m.ph.6 <- lmer(ResponseFDI ~ 1 + Condition + PowerC3H*(CosPhaseC3 + SinPhaseC3) + Condition:(CosPhaseC3 + SinPhaseC3) + (1+Condition|Subject), data=dataset, REML=FALSE)

lrtest(m.ph.1, m.ph.2)
lrtest(m.ph.2, m.ph.3)
lrtest(m.ph.3, m.ph.4)
lrtest(m.ph.4, m.ph.5) # Neither Adding condition
lrtest(m.ph.5, m.ph.6) # Nor adding the interaction of phase and condition is a signficant improvement.


# PhaseC3-PowerL interaction
MuPowerL <- seq(-3,3, by=1)
phase <- seq(-pi,pi, by=pi/36)
ggAll <- expand.grid(PowerC3H=MuPowerL, Phase=phase)
ggAll$CosPhaseC3 <- cos(ggAll$Phase)
ggAll$SinPhaseC3 <- sin(ggAll$Phase)
ggAll$Response  <-predict(m.ph.4,newdata=ggAll,re.form=NA)


pAll <- ggplot(ggAll, aes(x=Phase, y=Response, color=PowerC3H, group=PowerC3H)) + 
  geom_line() +
  xlim(-pi,pi) + ylim(-1, 1) +
  scale_x_continuous(breaks = seq(-0.75*pi,pi,pi/4), 
                     labels=c("early rising", "rising", "late rising", "peak", "early falling", "falling", "late falling", "trough"),
                     guide=guide_axis(angle = 90)) +
  theme_minimal()
pAll

pAll.polar <- pAll + 
  coord_polar(theta = "x", start = pi, direction=-1)
pAll.polar


ggsave(
  plot = pAll,
  filename = "phase-power-interaction.pdf",
  bg = "transparent", width=1000, height=700,units="px"
)

ggsave(
  plot = pAll.polar,
  filename = "phase-power-interaction-polar.pdf",
  bg = "transparent", width=2100, height=1000,units="px"
)

ggsave(
  plot = pAll.polar + coord_polar(theta = "x", start = pi, direction=-1, clip="off") + theme(axis.text.y = element_blank(),
                                                                                             axis.ticks.y = element_blank()),
  filename = "phase-power-interaction-polar-poster.pdf",
  bg = "transparent", width=1100, height=800,units="px"
)





lmer.sessionProgress <- lmer(ResponseFDI ~ 1 + Condition*PowerC3H*PowerC4H*ISI*SessionProgress + (1+Condition|Subject), 
                  data = dataset, REML=FALSE)
car::Anova(lmer.sessionProgress)