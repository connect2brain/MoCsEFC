#####
p.C3ph <- ggplot(data=dataset, aes(x=PosthocPhaseC3, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("phase (left)") +
  theme_bw()

p.C4ph <- ggplot(data=dataset, aes(x=PosthocPhaseC4, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  xlab("phase (left)") +
  theme_bw()

p.C3ph / p.C4ph



#####
# note: there was a 5ms offset; however, phastimate provides an offset-parameter 
# that allows us to correct for this offset.

dataset$CosPhaseC3 = cos(dataset$PosthocPhaseC3)
dataset$SinPhaseC3 = sin(dataset$PosthocPhaseC3)
dataset$CosPhaseC4 = cos(dataset$PosthocPhaseC4)
dataset$SinPhaseC4 = sin(dataset$PosthocPhaseC4)

#####
# See that C4 phase is unpredictive. Thus remove all C4-terms:
lmer.best.n.phase.bihemispheric <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL
                          + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL # + (CosPhaseC3 + SinPhaseC3):Condition # this does not add significant interactions
                          + (1+Condition|Subject),
                          data=dataset, REML=FALSE)
car::Anova(lmer.best.n.phase.bihemispheric)


anova_file <- file("p-values.log", "a")

## LR tests to check "adding phase" and "adding phase:power"
lmer.best   <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL
                    + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI 
                    + (1+Condition|Subject),
                    data=dataset, REML=FALSE)
car::Anova(lmer.best)
sink(anova_file,append=T)
cat("\n\n\n\n\n")
cat(rep('#', 50), sep="#", end="\n")
cat('\n P H A S E \n\n')
cat('\n______________________________\nSummary Model with only phase:\n')
lmer.best.w.phase <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL
                             + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI + CosPhaseC3 + SinPhaseC3
                             + (1+Condition|Subject),
                             data=dataset, REML=FALSE)
car::Anova(lmer.best.w.phase)
cat('\n______________________________\nSummary Model with phase, power interaction:\n')
lmer.best.w.phase.power <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL
                             + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL
                             + (1+Condition|Subject),
                             data=dataset, REML=FALSE)
car::Anova(lmer.best.w.phase.power)
cat('\n______________________________\nSummary Model with all interactions:\n')
lmer.best.w.all.interactions <- lmer(ResponseFDI ~ 1 + Condition + MuPowerR + MuPowerL
                             + Condition:MuPowerL + MuPowerR:MuPowerL + Condition:ISI + (CosPhaseC3 + SinPhaseC3)*MuPowerL + (CosPhaseC3 + SinPhaseC3):Condition
                             + (1+Condition|Subject),
                             data=dataset, REML=FALSE)
car::Anova(lmer.best.w.all.interactions)
cat(rep('_', 50), sep="_", end="\n")
anova(lmer.best.w.phase, lmer.best)
anova(lmer.best.w.phase.power, lmer.best.w.phase)
anova(lmer.best.w.all.interactions, lmer.best.w.phase.power)
sink()
close(anova_file)


# PhaseC3-PowerL interaction
MuPowerL <- seq(-3,3, by=1)
phase <- seq(-pi,pi, by=pi/36)

lmer.phase <- lmer(ResponseFDI ~ 1 + (CosPhaseC3 + SinPhaseC3)*MuPowerL + (1+Condition|Subject), data=dataset, REML=FALSE)
car::Anova(lmer.phase)

ggAll <- expand.grid(MuPowerL=MuPowerL, Phase=phase)
ggAll$CosPhaseC3 <- cos(ggAll$Phase)
ggAll$SinPhaseC3 <- sin(ggAll$Phase)
ggAll$Response  <-predict(lmer.phase,newdata=ggAll,re.form=NA)


pAll <- ggplot(ggAll, aes(x=Phase, y=Response, color=MuPowerL, group=MuPowerL)) + 
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

if(WRITE.PLOTS) {
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
  
}










# Individual phase-response

phase <- seq(-pi,pi, by=pi/36)
for(i in unique(dataset$Subject)) {
  dataset.ind <- data.frame(dataset)
  dataset.ind <- dataset.ind[dataset.ind$Subject == i,]
  lmer.phase.ind <- lm(ResponseFDI ~ 1 + (CosPhaseC3 + SinPhaseC3)*MuPowerL, data=dataset.ind, REML=FALSE)
  
  MuPowerL <- seq(-3,3, by=1)
  # OR:
  MuPowerL = seq(quantile(dataset.ind$MuPowerL, 0.10), quantile(dataset.ind$MuPowerL, 0.90), length.out=7)
  
  ggAll.ind <- expand.grid(MuPowerL=MuPowerL, Phase=phase)
  ggAll.ind$CosPhaseC3 <- cos(ggAll.ind$Phase)
  ggAll.ind$SinPhaseC3 <- sin(ggAll.ind$Phase)
  ggAll.ind$Response  <-predict(lmer.phase.ind,newdata=ggAll.ind,re.form=NA)
  
  
  pAll <- ggplot(ggAll.ind, aes(x=Phase, y=Response, color=MuPowerL, group=MuPowerL)) + 
    geom_line() +
    xlim(-pi,pi) + ylim(-1, 1) +
    scale_x_continuous(breaks = seq(-0.75*pi,pi,pi/4), 
                       labels=c("early rising", "rising", "late rising", "peak", "early falling", "falling", "late falling", "trough"),
                       guide=guide_axis(angle = 90)) +
    theme_minimal()
  invisible(cat(i, end="\n"))
  print(pAll)
  Sys.sleep(1)
  invisible(readline(prompt="Press [enter] to continue"))
}



# Quick note:
# The mu-powers differ a fair lot between participants.
ggplot(dataset, aes(x=MuPowerL, group=Subject, fill=Subject, color=Subject)) + geom_density(alpha=0.1) + theme_minimal()
# Interpretational danger: What holds on population level need not be true on the 
# individual level (individual may just be a sample from the population level results, 
# e.g. individuals with low mu power might show different effect than individuals 
# with high power)
# CHECK THIS.




















#####


# Conclusion: In none of the cases is the interaction of C3-Phase and FC significant!






##### 
# Plots:
# - Power-Condition interaction
# - FC-ISI interaction
# - Power-Power Interaction










# Power-Power plot:
MuPowerL <- seq(-5,5, by=0.5)
MuPowerR <- seq(-5,5, by=0.5)
ggLow <- expand.grid(MuPowerL=MuPowerL, MuPowerR=MuPowerR)
Condition <- rep("low",nrow(ggLow))
ISI <- rep(1.25, nrow(ggLow))
ggLow$Codition <- Condition
ggLow$ISI <- ISI

ggHigh <- expand.grid(MuPowerL=MuPowerL, MuPowerR=MuPowerR)
ggHigh$Condition <- rep("high", nrow(ggHigh))
ggHigh$ISI <- ISI

# prediction from the linear model
ggLow$Response <-predict(lmer.best,newdata=ggLow,re.form=NA)
ggHigh$Response <-predict(lmer.best,newdata=ggHigh,re.form=NA)



pLow <- ggplot(ggLow, aes(x=MuPowerR, y=MuPowerL, z=Response, fill=Response)) + 
  geom_raster(interpolate=TRUE) +
  geom_contour(colour = "white", binwidth = 0.1) +
  geom_contour(size=1, breaks=0, colour="white") +
  scale_fill_gradientn(colours=c("#1f114f", "#0000FF", "#388bff")) +
  theme_minimal() +
  ggtitle("Low")
pHigh <- ggplot(ggHigh, aes(x=MuPowerR, y=MuPowerL, z=Response, fill=Response)) + 
  geom_raster(interpolate=TRUE) +
  geom_contour(colour = "white", binwidth = 0.1) +
  geom_contour(size=1, breaks=0, colour="white") +
  scale_fill_gradientn(colours=c("#59142d", "#FF0000", "#fc822b")) +
  theme_minimal() +
  ggtitle("High")
pHigh / pLow

ggsave(
  plot = pHigh / pLow,
  filename = "PowerInteraction.png",
  bg = "transparent", width=1600, height=2800,units="px"
)







# Power-Power "animation" by phase
frame <- 0
for(ph in seq(-pi, pi, 0.1)) {


MuPowerL <- seq(-5,5, by=0.5)
MuPowerR <- seq(-5,5, by=0.5)
ggLow <- expand.grid(MuPowerL=MuPowerL, MuPowerR=MuPowerR)
Condition <- rep("low",nrow(ggLow))
ISI <- rep(1.25, nrow(ggLow))
phase <- rep(ph, nrow(ggLow)) #runif(nrow(ggLow), min=-pi, max=pi)
ggLow$CosPhaseC3 <- cos(phase)
ggLow$SinPhaseC3 <- sin(phase)
ggLow$Codition <- Condition
ggLow$ISI <- ISI

ggHigh <- expand.grid(MuPowerL=MuPowerL, MuPowerR=MuPowerR)
ggHigh$Condition <- rep("high", nrow(ggHigh))
ggHigh$ISI <- ISI
ggHigh$CosPhaseC3 <- cos(phase)
ggHigh$SinPhaseC3 <- sin(phase)

# prediction from the linear model
ggLow$Response <-predict(lmer.best.n.phase,newdata=ggLow,re.form=NA)
ggHigh$Response <-predict(lmer.best.n.phase,newdata=ggHigh,re.form=NA)



pLow <- ggplot(ggLow, aes(x=MuPowerR, y=MuPowerL, z=Response, fill=Response)) + 
  geom_raster(interpolate=TRUE) +
  geom_contour(colour = "white", binwidth = 0.1) +
  geom_contour(size=1, breaks=0, colour="white") +
  scale_fill_gradientn(colours=c("#1f114f", "#0000FF", "#388bff")) +
  theme_minimal() +
  ggtitle("Low")
pHigh <- ggplot(ggHigh, aes(x=MuPowerR, y=MuPowerL, z=Response, fill=Response)) + 
  geom_raster(interpolate=TRUE) +
  geom_contour(colour = "white", binwidth = 0.1) +
  geom_contour(size=1, breaks=0, colour="white") +
  scale_fill_gradientn(colours=c("#59142d", "#FF0000", "#fc822b")) +
  theme_minimal() +
  ggtitle("High")
pHigh / pLow

ggsave(
  plot = pHigh / pLow,
  filename = sprintf("powerInteractionAnimation/PowerInteraction_%0.3d.png", frame),
  bg = "transparent", width=1600, height=2800,units="px"
)

frame <- frame + 1
}
