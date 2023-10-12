setwd(file.path(path.expand(dirname("~")), "Downloads"))

dataset.raw <- read.delim("results___2023-06-13T171710.txt", sep=",")


#dataset.realtime.aug <- read.delim("realtime-cPLV.csv", sep=",")

dataset <- data.frame(dataset.raw)
dataset$MuPowerR <- log(dataset.raw$MuPowerR)
dataset$MuPowerL <- log(dataset.raw$MuPowerL)

dataset$Waittime <- dataset.raw$Waittime - min(dataset.raw$Waittime) # remove min ISI
dataset$Waittime <- sqrt(sqrt(dataset$Waittime))
dataset$Waittime[dataset.raw$Waittime == 10] <- NaN

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










dataset$Responder <- F

subjects <- unique(dataset$Subject)
for(i in 1:length(subjects)) {
  s <- subjects[i]
  dataset$Responder[dataset$Subject == s] <- ifelse(t.p.values[i] < 0.05, "responder", "nonresponder")
}

ggplot(data=dataset, aes(x=Responder, y=ResponseFDI, fill=Condition)) + geom_boxplot()

(l.power.by.responder <- ggplot(data=dataset, aes(x=MuPowerL, group=Responder, color=Responder, fill=Responder)) + 
  geom_density(alpha=.5) +
  scale_fill_manual(values=c("blue", "green")) +
  scale_color_manual(values=c("blue", "green")) + xlim(-5, 5) +
  theme_bw() + theme(legend.position = "none"))
wilcox.test(MuPowerL ~ Responder, dataset)
(r.power.by.responder <- ggplot(data=dataset, aes(x=MuPowerR, group=Responder, color=Responder, fill=Responder)) + 
  geom_density(alpha=.5) +
  scale_fill_manual(values=c("blue", "green")) +
  scale_color_manual(values=c("blue", "green")) +
  theme_bw())
wilcox.test(MuPowerR ~ Responder, dataset)

raw.response.by.responder <- ggplot(data=dataset, aes(x=ResponseFDI, group=Responder, color=Responder, fill=Responder)) + 
       geom_density(alpha=.5) +
       scale_fill_manual(values=c("blue", "green")) +
       scale_color_manual(values=c("blue", "green")) + xlim(1,11) +
       theme_bw()




ggsave(
  plot = l.power.by.responder | r.power.by.responder,
  filename = "Power-by-responder.pdf",
  bg = "transparent", width=2000, height=800,units="px"
)
ggsave(
  plot = l.power.by.responder | r.power.by.responder,
  filename = "Power-by-responder.png",
  bg = "transparent", width=2000, height=800,units="px"
)

ggsave(
  plot = raw.response.by.responder,
  filename = "RawResponse-by-responder.pdf",
  bg = "transparent", width=1200, height=800,units="px"
)





# Checking I/O-curve-results:
io.data <- read.delim('reviewer1_MEPIO.txt', sep=",")
io.data$Responder <- ifelse(io.data$Responder, "responder", "nonresponder")

#io.data$MEP <- io.data$MEP / max(io.data$MEP)
(io.curve.plot <- ggplot(data=io.data, aes(x=Intensity, y=MEP, color=Responder, fill=Responder)) + 
  geom_point(alpha=0.2,size=0.5) + 
  geom_smooth(method = glm, formula=y~x, se=T, method.args = list(family = "quasipoisson"), linewidth=0.5) + 
  scale_fill_manual(values=c("blue", "green")) +
  scale_color_manual(values=c("blue", "green")) + 
  theme_bw())

ggsave(
  plot = io.curve.plot,
  filename = "IOC-by-responder.pdf",
  bg = "transparent", width=1300, height=800,units="px"
)

(io.curve.individual <- ggplot(data=io.data, aes(x=Intensity, y=MEP, color=Subject, fill=Responder)) + 
       geom_smooth(method = glm, formula=y~x, se=T, method.args = list(family = "quasipoisson"), linewidth=0.1, alpha=0.2) + 
       scale_fill_manual(values=c("blue", "green")) +
       theme_bw() + theme(legend.position = "none"))
ggsave(
  plot = io.curve.individual,
  filename = "IOC_individual-by-responder.pdf",
  bg = "transparent", width=900, height=800,units="px"
)




# Stimulation intensities:
si.1 <- c(59,58,56,68,58,75,54,45,64,47,58,65,59,62,65)
si.2 <- c(53,62,54,56,59,76,53,51,66,46,63,67,62,65,74)
si.data <- data.frame(SI = c(si.1, si.2), Responder=ifelse(c(t.p.1, t.p.2) < 0.05, "responder", "nonresponder"))

(si.comparison <- ggplot(data=si.data, aes(x=Responder, y=SI,fill=Responder)) +
    geom_boxplot(alpha=0.5) + scale_fill_manual(values=c("blue", "green")) +
    theme_bw() + theme(axis.text.x=element_blank()))
ggsave(
  plot = si.comparison,
  filename = "SI-by-responder.pdf",
  bg = "transparent", width=1000, height=800,units="px"
)



io.curve.plot <- io.curve.plot + theme(legend.position="none")
raw.response.by.responder <- raw.response.by.responder + theme(legend.position="none")
(io.curve.individual / io.curve.plot) | (raw.response.by.responder / si.comparison)
ggsave(
  plot = (io.curve.individual / io.curve.plot) | (raw.response.by.responder / si.comparison),
  filename = "SI-responder-4x4.png",
  bg = "transparent", width=2000, height=1600,units="px"
)
