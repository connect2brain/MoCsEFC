setwd(file.path(path.expand(dirname("~")), "Downloads"))
library(tidyverse)
library(ggplot2)
library(qqplotr)
library(patchwork)
library(abind)
library(RColorBrewer)

WRITE.PLOTS = T


dataset <- read.delim("realtime-cPLV.csv", sep=",")




p.posthoc <- ggplot(data=dataset, aes(x=PLV, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() +
  theme(legend.position="none")
p.posthoc$theme$panel.background$fill  = "transparent"
p.posthoc$theme$plot.background$fill   = "transparent"
p.posthoc$theme$legend.background$fill = "transparent"

p.posthoc.im <- ggplot(data=dataset, aes(x=iPLV, group=Condition, fill=Condition, color=Condition)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw()
p.posthoc.im$theme$panel.background$fill  = "transparent"
p.posthoc.im$theme$plot.background$fill   = "transparent"
p.posthoc.im$theme$legend.background$fill = "transparent"



d.ph <- Arg(dataset$RealPart + 1i*dataset$ImagPar)
df <- data.frame(Condition = c(dataset$Condition, dataset$Condition, dataset$Condition), 
                 PhaseShift = c(d.ph - 2*pi, d.ph, d.ph+2*pi))


pl.joint <- ggplot(data=dataset, mapping=aes(color=Condition,x=RealPart, y=ImagPart,fill=Condition)) + 
  geom_point(alpha=0.03,size=1) + 
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() + theme(legend.position = "none")
pl.angle <- ggplot(data=df, mapping=aes(color=Condition,x=PhaseShift,fill=Condition)) + 
  geom_density(adjust=0.5,alpha=0.4) + 
  coord_cartesian(xlim=c(-pi, pi), ylim=c(0, 0.1), expand=0) + 
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() + theme(legend.position = "none")
pl.imag <- ggplot(data=dataset, mapping=aes(color=Condition,x=ImagPart,fill=Condition)) + 
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  coord_flip() +
  theme_bw()

(p.posthoc | p.posthoc.im) / (pl.angle | pl.joint | pl.imag)





if(WRITE.PLOTS) {
  ggsave(
    plot = (p.posthoc | p.posthoc.im) / (pl.angle | pl.joint | pl.imag),
    filename = "Phaseshift-inspection-realtime.pdf",
    bg = "transparent", width=4000, height=2400,units="px"
  )
  
  ggsave(
    plot = p.posthoc,
    filename = "rtPLV-by-condition.pdf",
    bg = "transparent", width=1700, height=800,units="px"
  )
  
  ggsave(
    plot = p.posthoc.im,
    filename = "rtImPLV-by-condition.pdf",
    bg = "transparent", width=1700, height=800,units="px"
  )
}


densities <- c()
d.at.zero <- c()

for(id in unique(dataset$Subject)) {
  dataset.ind <- dataset[dataset$Subject == id,]
  
  p.posthoc <- ggplot(data=dataset.ind, aes(x=PLV, group=Condition, fill=Condition, color=Condition)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values=c("red", "blue")) +
    scale_color_manual(values=c("red", "blue")) +
    theme_bw() +
    theme(legend.position="none")
  p.posthoc$theme$panel.background$fill  = "transparent"
  p.posthoc$theme$plot.background$fill   = "transparent"
  p.posthoc$theme$legend.background$fill = "transparent"
  
  p.posthoc.im <- ggplot(data=dataset.ind, aes(x=iPLV, group=Condition, fill=Condition, color=Condition)) +
    geom_density(adjust=1.5, alpha=.4) +
    scale_fill_manual(values=c("red", "blue")) +
    scale_color_manual(values=c("red", "blue")) +
    theme_bw()
  p.posthoc.im$theme$panel.background$fill  = "transparent"
  p.posthoc.im$theme$plot.background$fill   = "transparent"
  p.posthoc.im$theme$legend.background$fill = "transparent"
  
  
  
  d.ph <- Arg(dataset.ind$RealPart + 1i*dataset.ind$ImagPar)
  df <- data.frame(Condition = c(dataset.ind$Condition, dataset.ind$Condition, dataset.ind$Condition), 
                   PhaseShift = c(d.ph - 2*pi, d.ph, d.ph+2*pi))
  
  
  pl.joint <- ggplot(data=dataset.ind, mapping=aes(color=Condition,x=RealPart, y=ImagPart,fill=Condition)) + 
    geom_point(alpha=0.1,size=1) + 
    scale_fill_manual(values=c("red", "blue")) +
    scale_color_manual(values=c("red", "blue")) +
    theme_bw() + theme(legend.position = "none")
  pl.angle <- ggplot(data=df, mapping=aes(color=Condition,x=PhaseShift,fill=Condition)) + 
    geom_density(adjust=0.5,alpha=0.4) + 
    coord_cartesian(xlim=c(-pi, pi), ylim=c(0, 0.1), expand=0) + 
    scale_fill_manual(values=c("red", "blue")) +
    scale_color_manual(values=c("red", "blue")) +
    theme_bw() + theme(legend.position = "none")
  inspect <- ggplot_build(pl.angle);
  
  pl.imag <- ggplot(data=dataset.ind, mapping=aes(color=Condition,x=ImagPart,fill=Condition)) + 
    geom_density(alpha=0.4) +
    scale_fill_manual(values=c("red", "blue")) +
    scale_color_manual(values=c("red", "blue")) +
    coord_flip() +
    theme_bw()

  x.angle <- inspect$data[[1]]$x[inspect$data[[1]]$x > -pi & inspect$data[[1]]$x < pi & inspect$data[[1]]$group == 1];
  angle.density <- inspect$data[[1]]$density[inspect$data[[1]]$x > -pi & inspect$data[[1]]$x < pi & inspect$data[[1]]$group == 1];
  angle.density <- angle.density / sum(angle.density);
  d.at.zero <- c(d.at.zero, angle.density[which.min(abs(x.angle))])
  cat(sprintf("%02d : %f \n", id, angle.density[which.min(abs(x.angle))]))
  densities <- abind(densities, angle.density, along=2)
  
  if(WRITE.PLOTS) {
    ggsave(
      plot = (p.posthoc | p.posthoc.im) / (pl.angle | pl.joint | pl.imag),
      filename = sprintf("ps-i-realtime-%02d.pdf", id),
      bg = "transparent", width=4000, height=2400,units="px"
    )
  }
}

colors <- rgb(seq(0,255,255/15)/255, 0, rev(seq(0,255,255/15)/255))

pdf("angle-inclination-individual.pdf", width=8, height=5, bg="transparent")
plot(x=NULL, y=NULL,xlim=c(-pi, pi), ylim=c(0,0.01), xlab="Argument of cPLV", ylab="density")
i <- 1
for(id in order(d.at.zero)) {
    lines(x.angle,densities[,id],col=colors[i],lwd=3)
    i <- i+1;
}
dev.off()

# Based on density at zero (d.at.zero):
# thus removing 3, 14 and potentially 11 (biased towards 180deg), and checking if the statistics hold up