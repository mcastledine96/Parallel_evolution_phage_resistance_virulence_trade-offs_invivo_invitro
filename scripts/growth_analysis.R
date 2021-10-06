### analysis of growth rates
library(lme4)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

#vitro
vitro <- read.csv("phage_therapy/growth_curve/data/vitro_final.csv", header = T)

#dataset guide:
  #treat = combined treatment and phage resistance level
  #well = well each clone occupied on the plate
  #m.fit = model fit from the rolling regression giving the estimate of growth rate
  #sd.fit = standard deviation of each growth rate
  #tr2 = which treatment each in vitro isolate belonged to. Control = control. mix_add = phage added repeatedly. mix = phage added once. 
  #phage resistance level = control = control (later changed to susceptible as all remained phage susceptible). resist = resistant to both phage. resist_sing = resistant to one phage. rep = which treatment replicate each isolate was taken from (MA = mix_add, M = mix, C = control)

#analysis
vitro2 <- vitro
vitro2$res[vitro2$res == "control"] <- "susc"

vit_mod3.1 <- lmer(m.fit ~ res * tr2 + (1|rep), data = vitro2)
vit_mod3.2 <- lmer(m.fit ~ res + tr2 + (1|rep), data = vitro2)
anova(vit_mod3.1, vit_mod3.2)
vit_mod3.3 <- lmer(m.fit ~ tr2 + (1|rep), data = vitro2)
anova(vit_mod3.3, vit_mod3.2)
vit_mod3.4 <- lmer(m.fit ~ res + (1|rep), data = vitro2)
anova(vit_mod3.2, vit_mod3.4) 
emmeans::emmeans(vit_mod3.3, pairwise ~ tr2)

#vivo
vivo <- read.csv("phage_therapy/growth_curve/data/vivo_011119.csv", header = T)

#bootstrap analysis

library(boot)

meanfun <- function(data, i){
  d <- data[i,]
  return(mean(d))   
}

#bootstrap - pre-phage treatment
andat <- subset(vivo, treat=="anc") #subset data based on ancestral
andat <- andat$m.fit #extract growth rates
andat <- data.frame(andat) #convert to dataframe
data <- andat #rename to data so function can work

bo <- boot(data[, "andat", drop = FALSE], statistic=meanfun, R=5000)
bo #mean #0.3556
boot.ci(bo, conf=0.95, type="norm") #( 0.3420,  0.3691 )

#post-phage treatment
vivo_t1 <- subset(vivo, treat == "TP1") #subset data based on t1
vivo_t2 <- subset(vivo, treat == "TP2")
vivo_t12 <- merge(vivo_t1, vivo_t2, all = "T")
vivo_t12 <- vivo_t12$m.fit #extract growth rates
vivo_t12 <- data.frame(vivo_t12) #convert to dataframe
data <- vivo_t12 

bo <- boot(data[, "vivo_t12", drop = FALSE], statistic=meanfun, R=5000)
bo #mean 0.249
boot.ci(bo, conf=0.95, type="norm") #( 0.1949,  0.3030 )

#analyse resistance levels within post-phage treatment groups
clone_res <- function(x) {
  case_when(x == " E03" ~ "susceptible",
            x == " E04" ~ "resistant",
            x == " E05" ~ "resistant",
            x == " E06" ~ "resistant",
            x == " E07" ~ "susceptible",
            x == " E08" ~ "resistant",
            x == " E01" ~ "resistant",
            x == " E02" ~ "resistant",
            x == " F01" ~ "susceptible",
            x == " F02" ~ "resistant")
}

vivo_noa <- filter(vivo, ! treat == "anc")
clones <- as.character(vivo_noa$Well)
clones <- clone_res(clones)
an <- character(24)
an[1:24] <- "Ancestral"
clone <- c(an, clones)
vivo$res <- clone

#bootstrap of resistant clones
resist_vivo <- subset(vivo, res == "resistant") 
resist_vivo <- resist_vivo$m.fit #extract growth rates
resist_vivo <- data.frame(resist_vivo) #convert to dataframe
data <- resist_vivo 

bo <- boot(data[, "resist_vivo", drop = FALSE], statistic=meanfun, R=5000)
bo #mean
boot.ci(bo, conf=0.95, type="norm")

#bootstrap of susceptible clones
susc_vivo <- subset(vivo, res == "susceptible") #subset data based on t1
susc_vivo <- susc_vivo$m.fit #extract growth rates
susc_vivo <- data.frame(susc_vivo) #convert to dataframe
data <- susc_vivo

bo <- boot(data[, "susc_vivo", drop = FALSE], statistic=meanfun, R=5000)
bo #mean
boot.ci(bo, conf=0.95, type="norm")

#plots
clone_res <- function(x) {
  case_when(x == " E03" ~ "susceptible",
            x == " E04" ~ "resistant",
            x == " E05" ~ "resistant",
            x == " E06" ~ "single_resistant",
            x == " E07" ~ "susceptible",
            x == " E08" ~ "single_resistant",
            x == " E01" ~ "resistant",
            x == " E02" ~ "resistant",
            x == " F01" ~ "susceptible",
            x == " F02" ~ "resistant")
}

vivo2 <- vivo
clones <- vivo_noa$Well
clones <- clone_res(clones)
an <- character(24)
an[1:24] <- "susceptible"
clones <- c(an, clones)
vivo2$res <- clones

vivo2$treat[vivo2$treat == "TP1"] <- "pt"
vivo2$treat[vivo2$treat == "TP2"] <- "pt"

g1 <- ggplot(vivo2, aes(x = treat, y = m.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.1), size = 3, aes(fill = res), shape = 21)  +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_text(size = 16), axis.text.x = element_blank(), legend.text = element_text(size = 13), legend.title = element_text(size = 14), legend.position = "none", plot.title = element_text(hjust = 0, size = 16), axis.title.x = element_blank()) +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 3, labels = c("Resistant", "Resistant\n(one phage)", "Susceptible")) +
  scale_y_continuous(limits = c(0.04, 0.52)) +
  xlab("Treatment") +
  scale_x_discrete(labels = c("Pre-phage", "Post-phage")) +
  ylab(expression(paste("Growth rate ", "(", hr^-1, ")"))) +
  labs(fill = "Resistance", title = expression(paste("(a)")))
  
g2 <- ggplot(vitro2, aes(x = tr2, y = m.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.1), size = 3, aes(fill = res), shape = 21)  +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.title = element_text(size = 16), axis.text.x = element_blank(), legend.text = element_text(size = 13), legend.title = element_text(size = 14), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0, size = 16), axis.title.x = element_blank()) +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 3) +
  scale_y_continuous(limits = c(0.04, 0.52)) +
  scale_x_discrete(labels = c("Control", "Phage added\nonce", "Phage added\nrepeatedly")) +
  labs(fill = "Resistance", title = expression(paste("(b)"))) +
  xlab("Treatment") 

growps <- g1 + g2 + plot_layout(ncol = 2)
growps

