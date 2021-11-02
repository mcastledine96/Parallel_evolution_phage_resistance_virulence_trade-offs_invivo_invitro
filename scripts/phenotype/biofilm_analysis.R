##biofilm analysis

library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(MicrobioUoE)
library(tidyverse)
library(janitor)
library(patchwork)
library(brms)
library(tidybayes)

biofilm <- read.csv("phage_therapy/biofilm_assay/data/biofilm_wrangled.csv", header = T, stringsAsFactors = F)

#dataset guide
  #flor = fluorescence measure
  #res = phage resistance. res = resistant to both phage. sing_res = resistant to one phage. susc = susceptible
  #lab = whether isolates belong to in vitro or in vivo study. Media for blank controls
  #treat = condition with the in vivo or in vitro study. In vitro: mix = phage added once, mix_add = phage added repeatedly, control (no phage). In vivo: T1 = isolate from day 2. T2 = isolate from day 4. 
  #rep = technical replicate (n = 3)
  #well = Well within 96-well plate
  #tr_rep = refers to in vitro treatment reps clones were isolated from. 

biofilm$sp_rep <- interaction(biofilm$res, biofilm$lab, biofilm$ind_rep)

bio_means <- group_by(biofilm, res, lab, treat, tr_rep, ind_rep, sp_rep) %>%
  summarise(m_flor = mean(flor))

bio_means_noc <- filter(bio_means, ! lab == "media")

#analysis of the in vitro data

#analyse raw data
vit_biof <- filter(bio_means, lab == "vitro")
vit_biof$res[vit_biof$res == "control"] <- "susc"

vit_biof$random <- interaction(vit_biof$treat, vit_biof$tr_rep) #include a random effect to account for the fact some clones came from the same treatment replicate

vit_biof$res <- as.factor(vit_biof$res) #change to factors for analysis

#removing both outliers
vit_biof2 <- filter(vit_biof, ! m_flor == 711) #removing susceptible outlier
vit_biof3 <- filter(vit_biof, ! m_flor == 708.6) #removing control outlier
vit_biof4 <- filter(vit_biof3, ! m_flor == 711) #both removed

#analysis with both outliers removed
bmod7 <- lmer(log(m_flor) ~ res * treat + (1|random), data = vit_biof4)
bmod8 <- lmer(log(m_flor) ~ res + treat + (1|random), data = vit_biof4)
anova(bmod7, bmod8) 
bmod9 <- lmer(log(m_flor) ~ treat + (1|random), data = vit_biof4)
anova(bmod8, bmod9)
bmod10 <- lmer(log(m_flor) ~ res + (1|random), data = vit_biof4)
anova(bmod8, bmod10) 
bmod11 <- lmer(log(m_flor) ~ 1 + (1|random), data = vit_biof4)
anova(bmod10, bmod11) 
emmeans(bmod10, pairwise ~ res) 

#analysis with no outliers removed
bmod1 <- lmer(log(m_flor) ~ res * treat + (1|random), data = vit_biof)
bmod1a <- lmer(log(m_flor) ~ res + treat + (1|random), data = vit_biof)
anova(bmod1, bmod1a)
bmod1b <- lmer(log(m_flor) ~ treat + (1|random), data = vit_biof)
anova(bmod1a, bmod1b) 
bmod1c <- lmer(log(m_flor) ~ res + (1|random), data = vit_biof)
anova(bmod1a, bmod1c) 
bmod2 <- lmer(log(m_flor) ~ 1 + (1|random), data = vit_biof)
anova(bmod1c, bmod2) 
emmeans::emmeans(bmod1c, pairwise ~ res)

#outlier in the susceptible resistance group removed
bmod3 <- lmer(log(m_flor) ~ res * treat + (1|random), data = vit_biof2)
bmod4 <- lmer(log(m_flor) ~ res + treat + (1|random), data = vit_biof2)
anova(bmod3, bmod4) 
bmod3a <- lmer(log(m_flor) ~ treat + (1|random), data = vit_biof2)
anova(bmod4, bmod3a) 
bmod3b <- lmer(log(m_flor) ~ res + (1|random), data = vit_biof2)
anova(bmod4, bmod3b) 
emmeans(bmod4, pairwise ~ treat) 
emmeans(bmod4, pairwise ~ res) 

#outlier in the control group removed
bmod5 <- lmer(log(m_flor) ~ res * treat + (1|random), data = vit_biof3) 
bmod5a <- lmer(log(m_flor) ~ res + treat + (1|random), data = vit_biof3) 
anova(bmod5, bmod5a) 
bmod5b <- lmer(log(m_flor) ~ treat + (1|random), data = vit_biof3) 
anova(bmod5a, bmod5b) 
bmod5c <- lmer(log(m_flor) ~ res + (1|random), data = vit_biof3) 
anova(bmod5a, bmod5c) 
bmod6 <- lm(log(m_flor) ~ 1, data = vit_biof3)
anova(bmod5c, bmod6) 
emmeans(bmod5c, pairwise ~ res) 

#In vivo analysis - bootstrapping and comparing confidence intervals

library(boot)

#boostrap function
meanfun <- function(data, i){
  d <- data[i,]
  return(mean(d))   
}

##subset data
#In vivo
vivo_biof <- filter(bio_means, lab == "vivo")
#ancestral
anc_biof <- filter(vivo_biof, treat == "anc")
#resistant
res_vivo <- filter(vivo_biof, res == "res")
#susceptible
susc_vivo <- filter(vivo_biof, res == "susc")
#post-phage therapy
post_pt <- filter(vivo_biof, treat == "T1" | treat == "T2")

#Bootstrap 1 - pre-phage therapy
data <- log(anc_biof$m_flor) #rename to data so function can work
data <- data.frame(data)
summary(data)
bo <- boot(data[, "data", drop = FALSE], statistic=meanfun, R=5000)
bo #mean = 5.53, SE = 0.038
boot.ci(bo, conf = 0.95, type = "norm") #CIs of ancestral biofilm growth = 5.45 - 5.6

#Bootstrap 2 - post phage therapy
data <- log(post_pt$m_flor)
data <- data.frame(data)
summary(data)
bo <- boot(data[, "data", drop = FALSE], statistic=meanfun, R=5000)
bo #mean = 5.388, SE = 0.06
boot.ci(bo, conf = 0.95, type = "norm") #5.272 - 5.507

#Bootstrap 3 - resistant
data <- log(res_vivo$m_flor)
data <- data.frame(data)
summary(data)
bo <- boot(data[, "data", drop = FALSE], statistic=meanfun, R=5000)
bo #mean = 5.38, SE = 0.058
boot.ci(bo, conf = 0.95, type = "norm") #CIs of resistant in vivo isolates = 5.27 - 5.498

#Bootstrap 3 - susceptible
data <- log(susc_vivo$m_flor)
data <- data.frame(data)
summary(data)
bo <- boot(data[, "data", drop = FALSE], statistic=meanfun, R=5000)
bo #mean = 5.39, SE = 0.151
boot.ci(bo, conf = 0.95, type = "norm") #CIs of susceptible in vivo isolates = 5.098 - 5.689

###Plots

vivo_biof2 <- vivo_biof
vivo_biof2$treat[vivo_biof2$treat == "T1"] <- "pt"
vivo_biof2$treat[vivo_biof2$treat == "T2"] <- "pt"

#add in which ancestor and in vivo isolates were phage resistant to one phage
#clone 12 of ancestor is single res (to 14-1)
vivo_biof2$res[12] <- "sing_res"
vivo_biof2$res[vivo_biof2$res == "anc"] <- "susc"
vivo_biof2$res[28] <- "sing_res"
vivo_biof2$res[30] <- "sing_res"

b1 <- ggplot() +
  geom_boxplot(data = vivo_biof2, aes(x = treat, y = log(m_flor))) +
  geom_point(data = vivo_biof2, aes(x = treat, y = log(m_flor), fill = res), size = 3, shape = 21, position = position_jitter(0.1)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14, colour = "black"), axis.title = element_text(size = 16), axis.text.x = element_text(size = 14, colour = "black"), legend.position = "bottom", plot.title = element_text(hjust = 0, size = 16), legend.text = element_text(size = 13), legend.title = element_text(size = 14)) +
  ylab("Biofilm production (fluorescence log(560/590nm))") +
  xlab("Treatment") +
  scale_x_discrete(labels = c("Pre-phage", "Post-phage")) +
  labs(fill = "Resistance") +
  scale_y_continuous(limits = c(4.8, 6.8)) +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 3, labels = c("Resistant\n(both phages)", "Resistant\n(one phage)", "Susceptible")) +
  ggtitle("(c)")

b1

b2 <- ggplot(vit_biof4, aes(x = treat, y = log(m_flor))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.1), size = 3, aes(fill = res), shape = 21)  +
  theme_bw() +
  ylab("Fluorescence log(560/590nm)") +
  xlab("Treatment") +
  theme(axis.text.y = element_blank(), axis.title = element_text(size = 16), axis.text.x = element_text(size = 14, colour = "black"), legend.text = element_text(size = 13), legend.title = element_text(size = 14), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0, size = 16)) +
  palettetown::scale_fill_poke(pokemon = 'articuno', spread = 3, labels = c("Resistant\n(both phages)", "Resistant\n(one phage)", "Susceptible")) +
  labs(fill = "Resistance") +
  scale_y_continuous(limits = c(4.8, 6.8)) +
  scale_x_discrete(labels = c("Control", "Phage added\nonce", "Phage added\nrepeatedly")) +
  ggtitle("(d)")

b2

biofs <- b1 + b2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
biofs

#Save part 2 of Figure 3 (panels c and d)

ggsave("plots/Figure_3_cd.pdf", biofs, height = 6, width = 9)
ggsave("plots/Figure_3_cd.png", biofs, height = 6, width = 9)
