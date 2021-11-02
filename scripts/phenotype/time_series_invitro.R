##In vitro and in vivo time series assay 
library(lme4)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(emmeans)
library(patchwork)

vit_ts <- read.csv("data/phenotype/invitro_timeseries.csv", header = T)

#guide:
#treat = in vitro phage treatment. 
  #mix = phage added once
  #mix_add = phage added repeatedly
#rep = biological replicate within each treatment
#bact_clone = individual clones tested within each treatment replicate (a biological replicate nested within another biological replicate)
#bact_time = time point of bacterial evolution. 1 - 3 indicates days 4, 8 and 12 of coevolution respectively.
#phage_time = time point of phage evolution. 1 - 3 indicates days 4, 8 and 12 of coevolution respectively.
#infect. Binary dataframe indicating a successful infection of the phage to the bacteria (1) or resistance to phage (0)

#ensuring variables are read as characters
vit_ts$rep <- as.character(vit_ts$rep)
vit_ts$bact_time <- as.character(vit_ts$bact_time)
vit_ts$phage_time <- as.character(vit_ts$phage_time)

#need to gather proportions of infectivity across biological replicates
vit_ts_props <- group_by(vit_ts, treat, rep, bact_time, phage_time) %>%
  summarise(., sum_inf = sum(infect),
            pop = n(),
            prop = sum_inf/pop,
            prop_resist = (pop - sum_inf)/pop)

#tidy dataframe. This is just to generate summary statistics
vit_ts_props2 <- vit_ts_props[c(1:108), c(1:4,7)]

vit_ts_props2 <- na.omit(vit_ts_props2)

#Take mean proportion of bacteria infected across replicates within each treatment and bacteria / phage time point. This is used for generating plots later
vit_ts_means <- group_by(vit_ts_props2, treat, bact_time, phage_time) %>%
  summarise(., m.prop = mean(prop),
            sdp = sd(prop),
            sq = sqrt(6),
            se = sdp/sq)

#models
summary(vit_ts_props) #dataset used for analysis

#binomial model with success/fails - analysis used in manuscript - separate models for mix and mix_add because the model can't converge with everything (treatment, phage and bacteria time).
vit_mix_means <- filter(vit_ts_props, treat == "mix") 

inf <- cbind(vit_mix_means$sum_inf, vit_mix_means$pop - vit_mix_means$sum_inf) #successful infections / failed infections

m1.1 <- glmer(inf ~ bact_time * phage_time + (1|rep), family = binomial, data = vit_mix_means) 
m1.2 <- glmer(inf ~ bact_time + phage_time + (1|rep), family = binomial, data = vit_mix_means) 
anova(m1.1, m1.2) #no interaction
m1.3 <- glmer(inf ~ phage_time + (1|rep), family = binomial, data = vit_mix_means)
anova(m1.3, m1.2) #sig effect of bacterial time
m1.4 <- glmer(inf ~ bact_time + (1|rep), family = binomial, data = vit_mix_means)
anova(m1.2, m1.4) #no effect of phage time
m1.5 <- glmer(inf ~ 1 + (1|rep), family = binomial, data = vit_mix_means)
anova(m1.5, m1.4)
emmeans(m1.4, pairwise ~ bact_time, type = "response")

#mix add treatment

#binomial model with success/fails
vit_mixadd_means <- filter(vit_ts_props, treat == "mix_add")

inf <- cbind(vit_mixadd_means$sum_inf, vit_mixadd_means$pop - vit_mixadd_means$sum_inf) #successful infections / failed infections

m1.1 <- glmer(inf ~ bact_time * phage_time + (1|rep), family = binomial, data = vit_mixadd_means) 
m1.2 <- glmer(inf ~ bact_time + phage_time + (1|rep), family = binomial, data = vit_mixadd_means)
anova(m1.1, m1.2) #no interaction
m1.3 <- glmer(inf ~ phage_time + (1|rep), family = binomial, data = vit_mixadd_means)
anova(m1.3, m1.2) #sig effect of bacterial time
m1.4 <- glmer(inf ~ bact_time + (1|rep), family = binomial, data = vit_mixadd_means)
anova(m1.2, m1.4) #no effect of phage time
m1.5 <- glmer(inf ~ 1 + (1|rep), family = binomial, data = vit_mixadd_means)
anova(m1.5, m1.4)

emmeans(m1.4, pairwise ~ bact_time, type = "response")

#Plot - Figure 1

labs2 <- c("Once", "Continuously")

cols <- c('#377eb8', '#2AA747', '#9139EE')

vit_ts_means$treat[vit_ts_means$treat == "mix"] <- "(b) Once"
vit_ts_means$treat[vit_ts_means$treat == "mix_add"] <- "(c) Repeatedly"
vit_ts_means$bact_time[vit_ts_means$bact_time == "1"] <- "T1"
vit_ts_means$bact_time[vit_ts_means$bact_time == "2"] <- "T2"
vit_ts_means$bact_time[vit_ts_means$bact_time == "3"] <- "T3"
vit_ts_means$phage_time[vit_ts_means$phage_time == "1"] <- "T1"
vit_ts_means$phage_time[vit_ts_means$phage_time == "2"] <- "T2"
vit_ts_means$phage_time[vit_ts_means$phage_time == "3"] <- "T3"

names(vit_ts_means)[3] <- "Phage time"

vit_ts_props2 <- vit_ts_props
vit_ts_props2$bact_time[vit_ts_props2$bact_time == "1"] <- "T1"
vit_ts_props2$bact_time[vit_ts_props2$bact_time == "2"] <- "T2"
vit_ts_props2$bact_time[vit_ts_props2$bact_time == "3"] <- "T3"
vit_ts_props2$phage_time[vit_ts_props2$phage_time == "1"] <- "T1"
vit_ts_props2$phage_time[vit_ts_props2$phage_time == "2"] <- "T2"
vit_ts_props2$phage_time[vit_ts_props2$phage_time == "3"] <- "T3"

label_facets <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') ', string, sep = '')
  return(string)
}


vit_ts_props3 <- na.omit(vit_ts_props2)

names(vit_ts_props3)[4] <- "Phage time"

vit_ts_props3$treat <- as.character(vit_ts_props3$treat)
vit_ts_props3$treat[vit_ts_props3$treat == "mix"] <- "(b) Once"
vit_ts_props3$treat[vit_ts_props3$treat == "mix_add"] <- "(c) Repeatedly"


#Integrating resistance of in vitro bacteria (t3) to ancestral phage

vit_anc <- read.csv("phage_therapy/time_series_exp/data/end_resistance.csv", header = T)

vit_anc <- na.omit(vit_anc)

vit_anc_m <- group_by(vit_anc, treat, rep, phage) %>%
  summarise(s_infect = sum(infect),
            total = length(infect),
            prop = s_infect/total,
            se_inf = sd(infect)/sqrt(total))

vit_anc_m2 <- group_by(vit_anc, treat, phage) %>%
  summarise(s_infect = sum(infect),
            total = length(infect),
            m.prop = s_infect/total,
            se_inf = sd(infect)/sqrt(total))

vit_anc_m$phage[vit_anc_m$phage == "14_1"] <- "a14-1"
vit_anc_m$phage[vit_anc_m$phage == "PNM"] <- "aPNM"
vit_anc_m2$phage[vit_anc_m2$phage == "14_1"] <- "a14-1"
vit_anc_m2$phage[vit_anc_m2$phage == "PNM"] <- "aPNM"

names(vit_anc_m)[3] <- "phage_time"
names(vit_anc_m2)[2] <- "phage_time"
vit_anc_m <- mutate(vit_anc_m, bact_time = "T3")
vit_anc_m2 <- mutate(vit_anc_m2, bact_time = "T3")
vit_anc_m$treat[vit_anc_m$treat == "mix"] <- "(b) Once"
vit_anc_m2$treat[vit_anc_m2$treat == "mix"] <- "(b) Once"
vit_anc_m$treat[vit_anc_m$treat == "mix_add"] <- "(c) Repeatedly"
vit_anc_m2$treat[vit_anc_m2$treat == "mix_add"] <- "(c) Repeatedly"

vit <- ggplot(vit_ts_means, aes(group = bact_time, y = m.prop, x = `Phage time`)) +
  geom_point(data = vit_ts_props3, aes(x = `Phage time`, y = prop,  fill = bact_time, col = bact_time), position = position_jitterdodge(0.3), alpha = 0.3) +
  geom_point(position = position_dodge(0.8), aes(col = bact_time), size = 2) +
  geom_errorbar(aes(ymin = m.prop-se, ymax = m.prop+se, col = bact_time), width = 0.3, position = position_dodge(0.8)) +
  theme_bw() +
  facet_wrap(~treat) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14, hjust = 0), legend.text = element_text(size = 14), legend.position = "bottom", axis.text = element_text(size = 14, colour = "black"), axis.title = element_text(size = 16), legend.title = element_text(size = 16), axis.title.y = element_text(hjust = -0.2), plot.title = element_text(size = 15)) +
  annotate("rect", xmin = 0.7, xmax = 2.2, ymin = 0, ymax = 1, alpha = .1) + 
  scale_alpha(guide = 'none') +
  ylab("Proportion of bacterial clones susceptible to phage") +
  labs(fill = "Bacteria time", col = "Bacteria time", title = expression(italic("In vitro")~" phage added:")) +
  palettetown::scale_color_poke(pokemon = 'fearow', spread = 3) +
  geom_point(data = vit_anc_m, aes(x = phage_time, y = prop, col = bact_time), position = position_jitterdodge(0.3), alpha = 0.3) +
  geom_errorbar(data = vit_anc_m2, aes(x = phage_time, ymin = m.prop-se_inf, ymax = m.prop+se_inf, col = bact_time), position = position_dodge(0.3), width = 0.2) +
  geom_point(data = vit_anc_m2, aes(x = phage_time, y = m.prop, col = bact_time), position = position_dodge(0.3), size = 2) +
  geom_line(position = position_dodge(0.8), aes(col = bact_time)) +
  scale_x_discrete(labels = c('14-1\n(anc.)', 'PNM\n(anc.)', "T1", "T2", "T3")) 

vit

#running viv from time_series_invivo script on another tab to import plot object (run on same RStudio console then come back to run this line)
vivvit <- viv + vit + plot_layout(ncol = 1, heights = c(1,2))

ggsave("phage_therapy/Plots/Figure_1.pdf", vivvit, height = 8, width = 7)
ggsave("phage_therapy/Plots/Figure_1.png", vivvit, height = 8, width = 7)

## analysing probability of resistance of T3 bacteria in added once and repeatedly treatment to ancestral vs contemporary bacteria

vit_t3s <- filter(vit_ts_props, bact_time == 3 & phage_time == 3)

vit_t3s <- vit_t3s[,c(1:6)]

names(vit_t3s)[5] <- 's_infect'
names(vit_t3s)[6] <- 'total'

vit_anc3 <- vit_anc_m[, c(1:5, 8)]
vit_anc3$treat[vit_anc3$treat == '(b) Once'] <- "mix"
vit_anc3$treat[vit_anc3$treat == '(c) Repeatedly'] <- "mix_add"

t3s <- merge(vit_anc3, vit_t3s, all = T)

#analyse
infect <- cbind(t3s$s_infect, t3s$total - t3s$s_infect)

mod1 <- glmer(infect ~ treat * phage_time + (1|rep), data = t3s, family = binomial)
mod2 <- glmer(infect ~ treat + phage_time + (1|rep), data = t3s, family = binomial)
anova(mod1, mod2) #sig interaction
emmeans::emmeans(mod1, pairwise ~ phage_time | treat)
emmeans::emmeans(mod1, pairwise ~ treat | phage_time)
