library(lme4)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(emmeans)

##In vivo timeseries

viv_ts <- read.csv("data/phenotype/invivo_timeseries.csv", header = T)

#this data pertains to a time-series assay in which in vivo isolates from time points 1 and 2 are assayed against phage from timepoints 1, 2 and 3. Phage at timepoint 3 are post phage therapy. Isolates also tested against ancestral phage individually (PVM and 14_1)
#1 = bacteria susceptible
#0 = bacteria resistant

viv_ts_props <- group_by(viv_ts, bact_time, phage_time, phage_type) %>%
  summarise(., sum_inf = sum(infect),
            pop = length(bact_clone),
            prop = sum_inf/pop,
            sd_p = sd(prop)) %>%
  ungroup()

viv_ts_props_t1 <- filter(viv_ts_props, bact_time == "T1")

viv_ts_props_t2 <- filter(viv_ts_props, bact_time == "T2")

viv_ts_props <- merge(viv_ts_props_t1, viv_ts_props_t2, all = T)

##plot

anc_res <- read.csv("data/phenotype/anc_resistance.csv", header = T)

labs <- c(expression('14-1 (anc.)'), expression('PNM (anc.)'), "T1", "T2", "T3")

viv <- ggplot(data = viv_ts_props, aes(x = phage_time, y = prop, group = bact_time)) +
  geom_point(data = viv_ts_props, aes(x = phage_time, y = prop, group = bact_time, col = bact_time), size = 2, position = position_dodge(0.3)) +
  geom_point(data = anc_res, aes(x = phage_time, y = prop, group = bact_time), size = 2) +
  geom_line(data = viv_ts_props, position = position_dodge(0.3), aes(col = bact_time, group = bact_time)) +
  theme_bw() +
  scale_alpha(guide = "none") +
  labs(col = "Clone") +
  theme(plot.title = element_text(size = 15), axis.text = element_text(size = 13, colour = "black"), axis.title = element_text(size = 16), axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none") +
  scale_x_discrete(labels = labs)  +
  labs(title = expression("(a)"~italic("In vivo"))) +
  palettetown::scale_color_poke(pokemon = 'fearow', spread = 3)  +
  geom_line(data = anc_res, position = position_dodge(0.3), col = "black") +
  geom_label(label = "Anc.", y = 1, x = 5.2)

viv #add this to the plot in the in vitro time series script

#Analysis

vivo_only <- filter(viv_ts, phage_type == "vivo")
vivo_t1 <- filter(vivo_only, bact_time == "T1")

viv_tsm2.1 <- glmer(infect ~ phage_time + (1|bact_clone), data = vivo_t1, family = binomial)
RVAideMemoire::overdisp.glmer(viv_tsm2.1) 

viv_tsm2.2 <- glmer(infect ~ 1 + (1|bact_clone), data = vivo_t1, family = binomial)
anova(viv_tsm2.1, viv_tsm2.2) #sig effect of phage time
emmeans(viv_tsm2.1, pairwise ~ phage_time) #all non-sig
