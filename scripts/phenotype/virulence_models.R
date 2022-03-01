#---------------------------------------------------------------------#
# survival analyses of bacteria isolated from phage therapy trials ####
#---------------------------------------------------------------------#

# load in packages ####
library(survival)
library(rstanarm) # remotes::install_github('jburos/rstanarm@bugfix/stansurv-formula') - This may not create documentation for each function though
# installing the development branch of rstanarm is very difficult at the moment. But the amazing functionality means that we have persevered with it.
library(bayesplot)
library(tidyverse)
library(tidybayes)
library(patchwork)
library(palettetown)

# load in data
vir_dat <- read.csv("data/phenotype/virulence_data.csv", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)

#guide:
  #clone: bacterial clone ID
  #treat: whether bacterial clones were isolated from in vivo or in vitro (mix (phage added once), or mix_add (phage added repeatedly)). Time in vivo denoted by t1 (day 2) or t2 (day 4)
  #rep: similar to clonal ID but should be consistent across experimental datasets for virulence, growth rate and biofilm production
  #res: resistance to phage. susc = susceptible. sing_res = resistant to one phage. res = resistant to both phage
  #time_of_death: time at which a galleria is dead (hrs). If a galleria is not dead, it's time is displayed as time of final measurement with its status indicated as alive
  #status: 1 = dead, 0 = alive
  #galleria_ind1: continuous ID for every galleria in this experiment (total technical replicates)
  #galleria_ind2: 1 - 20 for each galleria used (technical replicates) for each bacterial clone

vir_dat$clone %>% unique()
unique(vir_dat$treat)

# split data into invitro and invivo treatments
d_invivo <- filter(vir_dat, treat %in% c('anc', 'vivo_t1', 'vivo_t2'))
d_invitro <- filter(vir_dat, treat %in% c('control_nop', 'mix', 'mix_add'))

#----------------#
# invivo data ####
#----------------#

# first name resistance based on the wells of the data
d_invivo <- mutate(d_invivo, res = case_when(clone == "E3" ~ "susceptible",
                                             clone == "E4" ~ "resistant",
                                             clone == "E5" ~ "resistant",
                                             clone == "E6" ~ "resistant",
                                             clone == "E7" ~ "susceptible",
                                             clone == "E8" ~ "resistant",
                                             clone == "E1" ~ "resistant",
                                             clone == "E2" ~ "resistant",
                                             clone == "F1" ~ "susceptible",
                                             clone == "F2" ~ "resistant",
                                             clone == "G1" ~ "ancestral",
                                             clone == "G2" ~ "ancestral",
                                             clone == "G3" ~ "ancestral",
                                             clone == "G4" ~ "ancestral",
                                             clone == "G5" ~ "ancestral",
                                             clone == "G6" ~ "ancestral",
                                             clone == "G7" ~ "ancestral",
                                             clone == "G8" ~ "ancestral",
                                             clone == "G9" ~ "ancestral",
                                             clone == "G10" ~ "ancestral",
                                             clone == "G11" ~ "ancestral",
                                             clone == "G12" ~ "ancestral",
                                             clone == "H1" ~ "ancestral",
                                             clone == "H2" ~ "ancestral",
                                             clone == "H3" ~ "ancestral",
                                             clone == "H4" ~ "ancestral",
                                             clone == "H5" ~ "ancestral",
                                             clone == "H6" ~ "ancestral",
                                             clone == "H7" ~ "ancestral",
                                             clone == "H8" ~ "ancestral",
                                             clone == "H9" ~ "ancestral",
                                             clone == "H10" ~ "ancestral",
                                             clone == "H11" ~ "ancestral",
                                             clone == "H12" ~ "ancestral"))

d_invivo <- mutate(d_invivo, res2 = case_when(res == 'ancestral' & clone != 'G12' ~ 'susceptible',
                                              clone == 'G12' ~ 'resistant',
                                              TRUE ~ res),
                   treat2 = ifelse(res == 'ancestral', 'before_phage_therapy', 'after_phage_therapy'))

# lets look at how the data is organised
select(d_invivo, clone, treat, res2) %>%
  distinct() %>%
  group_by(treat) %>%
  count(res2)

filter(d_invivo, treat2 == 'after_phage_therapy') %>%
  select(clone, time_of_death, status) %>%
  group_by(clone) %>%
  summarise(mean = mean(time_of_death),
            n = sum(status))

# run a parametric proportional hazards model
# follows instructions from "Bayesian Survival Analysis Using the rstanarm R Package, 2020, Brilleman et al. arXiv"
# run an m-splines model that models changes in the baseline hazard across time
# use proportional hazards because we dont expect any change in the difference of death rate to occur between groups apart from the treatment. no time varying effect.

# Adjust models so we only pre and post phage therapy as explanatory factors

mod_invivo <- stan_surv(Surv(time_of_death, status) ~ treat2 + (1|clone),
                  data = d_invivo,
                  chains = 3,
                  cores = 3,
                  seed = 42,
                  iter = 3000)

post_phage_therapy_clones <- filter(d_invivo, treat2 == 'after_phage_therapy') %>%
  pull(clone) %>%
  unique(.)

mod_invivo

# the hazards ratio is time independent
# it is the exp() of the estimated hazard, but is therefore linked to the intercept

summary(mod_invivo)

# get a list of the variables in the model
tidybayes::get_variables(mod_invivo)

to_plot <- tidybayes::get_variables(mod_invivo)[1:2]

ranefs <- str_subset(tidybayes::get_variables(mod_invivo), paste(post_phage_therapy_clones, collapse = '|'))

# check key mcmc trace plots
mcmc_trace(mod_invivo, pars = to_plot)
# fuzzy caterpillars

# extract draws
# the exp() of the difference between two factors is the hazards ratio
params_invivo <- spread_draws(mod_invivo, !!!syms(c(to_plot, ranefs))) %>%
  janitor::clean_names()

# calculate the hazard for each treatment by adding things onto the intercept as would be normal for a summary(model) table in R
params_invivo <- mutate(params_invivo, beforephage = intercept + treat2before_phage_therapy,
         postphage = intercept) %>%
  # calculate a hazard ratio for the postphage treatment
  mutate(hazard_post_vs_pre = exp(postphage - beforephage))

# calculate invivo hazard for each clone
invivo_hazards <- mutate(params_invivo, T1_3 = postphage + b_intercept_clone_e3,
                         T1_4= postphage + b_intercept_clone_e4,
                         T1_5 = postphage + b_intercept_clone_e5,
                         T1_6 = postphage + b_intercept_clone_e6,
                         T1_7 = postphage + b_intercept_clone_e7,
                         T1_8 = postphage + b_intercept_clone_e8,
                         T2_1 = postphage + b_intercept_clone_f1,
                         T2_2 = postphage + b_intercept_clone_f2) %>%
  select(chain, iteration, draw, postphage, starts_with('T', ignore.case = FALSE)) %>%
  pivot_longer(., cols = starts_with('T', ignore.case = FALSE), names_to = 'clone', values_to = 'val') %>%
  mutate(hazard_ratio = exp(val - postphage)) %>%
  group_by(clone) %>%
  median_qi(hazard_ratio)

# save out this file
write.csv(select(invivo_hazards, clone, virulence = hazard_ratio), 'data/phenotype/invivo_hazard_ratios.csv', row.names = FALSE)

# calculate credible intervals of hazard ratios
params2_invivo <- select(params_invivo, beforephage:hazard_post_vs_pre) %>%
  pivot_longer(cols = everything(), names_to = 'variable', values_to = 'estimate') %>%
  group_by(variable) %>%
  median_qi() %>%
  ungroup() %>%
  rename(treat = variable)

# hazard ratio is 0.0450, which means the clones isolated after phage therapy were 95% less likely to cause the death of a galleria at any time point during the study

cols <- c("dark grey", 'black')

# calculate survival curves
# predict over population-level estimates
d_preds <- select(d_invivo, treat2, clone) %>%
  distinct() %>%
  mutate(id = 1:n(),
         id2 = group_indices(., treat2),
         treat = treat2) %>%
  nest_legacy(-c(id2, treat)) %>%
  mutate(., preds = map(data, ~posterior_survfit(mod_invivo, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE, dynamic = TRUE)))

d_preds <- unnest(d_preds, preds) %>%
  select(-data) %>%
  mutate(treat = ifelse(treat == 'after_phage_therapy', 'After phage', 'Before phage'))

d_preds_ancest <- filter(d_preds, treat == 'Before phage')

vivo_plot <- ggplot(d_preds, aes(time, median, fill = treat)) +
  geom_line(aes(col = treat)) +
  geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), alpha = 0.2) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(axis.text.y = element_text(size = 18, colour = "black"), axis.title = element_text(size = 20), axis.text.x = element_text(size = 18, colour = "black"), strip.background = element_blank(), plot.title = element_text(size = 18, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 17), legend.title = element_blank()) +
  ylab("Survival probability") +
  xlab("") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ggtitle('(a)') +
  scale_x_continuous(n.breaks = 5) 

vivo_plot <- vivo_plot + geom_label(x = 22, y = 0.96, label = "Proportion that died: 0.99  \nMean time to death: 15 hours", hjust = 0, fill = 'white', label.size = NA, size = 5, col = 'white') +
  annotate(x = 22, y = 0.96, label = "Proportion that died: 0.99  \nMean time to death: 15 hours", size = 5, geom = 'text', hjust = 0) + 
  geom_label(x = 22, y = 0.84, label = "Proportion that died: 0.7  \nMean time to death: 20.5 hours", hjust = 0, fill = 'white', label.size = NA, size = 5, col = 'white') +
  annotate(x = 22, y = 0.84, label = "Proportion that died: 0.7  \nMean time to death: 20.5 hours", hjust = 0, col = 'dark grey', size = 5, geom = 'text')

vivo_plot
  
# calculate mean time of death for each group
group_by(d_invivo, treat2, status) %>%
  summarise(mean = mean(time_of_death),
            n = n())
479/480
112/(112+48)

#-----------------#
# invitro data ####
#-----------------#

# lets look at how the data is organised
select(d_invitro, clone, treat, res) %>%
  distinct() %>%
  group_by(treat) %>%
  count(res)

select(d_invitro, clone, treat, res, status) %>%
  group_by(treat, res) %>%
  count(status)

# try and run a model using rstanarm - this uses the survival branch of the rstanarm package
# run a parametric proportional hazards model
# follows instructions from "Bayesian Survival Analysis Using the rstanarm R Package, 2020, Brilleman et al. arXiv"

# run model for survival curve with random effect of clone and replicate
# uses the defaut m-splines model
mod_invitro <- stan_surv(Surv(time_of_death, status) ~ res * treat + (1|clone) + (1|rep),
                  data = d_invitro,
                  chains = 3,
                  cores = 3,
                  seed = 42,
                  iter = 3000)

mod_invitro
summary(mod_invitro)

# get a list of the variables in the model
tidybayes::get_variables(mod_invitro)

# plot parameters
to_plot <- tidybayes::get_variables(mod_invitro)[1:9]

# check key mcmc trace plots
mcmc_trace(mod_invitro, pars = to_plot)

# extract draws - these are the hazard ratios
params_invitro <- spread_draws(mod_invitro, !!!syms(to_plot)) %>%
  janitor::clean_names() %>%
  mutate(control_np.susc = intercept + ressusc,
         mix.resist = intercept + treatmix,
         mix.resist_sing = intercept + resresist_sing + treatmix + resresist_sing_treatmix,
         mix.susc = intercept + treatmix + ressusc + ressusc_treatmix,
         mix_add.resist = intercept + treatmix_add,
         mix_add.resist_sing = intercept + resresist_sing + treatmix_add + resresist_sing_treatmix_add,
         mix_add.susc = intercept + treatmix_add + ressusc + ressusc_treatmix_add,
         mix = (mix.susc + mix.resist + mix.resist_sing)/3,
         mix_add = (mix_add.susc + mix_add.resist + mix_add.resist_sing)/3,
         susc = (control_np.susc + mix.susc + mix_add.susc)/3,
         resist_single = (mix.resist_sing + mix_add.resist_sing)/2,
         resist = (mix_add.resist + mix.resist)/2)

# the exp() of the difference between two factors is the hazards ratio
# calculate the hazard for each treatment by adding things onto the intercept as would be normal for a summary(model) table in R
params_invitro <- mutate(params_invitro, hazard_mix_vs_control = exp(mix - control_np.susc),
         hazard_mixadd_vs_control = exp(mix_add - control_np.susc),
         hazard_mixadd_vs_mix = exp(mix_add - mix),
         hazard_resist_susc = exp(resist - susc),
         hazard_resist_resist_one = exp(resist - resist_single),
         hazard_resist_one_susc = exp(resist_single - susc),
         hazard_suscept_mix_control = exp(mix.susc - control_np.susc),
         hazard_suscept_mix_add_control = exp(mix_add.susc - control_np.susc),
         hazard_suscept_mix_add_mix = exp(mix_add.susc - mix.susc))

# calculate credible intervals of hazard ratios
params2_invitro <- select(params_invitro, contains('hazard')) %>%
  pivot_longer(cols = everything(), names_to = 'variable', values_to = 'estimate') %>%
  group_by(variable) %>%
  median_qi() %>%
  ungroup() %>%
  rename(treat = variable)

# calculate survival curves
# predict over population-level estimates

# standardise over treatments only
d_preds <- select(d_invitro, res, treat, clone, rep) %>%
  distinct() %>%
  mutate(id = 1:n(),
         id2 = group_indices(., res, treat),
         res2 = res, 
         treat2 = treat) %>%
  nest_legacy(-c(treat2)) %>%
  mutate(., preds = map(data, ~posterior_survfit(mod_invitro, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE)))

d_preds <- unnest(d_preds, preds) %>%
  select(-data)

# standardise over resistance levels only
d_preds2 <- select(d_invitro, res, treat, clone, rep) %>%
  distinct() %>%
  mutate(id = 1:n(),
         id2 = group_indices(., res, treat),
         res2 = res, 
         treat2 = treat) %>%
  nest_legacy(-c(res2)) %>%
  mutate(., preds = map(data, ~posterior_survfit(mod_invitro, newdata = .x, times = 0, standardise = TRUE, extrapolate = TRUE)))

d_preds2 <- unnest(d_preds2, preds) %>%
  select(-data)

# plot treatment only effect

cols <- c('#fdcc8a', '#fc8d59', '#d7301f')

invitro <- ggplot(d_preds, aes(time, median, fill = treat2)) +
  geom_label(x = 0, y = 0.3, label = "Proportion that died: 0.8  \nMean time to death: 19 hours", hjust = 0, fill = "white", label.size = NA, size = 5, col = "white") + 
  annotate(x = 0, y = 0.3, label = "Proportion that died: 0.8  \nMean time to death: 19 hours", hjust = 0, size = 5, col = "#fdcc8a", geom = 'text') +
  geom_label(x = 0, y = 0.175, label = "Proportion that died: 0.33  \nMean time to death: 36.7 hours", hjust = 0, fill = "white", label.size = NA, size = 5, col = "white") +
  annotate(x = 0, y = 0.175, label = "Proportion that died: 0.33  \nMean time to death: 36.7 hours", hjust = 0, size = 5, col = "#fc8d59", geom = 'text') +
  geom_label(x = 0, y = 0.05, label = "Proportion that died: 0.41  \nMean time to death: 35.4 hours", hjust = 0, fill = 'white', label.size = NA, size = 5, col = "white") +
  annotate(x = 0, y = 0.05, label = "Proportion that died: 0.41  \nMean time to death: 35.4 hours", hjust = 0, size = 5, col = "#d7301f", geom = 'text') +
  geom_line(aes(col = treat2)) +
  geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), alpha = 0.2) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title = element_text(size = 20), axis.text.x = element_text(size = 18, colour = "black"), strip.background = element_blank(), plot.title = element_text(size = 18, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 17), legend.title = element_blank()) +
  ylab("Survival probability") +
  xlab("Time since infection (hours)") +
  ggtitle('(b)') +
  scale_color_manual(values = cols, labels = c('Control', 'Phage added\nonce', 'Phage added\nrepeatedly')) +
  scale_fill_manual(values = cols, labels = c('Control', 'Phage added\nonce', 'Phage added\nrepeatedly'))

invitro2 <- ggplot(d_preds2, aes(time, median, fill = res2)) +
  geom_label(x = 0, y = 0.3, label = "Proportion that died: 0.08  \nMean time to death: 35 hours", hjust = 0, fill = "white", label.size = NA, size = 5, col = "white") + 
  annotate(x = 0, y = 0.3, label = "Proportion that died: 0.08  \nMean time to death: 35 hours", hjust = 0,  size = 5, col = "#6090C8", geom = 'text') +
  geom_label(x = 0, y = 0.175, label = "Proportion that died: 0.47  \nMean time to death: 24.4 hours", hjust = 0, fill = "white", label.size = NA, size = 5, col = "white") +
  annotate(x = 0, y = 0.175, label = "Proportion that died: 0.47  \nMean time to death: 24.4 hours", hjust = 0,  size = 5, col = "#404058", geom = 'text') +
  geom_label(x = 0, y = 0.05, label = "Proportion that died: 0.65  \nMean time to death: 32.5 hours", hjust = 0, fill = "white", label.size = NA, size = 5, col = "white") +
  annotate(x = 0, y = 0.05, label = "Proportion that died: 0.65  \nMean time to death: 32.5 hours", hjust = 0,  size = 5, col = "#B82838", geom = 'text') +
  geom_line(aes(col = res2)) +
  geom_ribbon(aes(time, ymin = ci_lb, ymax = ci_ub), alpha = 0.2) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title = element_text(size = 20), axis.text.x = element_text(size = 18, colour = "black"), strip.background = element_blank(), plot.title = element_text(size = 18, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 17)) +
  ylab("Survival probability") +
  xlab("") +
  ggtitle('(c)') +
  scale_color_poke('', pokemon = 'articuno', spread = 3, labels = c("Resistant\n(both phages)", "Resistant\n(one phage)", "Susceptible")) +
  scale_fill_poke('', pokemon = 'articuno', spread = 3, labels = c("Resistant\n(both phages)", "Resistant\n(one phage)", "Susceptible"))

# calculate mean time of death for each group
group_by(d_invitro, treat, status) %>%
  summarise(mean = mean(time_of_death),
            n = n())
96/(96+24)
117/(117+243)
107/(107+153)

group_by(d_invitro, res, status) %>%
  summarise(mean = mean(time_of_death),
            n = n())

total_plot <- vivo_plot + invitro + invitro2

ggsave('plots/Figure_2.pdf', total_plot, height = 7, width = 18)
ggsave('plots/Figure_2.png', total_plot, height = 7, width = 18)
