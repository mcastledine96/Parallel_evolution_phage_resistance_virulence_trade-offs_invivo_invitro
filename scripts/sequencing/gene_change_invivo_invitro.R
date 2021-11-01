#----------------------------------------------------#
# Processing in vivo clone vcf files and plotting ####
#----------------------------------------------------#

# summary of what this script does

# 1. reads in clone in vivo vcf files and processes them
# 2. makes Figure 4
# 3. looks at distribution of SNPs/indels in in vivo clones and looks for patterns across phenotypes

# load packages
library(tidyverse)
library(vcfR)
library(palettetown)
library(ggvegan)
library(patchwork)

# function to clean up freebayes vcf files

tidy_freebayes <- function(freebayes_vcf){
  temp <- vcfR::read.vcfR(freebayes_vcf, verbose = FALSE)
  temp <- vcfR::vcfR2tidy(temp, single_frame = TRUE) %>%
    .$dat %>%
    janitor::clean_names() %>%
    dplyr::filter(dp > mean(dp) - 2*sd(dp) & dp < mean(dp) + 2*sd(dp)) %>%
    dplyr::select(., -c(id, filter, ns))
  temp$file <- tools::file_path_sans_ext(basename(freebayes_vcf))
  return(temp)
}

#-------------------------------------#
# first process the invivo samples ####
#-------------------------------------#

files_all <- list.files('data/sequencing/vcf/clone', full.names = TRUE, pattern = '.vcf')

# map ancestral reads
d_invivo <- files_all %>%
  map_df(., tidy_freebayes)

# split up file name
d_invivo <- separate(d_invivo, file, c('context', 'clone', 'blah', 'mapping_tool', 'ref_genome', 'blah1', 'blah2', 'blah3'), sep = '_')

# remove some columns to make the dataframe more manageable
d_invivo <- select(d_invivo, -c(blah1, blah2, blah3, ro, ao, pro, pao, qr, qa, pqr, pqa, srf, srr, sar, srp, ab, abp, run, rpp, rppr, rpl, rpr, eppr, dpra, odds, gti, cigar))

# filter SNPs and indels based on depth, quality, and other metrics
# recommended by the developer of freebayes
d2 <- filter(d_invivo, dp > 2) %>%
  filter(saf > 20) %>%
  filter(sap > 20) %>%
  filter(epp > 20) %>%
  filter(qual / dp > 2)

# read in data to map clone number to time point
d_meta <- read.csv('data/sequencing/info_for_renaming_ena_files.csv') %>%
  filter(context == 'invivo') %>%
  mutate(clone = parse_number(sample_name) %>% as.character()) %>%
  select(clone, time_point)

# add this info to d2
d2 <- left_join(d2, d_meta)

# make a bunch of key columns numeric and choose only the columns that are needed
# keep only variants that are at fixation
d2 <- select(d2, clone, time_point, mapping_tool, ref_genome, chrom, pos, ref, alt, qual, dp, af, numalt) %>%
  mutate(., across(qual:numalt, as.numeric)) %>%
  filter(af == 1)

# load in 24 ancestral clone data
d_ancestors <- readRDS('data/sequencing/ancestors.rds') %>%
  rename(ancest_alt = alt, ancest_prop = af_mean) %>%
  select(-ref)

d_ancest2 <- group_by(d_ancestors, pos) %>%
  filter(n() > 1) %>%
  mutate(num = 1:n()) %>%
  select(-ancest_prop) %>%
  pivot_wider(names_from = num, values_from = ancest_alt) %>%
  unite(ancest_alt, 2:ncol(.), sep = '_', na.rm = TRUE)

d_ancestors <- filter(d_ancestors, ! pos %in% d_ancest2$pos) %>%
  select(-ancest_prop) %>%
  bind_rows(., d_ancest2)

d2 <- left_join(d2, d_ancestors) %>%
  mutate(., ancest_pres = ifelse(is.na(ancest_alt), 'no', 'yes'))

# look at number of shared genetic variants
d_sum <- group_by(d2, alt, pos, ancest_pres) %>%
  tally()

# look at distribution of variants across all the clones
ggplot(d_sum, aes(n)) +
  geom_histogram(fill = 'white', col = 'black', bins = length(unique(d_sum$n))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(n.breaks = 10, limits = c(0,10)) +
  labs(y = 'count', x = 'number of clones',
       title = 'How many clones does each genetic variant occur in?') +
  facet_wrap(~ ancest_pres)

#----------------------------------------#
# load in SNPs from the invitro study ####
#----------------------------------------#

# create data frame of unique SNPs from in vitro data
invitro_snps <- read.csv('data/sequencing/invitro_snps_indels.csv') %>%
  filter(!is.na(ref)) %>%
  dplyr::select(-n) %>%
  select(., gene_name, ref, alt, id3, pos, change2) %>%
  distinct() %>%
  rename(ref_invitro = ref, alt_invitro = alt)

# load in invitro data
d_invitro <- read.csv('data/sequencing/invitro_snps_indels.csv', header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(-n) %>%
  select(., af, gene_name, id3, change2, id2, pos) %>%
  mutate(id4 = 'invitro populations',
         id2 = 1+id2)

# add ancestors to dataset
ancestor <- read.csv('data/sequencing/invitro_snps_indels.csv', header = TRUE, stringsAsFactors = FALSE) %>%
  select(., ancestral_prop, pos, gene_name, ref, alt, id3) %>%
  distinct(., .keep_all = TRUE) %>%
  pivot_longer(., ancestral_prop, names_to = 'treat', values_to = 'af') %>%
  mutate(., id4 = 'invitro populations',
         id2 = 1)

d_invitro <- bind_rows(d_invitro, ancestor)

# search for these SNPs in the in vivo clones
# create a column for where the in vivo clones will be in the plot
d3 <- filter(d2, pos %in% invitro_snps$pos) %>%
  left_join(invitro_snps) %>%
  group_by(time_point, clone) %>%
  mutate(id2 = cur_group_id() + 19,
         id4 = 'invivo clones') %>%
  ungroup() %>%
  select(., af, gene_name, id3, change2, id2, pos, id4)

# calculate proportion of each SNP in the in vivo clones as a whole (10 clones)
d3 <- d3 %>%
  group_by(gene_name, id3, change2, pos, id4) %>%
  summarise(af = sum(af)/10, .groups = 'drop') %>%
  mutate(id2 = 20)

d3 <- bind_rows(d_invitro, d3)

#------------------#
# make Figure 4 ####
#------------------#

# check number of SNPs/indels is right
unique(d3$id3) %>% length()

# create data from for where lines are
gene_lines <- dplyr::select(d3, gene_name, pos, id3) %>%
  distinct(.keep_all = TRUE) %>%
  arrange(desc(gene_name)) %>%
  mutate(., n = 1:n())

gene_lines2 <- gene_lines %>%
  group_by(gene_name) %>%
  top_n(-1) %>%
  ungroup()

d3 <- left_join(d3, gene_lines)

d4 <- group_by(d3, id2, id4) %>%
  do(data.frame(pos = unique(d3$pos))) %>%
  ungroup() %>%
  left_join(., gene_lines)

d4 <- left_join(d4, d3) %>%
  mutate(af = ifelse(is.na(af), 0, af),
         change2 = ifelse(is.na(change2), 'existing', change2))

xlines <- tibble(x = c(1.5, 7.5, 13.5, 19.5), id4 = c('invitro populations', 'invitro populations', 'invitro populations', 'invivo clones'))

head(d4)

# create plot three times to get the legends as we want

p1 <- ggplot(d4, aes(x = id2, y = interaction(gene_name, n))) +
  geom_tile(aes(alpha = af, fill = change2), col = 'black', show.legend = FALSE) +
  ylab("Gene") +
  theme_bw(base_size = 14) +
  facet_grid(. ~id4, scales = 'free_x', drop = TRUE, space = 'free_x') +
  scale_y_discrete(labels = function(x)(gsub('\\..*', '', x)), expand = c(0,0)) +
  scale_x_continuous(breaks = c(1, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 20), labels = c('Ancestor', '', 'Control', '', 'Phage added\nonce', '', 'Phage added\nrepeatedly', 'In vivo\nclones'), name = 'Treatment', expand = c(0, 0)) +
  scale_alpha_continuous(range = c(0,1),
                         guide = guide_legend(override.aes = list(fill = "#7090A0"))) +
  geom_vline(aes(xintercept = x), size = 1, xlines) +
  geom_hline(aes(yintercept = n - 0.5), size = 0.75, gene_lines2) +
  scale_fill_manual(values = c(rev(ichooseyou('vileplume', spread = 2))), labels = c("Existing", "Gain")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=3, fill = NA),
        panel.ontop = TRUE,
        plot.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  labs(fill = "Change type", alpha = "Frequency")

p2 <- ggplot(d4, aes(x = id2, y = interaction(gene_name, n))) +
  geom_tile(aes(alpha = af, fill = change2), col = 'black') +
  ylab("Gene") +
  theme_bw(base_size = 14) +
  facet_grid(. ~id4, scales = 'free_x', drop = TRUE, space = 'free_x') +
  scale_y_discrete(labels = function(x)(gsub('\\..*', '', x)), expand = c(0,0)) +
  scale_x_continuous(breaks = c(1, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 20), labels = c('Ancestor', '', 'Control', '', 'Phage added\nonce', '', 'Phage added\nrepeatedly', 'In vivo\nclones'), name = 'Treatment', expand = c(0, 0)) +
  scale_alpha_continuous(range = c(0,1),
                         guide = guide_legend(override.aes = list(fill = "#7090A0"),
                                              title.hjust = 0.4)) +
  geom_vline(aes(xintercept = x), size = 1, xlines) +
  geom_hline(aes(yintercept = n - 0.5), size = 0.75, gene_lines2) +
  scale_fill_manual(values = c(rev(ichooseyou('vileplume', spread = 2))), labels = c("Existing", "Gain"), guide = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=3, fill = NA),
        panel.ontop = TRUE,
        plot.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.justification = 'center') +
  labs(alpha = "Ancestral/Existing\nvariant frequency")

p3 <- ggplot(d4, aes(x = id2, y = interaction(gene_name, n))) +
  geom_tile(aes(alpha = af, fill = change2), col = 'black') +
  ylab("Gene") +
  theme_bw(base_size = 14) +
  facet_grid(. ~id4, scales = 'free_x', drop = TRUE, space = 'free_x') +
  scale_y_discrete(labels = function(x)(gsub('\\..*', '', x)), expand = c(0,0)) +
  scale_x_continuous(breaks = c(1, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 20), labels = c('Ancestor', '', 'Control', '', 'Phage added\nonce', '', 'Phage added\nrepeatedly', 'In vivo\nclones'), name = 'Treatment', expand = c(0, 0)) +
  scale_alpha_continuous(range = c(0,1),
                         guide = guide_legend(override.aes = list(fill = "#F87850"),
                                              title.hjust = 0.4)) +
  geom_vline(aes(xintercept = x), size = 1, xlines) +
  geom_hline(aes(yintercept = n - 0.5), size = 0.75, gene_lines2) +
  scale_fill_manual(values = c(rev(ichooseyou('vileplume', spread = 2))), labels = c("Existing", "Gain"), guide = NULL) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=3, fill = NA),
        panel.ontop = TRUE,
        plot.title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.justification = 'center') +
  labs(alpha = "New variant\nfrequency")

legend1 <- cowplot::get_legend(p2)
legend2 <- cowplot::get_legend(p3)
legends <- cowplot::plot_grid(legend1, legend2, ncol = 1, align = 'v')

p1 +  (plot_spacer()/legends + plot_layout(heights = c(0.6, 0.4))) + plot_layout(widths = c(0.8, 0.2))

ggsave('plots/Figure_4.pdf', last_plot(), height = 10, width = 8.5)
ggsave('plots/Figure_4.png', last_plot(), height = 10, width = 8.5)