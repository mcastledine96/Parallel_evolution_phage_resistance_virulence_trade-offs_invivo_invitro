#--------------------------------------#
# analysis of invitro genomics data ####
#--------------------------------------#

# summary of what this script does


# clear workspace (ideally should start a new R session at the beginning of each script)
rm(list = ls())

# load packages
library(tidyverse)
library(vegan)
library(ggvegan)
library(MicrobioUoE) # remotes::install_github('padpadpadpad/MicrobioUoE')
library(patchwork)
library(palettetown)
# also need concaveman, ggrepel, and ggforce 

# function to assign a SNP pos to its position in the genome
# returns the downstream, current, and upstream gene, and the distance to those
get_gene_info <- function(SNP_pos, d_genes, types_to_keep = c('gene'), pos_to_keep = c('downstream', 'current_pos', 'upstream')){
  
  # make numeric just in case it is not - needed for use of tidyr and nest()
  SNP_pos <- as.numeric(SNP_pos)
  
  # keep just gene names in d_genes
  d_genes <- dplyr::filter(d_genes, type %in% types_to_keep)
  
  # get gene name
  temp <- dplyr::filter(d_genes, start <= SNP_pos & stop >= SNP_pos) %>%
    dplyr::pull(., name)
  
  # if it is not there then return NA
  if(length(temp) == 0) temp <- NA
  
  # if is not NA, delete gene name from d_genes
  if(length(temp) > 0){d_genes <- dplyr::filter(d_genes, ! name %in% temp)}
  
  # some repeat regions overlap so just put return the first output of this
  # I dont think any of the gene positions overlap
  if(length(temp > 1)) temp <- temp[1]
  
  # get downstream gene - later position
  d_down <- filter(d_genes, start >= SNP_pos) %>%
    arrange(., start) %>%
    slice(seq_len(1))
  if(nrow(d_down) == 0){
    d_down <- data.frame(name = NA,
                         start = SNP_pos)
  }
  
  # get upstream_gene - earlier position
  d_up <- filter(d_genes, stop <= SNP_pos) %>%
    arrange(., desc(start)) %>%
    slice(seq_len(1))
  if(nrow(d_up) == 0){
    d_up <- data.frame(name = NA,
                       stop = SNP_pos)
  }
  
  # get upstream pos
  d_temp <- data.frame(gene_type = c('current_pos', 'downstream', 'upstream'),
                       gene_name = c(temp, d_down$name, d_up$name),
                       distance = c(0, abs(d_down$start - SNP_pos), abs(d_up$stop - SNP_pos)))
  
  # filter based on what to keep
  d_temp <- dplyr::filter(d_temp, gene_type %in% pos_to_keep)
  
  return(d_temp)
  
}

# get distance from 00
dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}

#-----------------#
# load in data ####
#-----------------#

# read in data
d <- readRDS('data/sequencing/invitro_pool.rds')

#-------------------#
# data wrangling ####
#-------------------#

# look at per sample how many new snps there are and how changes are lost and how many are shared
d_summary <- group_by(d, treat, rep, tech_rep) %>%
  mutate(., nrow = n()) %>%
  summarise(., lost = sum(is.na(af)),
            gain = sum(is.na(ancestral_prop)),
            same = unique(nrow) - lost - gain) %>%
  ungroup()

# pool over technical replicates, get rid of things that are not in both
d <- group_by(d, treat, rep, pos, ref, alt) %>%
  mutate(num = n()) %>%
  ungroup() %>%
  arrange(treat, rep, pos, ref, alt)

group_by(d, num) %>% tally()

# filter out ones that are only there once
d <- filter(d, num == 2)

# replace NAs for af (proportion of variant) with 0s
d <- mutate(d, af = replace_na(af, 0),
            ancestral_prop = replace_na(ancestral_prop, 0))

# calculate mean frequency and then filter when difference between technical replicates is > 0.3
d_sum <- group_by(d, treat, rep, pos, ref, alt, ancestral_prop) %>%
  summarise(., diff_af = max(af) - min(af),
            af = mean(af)) %>%
  ungroup()

hist(d_sum$diff_af)

# remove changes that are bigger than 0.3 difference between technical replicates
d_sum <- filter(d_sum, diff_af < 0.3)

# round frequency and ancestral prop to nearest 1/24
# calculate change to the ancestor
d2 <- mutate(d_sum,
             af = plyr::round_any(af, 1/24),
             ancestral_prop = plyr::round_any(ancestral_prop, 1/24),
             diff = af - ancestral_prop)

# filter out where both af and ancestral_prop are the same - i.e. no difference to the ancestor
d2 <- filter(d2, abs(diff) > 1/24)

# create a column for lost or novel variants, and when an existing variant has changed frequency
d2 <- mutate(d2, change = case_when(af == 0 ~ 'lost',
                                    ancestral_prop == 0 ~ 'gain',
                                    TRUE ~ 'existing'))

#---------------------------------------------------------------#
# filter out more snps to identify important genetic changes ####
#---------------------------------------------------------------#

# 1. filter to keep only genetic variants that are present in at least 3 replicates (of 6) of a treatment

d_filt <- filter(d2, af > 0) %>%
  mutate(genetic_change = group_indices(., pos, ref, alt)) %>%
  group_by(genetic_change, treat) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  group_by(pos, ref, alt) %>%
  mutate(genetic_change = cur_group_id()) %>%
  ungroup()

# get all positions of genetic changes
d_genetic_change <- select(d_filt, pos, ref,alt, genetic_change, ancestral_prop) %>% 
  distinct(.keep_all = TRUE)

# there are 960
unique(d_genetic_change$genetic_change) %>% length()

# need to create a dataframe of all instances of genetic changes - for when it isnt present becomes a 0
d_add <- group_by(d_filt, treat, rep) %>%
  do(data.frame(genetic_change = d_genetic_change$genetic_change)) %>%
  ungroup()

d_add <- left_join(d_add, d_genetic_change)

# add those into the dataframe, replace NAs with 0s
d_filt <- left_join(d_add, d_filt) %>%
  mutate(af = replace_na(af, 0)) %>%
  mutate(change = case_when(af == 0 ~ 'lost',
                            ancestral_prop == 0 ~ 'gain',
                            TRUE ~ 'existing'))

# 2. filter out only snps / indels that change significantly between groups
# perform paired wilcoxon test between control and mix_add to see if distribution of snps is different
d_changes <- d_filt %>%
  group_by(., pos, ref, alt, genetic_change) %>%
  nest() %>%
  mutate(., kruskal = purrr::map(data, ~kruskal.test(af ~ treat, .x)))

d_changes <- d_changes %>%
  mutate(tidy_model = map(kruskal, broom::tidy)) %>%
  unnest(tidy_model) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = 'BH'))

# look at the pvalue histogram
ggplot(d_changes, aes(x = padj)) +
  geom_histogram(boundary = 0)

unique(d_changes$p.value) %>% length()

# filter out SNPs wth a significant change
d_changes_to_keep <- filter(d_changes, padj < 0.05)

# filter for those genetic positions
d_filt <- filter(d_filt, genetic_change %in% d_changes_to_keep$genetic_change)

#-----------------------------#
# assign changes to a gene ####
#-----------------------------#

# read in gff file
# taken directly from https://cran.r-project.org/web/packages/rmonad/vignettes/gff-processing.html
gff <- readr::read_tsv('data/sequencing/PA14_gff.gff',
                col_names = c(
                  "seqid",
                  "source",
                  "type",
                  "start",
                  "stop",
                  "score",
                  "strand",
                  "phase",
                  "attr"
                ),
                na        = ".",
                comment   = "#",
                col_types = "ccciidcic"
)

# some hella long attributes section here
gff$attr[2]

# clean up the gff output
gff <- # remove first line (type == region) as it described the whole genome
  filter(gff, type != 'region' & type != 'rRNA' & type != 'CDS') %>%
  # some funky regex
  mutate(.,
         attr_id = gsub(pattern = "(.*ID=)(.*?)(;.*)", replacement = "\\2", attr),
         name = gsub(pattern = "(.*Name=)(.*?)(;.*)", replacement = "\\2", attr),
         # make instances where no name is present the type by id
         name = ifelse(grepl('=', name), paste(attr_id, type, sep = '_'), name),
         locus_tag = gsub(pattern = "(.*old_locus_tag=)(.*?)()", replacement = "\\2", attr)) %>%
  select(., seqid, type, start, stop, strand, attr_id, name, locus_tag)

# assign each SNP/indel to a gene
d_filt <- group_by(d_filt) %>%
  nest(pos) %>%
  mutate(., gene_info = purrr::map(data, ~get_gene_info(.x, gff))) %>%
  unnest_legacy(data) %>%
  unnest_legacy(gene_info)

# look at NAs for no gene info
filter(d_filt, gene_type == 'current_pos') %>%
  filter(is.na(gene_name)) %>%
  distinct(pos, .keep_all = TRUE)
# 29 changes where there is no gene info
# delete these

# merge with old locus tags
d_filt2 <- filter(d_filt, gene_type == 'current_pos') %>%
  filter(!is.na(gene_name)) %>%
  semi_join(., select(gff, gene_name = name, locus_tag))
  
d_filt2 <- mutate(d_filt2, af = replace_na(af, 0),
                   diff = af - ancestral_prop) %>%
  select(-diff_af) %>%
  mutate(., change = case_when(ancestral_prop > 0 & af == 0 ~ 'lost',
                               ancestral_prop == 0 ~ 'gain',
                               TRUE ~ 'existing'),
         change2 = ifelse(ancestral_prop == 0, 'gain', 'existing'))

# split into known and unknown genes
d_genes_unknown <- filter(d_filt2, stringr::str_detect(gene_name, 'PA14'))
d_genes_known <- filter(d_filt2, !stringr::str_detect(gene_name, 'PA14'))
unique(d_genes_unknown$genetic_change) %>% length()
unique(d_genes_known$genetic_change) %>% length()

ancestor <- select(d_genes_known, ancestral_prop, pos, gene_name, ref, alt) %>%
  distinct(., .keep_all = TRUE) %>%
  pivot_longer(., ancestral_prop, names_to = 'treat', values_to = 'af')

d_genes_known <- unite(d_genes_known, 'id', c(treat, rep), na.rm = TRUE, remove = FALSE) %>%
  mutate(., id2 = as.numeric(as.factor(id)),
         change2 = ifelse(is.na(change2), 'ancestor', change2))

d_genes_known <- unite(d_genes_known, 'id3', c(gene_name, pos), na.rm = TRUE, remove = FALSE)

# save out as a table
write_csv(d_genes_known, 'data/sequencing/invitro_snps_indels.csv')

#------------------#
# do clustering ####
#------------------#

# cluster just on the 50 significant changes where gene functions are known
head(d_genes_known)

# setup data for clustering
d_ancest = mutate(ancestor, treat = 'ancestor',
                  rep = '')

d_clust <- select(d_genes_known, treat, rep, pos, ref, alt, af, gene_name) %>%
  bind_rows(., d_ancest) %>%
  unite(., 'id', c(treat, rep), sep ='.') %>%
  spread(., id, af) %>%
  arrange(pos) %>%
  group_by(gene_name) %>%
  mutate(n = 1:n()) %>%
  ungroup() %>%
  mutate(gene_name = ifelse(n > 1, paste(gene_name, n, sep = '_'), gene_name)) %>% 
  select(., -c(ref, alt, pos, n))

d_clust <- tibble::column_to_rownames(d_clust, 'gene_name')

# transpose the dataframe
d_clust <- t(d_clust)

# do nmds
nmds <- vegan::metaMDS(d_clust, distance = 'euclidean')

stressplot(nmds)
plot(nmds)

# vars
d_vars <- select(d_genes_known, gene_name, change2, change, pos, alt, ref, alt, ancestral_prop) %>%
  distinct() %>%
  group_by(gene_name) %>%
  mutate(n = 1:n()) %>%
  ungroup() %>%
  mutate(gene_name = ifelse(n > 1, paste(gene_name, n, sep = '_'), gene_name))

# get data from nmds
d_nmds <- fortify(nmds) %>%
  janitor::clean_names() %>%
  mutate_if(., is.factor, as.character)

# wrangle sites
d_sample <- filter(d_nmds, score== 'sites') %>%
  mutate(label = gsub('mix.add', 'mix_add', label)) %>%
  separate(., label, c('treat', 'rep'), sep = '\\.')

# wrangle species
d_gene <- filter(d_nmds, score == 'species') %>%
  rename(., gene_name = label) %>%
  mutate(dist = dist_from_00(nmds1, nmds2)) %>%
  merge(., d_vars, by = 'gene_name') %>%
  mutate(gene_name = gsub('_.*', '', gene_name))

d_gene2 <- select(d_gene, gene_name, nmds1, nmds2, dist, change2) %>%
  distinct()

#------------------#
# make Figure 5 ####
#------------------#

p_nmds <- ggplot() +
  geom_segment(aes(x = 0, y = 0, yend = nmds2, xend = nmds1, group = score, col = change2), d_gene, arrow = arrow(length = unit(0.01, "npc"))) +
  geom_label(aes(nmds1, nmds2, label = gene_name, hjust = -sign(nmds1), vjust = -0.3*sign(nmds2), col = change2), filter(d_gene2, nmds1 < 0 & nmds2 > 0), show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(nmds1, nmds2, label = gene_name, col = change2), 
                            filter(d_gene2, nmds1 < 0 & nmds2 < 0), 
                            segment.colour = "grey", force_pull = 1, nudge_y = -0.2,
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.5, "lines"),
                            max.overlaps = 10, show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(nmds1, nmds2, label = gene_name, col = change2), 
                            filter(d_gene2, nmds1 > 0 & nmds2 < 0), 
                            segment.colour = "grey", force_pull = 0.00001, nudge_y = -0.7, xlim = c(0, 6), ylim = c(1,-3.5),
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.5, "lines"),
                            max.overlaps = 10, show.legend = FALSE) +
  ggrepel::geom_label_repel(aes(nmds1, nmds2, label = gene_name, col = change2), 
                            filter(d_gene2, nmds1 > 0 & nmds2 > 0), 
                            segment.colour = "grey", force_pull = 0.00001, nudge_y = 0.5, nudge_x = 1, xlim = c(0, 6), ylim = c(0,4),
                            box.padding = unit(0.35, "lines"),
                            point.padding = unit(0.5, "lines"),
                            max.overlaps = 10, show.legend = FALSE) +
  ggforce::geom_mark_hull(aes(nmds1, nmds2, group = treat), filter(d_sample, treat != 'ancestor')) +
  geom_point(aes(nmds1, nmds2, shape = treat), size = 5, fill = 'white', data = filter(d_sample, treat != 'ancestor')) +
  geom_point(aes(nmds1, nmds2, shape = treat), size = 6, fill = '#7090A0', data = filter(d_sample, treat == 'ancestor')) +
  scale_x_continuous(expand = c(.25, .25)) +
  scale_y_continuous(expand = c(.25, .25))  +
  labs(title = '(a) NMDS of euclidean distances between genome changes of known function',
       x = 'nmds 1',
       y = 'nmds 2',
       col = "Change type",
       shape = "Treatment") +
  theme_bw(base_size = 14) +
  scale_color_manual(name = NULL, values = rev(ichooseyou('vileplume', spread = 2)), labels = c(expression(atop("Existing","variant")), expression(atop(italic(de~novo),"variant"))), guide = 'legend') +
  scale_shape_manual(name = NULL, values = c(21:24), labels = c("Ancestor", "Control", "Phage added\nonce", "Phage added\nrepeatedly")) +
  theme(legend.position = "bottom", legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 1.2, linetype = 1)),
         shape = guide_legend(override.aes = list(fill = c('#7090A0', 'white', 'white', 'white'))))
  
p_nmds

#----------------------------------#
# do analyses on genomic traits ####
#----------------------------------#

# 1. calculate distance from ancestor
d_snp_summary <- d_genes_known %>%
  group_by(treat, rep) %>%
  summarise(dist_ancest = sum(abs(diff)), .groups = 'drop')
d_snp_summary2 <- d_genes_known %>%
  filter(change2 != 'existing') %>% 
  filter(af > 0) %>%
  group_by(treat, rep) %>%
  tally() %>%
  select(treat, rep, new_snps_indels = n) %>%
  ungroup()
d_snp_summary <- full_join(d_snp_summary, d_snp_summary2) %>%
  mutate(across(c(new_snps_indels, dist_ancest), replace_na, 0))

p1 <- ggplot(d_snp_summary, aes(treat, dist_ancest)) +
  geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(aes(shape = treat), fill = 'white', size = 5, position = position_jitter(width = 0.1)) +
  ylab('genotypic distance from ancestor') +
  ggtitle('(b)') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_shape_manual(values = c(22:24)) +
  scale_x_discrete(labels = c('Control', 'Phage\nadded\nonce', 'Phage\nadded\nrepeatedly')) +
  scale_y_continuous(limits = c(0, 30), n.breaks = 4)

p2 <- ggplot(d_snp_summary, aes(treat, new_snps_indels)) +
  geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(aes(shape = treat), fill = 'white', size = 5, position = position_jitter(width = 0.1, height = 0)) +
  ylab('number of SNPs / indels') +
  ggtitle('(c)') +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12)) +
  scale_shape_manual(values = c(22:24)) +
  scale_x_discrete(labels = c('Control', 'Phage\nadded\nonce', 'Phage\nadded\nrepeatedly')) +
  scale_y_continuous(limits = c(0, 25), n.breaks = 4)

# alpha diversity - heterozygosity

# calculate Hardy Weinberg equilibrium for each SNP
alpha_div <- mutate(d_genes_known, p = af,
                    q = 1-p,
                    div = 1 - p^2 - q^2) %>%
  group_by(., treat, rep) %>%
  summarise(., diversity = sum(div), .groups = 'drop')

# plot alpha_div ####
p3 <- ggplot(alpha_div, aes(treat, diversity)) +
  geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(aes(shape = treat), fill = 'white', position = position_jitter(width = 0.1), size = 5) +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  ggtitle('(d)') +
  xlab('') +
  ylab('alpha diversity') +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 12)) +
  scale_shape_manual(values = c(22:24)) +
  scale_x_discrete(labels = c('Control', 'Phage\nadded\nonce', 'Phage\nadded\nrepeatedly')) +
  scale_y_continuous(limits = c(0, 7.5), n.breaks = 4)

patchwork <- p_nmds + (p1 + p2 + p3)  +
  plot_layout(ncol = 1, heights = c(0.7, 0.3))

patchwork

ggsave('plots/Figure_5.pdf', patchwork, width = 11, height = 10)
ggsave('plots/Figure_5.png', patchwork, width = 11, height = 10)

#---extra analyses not used in the manuscript---#

# data analysis comparing SNPs / distance from ancestor among groups
gen_m <- lm(new_snps_indels ~ treat, data = d_snp_summary)
gen_m0 <- lm(new_snps_indels ~ 1, data = d_snp_summary)
anova(gen_m, gen_m0, test = "F") #significant effect of treatment
emmeans::emmeans(gen_m, pairwise ~ treat) #all comparisons significant

dist1 <- lm(dist_ancest ~ treat, data = d_snp_summary)
dist0 <- lm(dist_ancest ~ 1, data = d_snp_summary)
anova(dist1, dist0, test = "F")
emmeans::emmeans(dist1, pairwise ~ treat) #all comparisons significant

alph_m1 <- lm(diversity ~ treat, data = alpha_div)
alph_m0 <- lm(diversity ~ 1, data = alpha_div)
anova(alph_m1, alph_m0, test = "F")
emmeans::emmeans(alph_m1, pairwise ~ treat) #no difference between control and mix but sig differences between control and mix_add and mix and mix_add

# distance matrices
euclid_matrix <- dist(d_clust)

# make some character vectors factors
d_vars <- mutate_at(d_sample, c('treat'), as.factor)
row.names(d_vars) <- paste(d_vars$treat, d_vars$rep, sep ='.')

# adonis
mod_euclid <- adonis(euclid_matrix ~ treat, data = d_vars, permutations = 9999)

# Looking at variance across groups
# can only do one factor so lets use id for every combination
mod_dispers <- betadisper(euclid_matrix, d_vars$treat)

# plot of model
plot(mod_dispers)
boxplot(mod_dispers)

# anova
anova(mod_dispers)

# Permutation test for F
pmod <- permutest(mod_dispers, pairwise = TRUE)

# Tukey's Honest Significant Differences
T_HSD <- TukeyHSD(mod_dispers)
