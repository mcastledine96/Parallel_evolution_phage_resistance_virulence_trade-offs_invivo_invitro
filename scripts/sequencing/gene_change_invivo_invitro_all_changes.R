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
invitro_snps <- bind_rows(read.csv('data/sequencing/invitro_snps_indels.csv'),
                          read.csv('data/sequencing/invitro_snps_indels_unknown_genes.csv')) %>%
  filter(!is.na(ref)) %>%
  dplyr::select(-n) %>%
  select(., gene_name, ref, alt, id3, pos, change2) %>%
  distinct() %>%
  rename(ref_invitro = ref, alt_invitro = alt)

# load in invitro data
d_invitro <- bind_rows(read.csv('data/sequencing/invitro_snps_indels.csv'),
                       read.csv('data/sequencing/invitro_snps_indels_unknown_genes.csv')) %>%
  dplyr::select(-n) %>%
  select(., af, gene_name, id3, change2, id2, pos) %>%
  mutate(id4 = 'invitro populations',
         id2 = 1+id2)

# add ancestors to dataset
ancestor <- bind_rows(read.csv('data/sequencing/invitro_snps_indels.csv'),
                      read.csv('data/sequencing/invitro_snps_indels_unknown_genes.csv')) %>%
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

# add column for known vs unknown gene
d3 <- mutate(d3, gene_status = ifelse(str_detect(id3, 'PA14'), 'unknown', 'known'))
d_invitro <- mutate(d_invitro, gene_status = ifelse(str_detect(id3, 'PA14'), 'unknown', 'known'))

# have a look at how many were of the invivo snps are present invitro
select(d_invitro, id3, gene_status) %>%
  distinct() %>%
  group_by(gene_status) %>%
  tally()
select(d3, id3, gene_status) %>%
  distinct() %>%
  group_by(gene_status) %>%
  tally()

# unknown vs known gene % present
31/50
153/234

# have a look at how many of the invivo snps are present invitro
select(d_invitro, gene_name, gene_status) %>%
  distinct() %>%
  group_by(gene_status) %>%
  tally()
select(d3, gene_name, gene_status) %>%
  distinct() %>%
  group_by(gene_status) %>%
  tally()

# unknown vs known gene % present
16/26
85/99

group_by(d3, gene_status) %>%
  summarise(mean_af = mean(af))

select(d_invitro, gene_status, change2, id3) %>%
  distinct() %>%
  filter(!is.na(change2)) %>%
  group_by(gene_status, change2) %>%
  tally()

select(d3, gene_status, change2, id3) %>%
  distinct() %>%
  group_by(gene_status, change2) %>%
  tally()

24/(24+7)
71/(71+82)

# look at correlation between average proportion in treatments in vitro and those seen in vivo
d_invitro_summary <- bind_rows(read.csv('data/sequencing/invitro_snps_indels.csv'),
                                            read.csv('data/sequencing/invitro_snps_indels_unknown_genes.csv')) %>%
  group_by(treat, id3) %>%
  summarise(mean_af = mean(af), .groups = 'drop') %>%
  mutate(., gene_status = ifelse(str_detect(id3, 'PA14'), 'unknown', 'known'))

d_invitro_summary2 <- bind_rows(read.csv('data/sequencing/invitro_snps_indels.csv'),
                               read.csv('data/sequencing/invitro_snps_indels_unknown_genes.csv')) %>%
  group_by(id3, change2, treat) %>%
  summarise(mean_af = mean(af), .groups = 'drop') %>%
  mutate(., gene_status = ifelse(str_detect(id3, 'PA14'), 'unknown', 'known'))

d_correlation <- d_invitro_summary2 %>%
  left_join(., select(d3, id3, invivo_af = af)) %>%
  mutate(invivo_af = replace_na(invivo_af, 0)) %>%
  filter(mean_af > 0)

ggplot(d_correlation, aes(mean_af, invivo_af)) +
  geom_point() +
  facet_wrap(~treat)

unique(d_invitro$id3) %>% length()
284*3
nrow(d_invitro_summary)


