#----------------------------------------------------#
# clean up pool seq in vitro sequencing vcf files ####
#----------------------------------------------------#

# summary of what this script does

# 1. reads in ancestor vcf file and calculates proportion of each variant
# 2. reads in in vitro vcf files and calculates proportion of each variant

# load packages
library(vcfR)
library(tidyverse)

# function to clean up freebayes vcf files
# filter depth to be no more than 2 s.d. from the mean
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

#---------------------------------------#
# first process the ancestor samples ####
#---------------------------------------#

# list files
files_all <- list.files('data/sequencing/vcf/pool', full.names = TRUE, pattern = '.vcf')
files_ancest <- files_all[grepl('ancestors', files_all)]

# map ancestral reads
d_ancest <- files_ancest %>%
  map_df(., tidy_freebayes)

# split up file
d_ancest <- separate(d_ancest, file, c('treat', 'blah', 'tech_rep', 'mapping_tool', 'ref_genome', 'blah1', 'blah2'), sep = '_')

# remove gt_ad as it is not needed
d_ancest <- select(d_ancest, -c(gt_ad, blah, blah1, blah2))

# filter samples that have multiple variants
d_ancest_mult <- filter(d_ancest, numalt != '1')

# split these up into multiple rows
cols_to_stack <- stringr::str_detect(d_ancest_mult[1,], ',') %>%
  names(d_ancest_mult)[.] %>%
  .[!is.na(.)]
d_ancest_mult <- separate_rows(d_ancest_mult, cols_to_stack, sep = ',')

# rebind with the other samples
d_ancest <- bind_rows(filter(d_ancest, numalt == '1'), d_ancest_mult)

# make a load of columns numeric
d_ancest <- mutate_at(d_ancest, vars(ac,af,ao,pao,qa,pqa,saf,sap,sar,ab,abp,run,rpp,rpl,rpr,epp,dpra,mqm), as.numeric)

# filter on quality score
d_ancest_filt <- filter(d_ancest, qual > 30)

# conservatively calling SNPs on the ancestors would be when they are there in either one of the lanes
d_ancest_snps <- group_by(d_ancest_filt, pos) %>%
  mutate(., n = n()) %>%
  ungroup() %>%
  mutate(., number = ifelse(n == 1, 'single_lane', 'both')) 

# look at how many are present in one of the lanes only
group_by(d_ancest_snps, number) %>% tally()

# calculate average proportion of snp across both samples
d_ancest_snps <- filter(d_ancest_snps, number == 'both') %>%
  group_by(., pos, chrom, ref, alt, treat) %>%
  summarise_at(., vars(dp, af, ac, ao), .funs = list(function(x) mean(x, na.rm = TRUE), sd)) %>%
  ungroup() %>%
  rename_at(., vars(ends_with('_fn1')), function(x) gsub('_fn1', '_mean', x)) %>%
  rename_at(., vars(ends_with('_fn2')), function(x) gsub('_fn2', '_sd', x))

head(d_ancest_snps)

saveRDS(select(d_ancest_snps, pos, ref, alt, af_mean), 'data/sequencing/ancestors.rds')

#-----------------------------------#
# read in invitro pooled samples ####
#-----------------------------------#

files_invitro <- files_all[!files_all %in% files_ancest]

# stack all data together
d_invitro <- files_invitro %>%
  map_df(., tidy_freebayes)

# create treatment columns
d_invitro <- mutate(d_invitro, file = gsub('mix_add', 'mixadd', file)) %>%
  separate(., file, c('treat', 'rep', 'blah', 'tech_rep', 'mapping_tool', 'ref_genome', 'blah1', 'blah2'), sep = '_') %>%
  select(., -starts_with('blah'))

# remove gt_ad as it is not needed
d_invitro <- select(d_invitro, -gt_ad)

# filter samples that have multiple variants
d_invitro_mult <- filter(d_invitro, numalt != '1')

# split these up into multiple rows
cols_to_stack <- stringr::str_detect(d_invitro_mult[1,], ',') %>%
  names(d_invitro_mult)[.] %>%
  .[!is.na(.)]
d_invitro_mult <- separate_rows(d_invitro_mult, cols_to_stack, sep = ',')

# rebind with the other samples
d_invitro <- bind_rows(filter(d_invitro, numalt == '1'), d_invitro_mult)

# make a load of columns numeric
d_invitro <- mutate_at(d_invitro, vars(ac,af,ao,pao,qa,pqa,saf,sap,sar,ab,abp,run,rpp,rpl,rpr,epp,dpra,mqm), as.numeric)

# filter on quality score
d_invitro_filt <- filter(d_invitro, qual > 30)

# want to merge with the ancestor, keep everything
# number of positions / snps which are not in the ancestor, but in the sample
# number of positions / snps which are in the ancestor, but not in the sample
# number which are in both

# bindable dataframe for ancestral samples
d_ancest_bind <- select(d_ancest_snps, pos, ref, alt, ancestral_prop = af_mean)

# bind with sample
d_invitro_merge <- select(d_invitro_filt, pos, ref, alt, af, treat, rep, tech_rep) %>%
  group_by(., treat, rep, tech_rep) %>%
  nest_legacy() %>%
  mutate(., new_d = purrr::map(data, ~merge(.x, d_ancest_bind, by = c('pos', 'ref', 'alt'), all = TRUE))) %>%
  select(-data) %>%
  unnest_legacy(new_d)

# save data out
saveRDS(d_invitro_merge, 'data/sequencing/invitro_pool.rds')
