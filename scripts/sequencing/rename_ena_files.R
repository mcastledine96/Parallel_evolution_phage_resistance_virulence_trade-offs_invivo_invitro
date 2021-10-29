#---------------------------------------------------------#
# script for renaming ena files after downloading them ####
#---------------------------------------------------------#

# only run this after downloading files using ena_download.sh

# load in packages
library(tidyverse)
library(progress)

# load in data that encodes information on files
d <- read.csv('data/info_for_renaming_ena_files.csv', stringsAsFactors = FALSE)

# list files where you downloaded them from ena
files <- list.files('~/PRJEB47945', recursive = TRUE, full.names = TRUE)

# make new directory for clone and pool-seq files
dir.create('~/PRJEB47945/clone')
dir.create('~/PRJEB47945/pool')

# set up progress bar
pb <- progress_bar$new(total = length(files))

# rename all files into their respective folders (clone or pool-seq)
# rename them as their sample name
for(i in 1:length(files)){
  pb$tick()
  file = files[i]
  file2 = basename(file)
  id <- dirname(file) %>% basename()
  fwd_rev <- gsub('.*_', '', file2) %>% parse_number() %>% paste('R', ., sep = '')
  
  id2 <- filter(d, run_accession == id)
  
  if(str_detect(id2$sample_name, 'clone')){
    file.copy(file, paste('~/PRJEB47945/clone/', id2$sample_name, '.', fwd_rev, '.fastq.gz', sep = ''))
  }
  if(str_detect(id2$sample_name, 'pool')){
    file.copy(file, paste('~/PRJEB47945/pool/', id2$sample_name, '.', fwd_rev, '.fastq.gz', sep = ''))
  }
  file.remove(file)
}
