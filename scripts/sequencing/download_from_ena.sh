#------------------------------------------------#
# example code for downloading files from ena ####
#------------------------------------------------#

# first download latest version of enaBrowserTools
wget https://github.com/enasequence/enaBrowserTools/archive/refs/tags/v0.0.3.zip
unzip v0.0.3.zip

# follow the instructions on enaBrowserTools GitHub
# https://github.com/enasequence/enaBrowserTools

# download all fastq files
enaGroupGet -f fastq -g read PRJEB47945

# Now move on to rename_ena_files.R