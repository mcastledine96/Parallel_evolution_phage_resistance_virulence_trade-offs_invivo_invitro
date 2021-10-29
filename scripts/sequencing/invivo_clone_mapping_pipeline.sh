#------------------------------------------------------------------------------#
# processing and variant calling of clonal data from in vivo experiments ####
#------------------------------------------------------------------------------#

# prerequisites
# Work through ena_download.sh
# Ran rename_ena_files.R
# Once there is a PRJEB47945/pool with the raw files in - you can proceed
# download reference genome and place it in ~/PRJEB47945/ref
# used https://www.ncbi.nlm.nih.gov/assembly/GCF_000014625.1
# Pseudomonas aeruginosa UCBPP-PA14

mkdir -p ~/PRJEB47945/ref
# save as PA14_reference_genome.fna

#------------------#
# software used ####
#------------------#

# freebayes
# fastqc
# bcftools
# samtools
# bwa
# samblaster
# vcflib
# trim-galore

#--------------------------#
# example code workflow ####
#--------------------------#

# first install mamba as its faster for installing than conda
conda install -c conda-forge mamba

# create an environment
mamba create -n sequencing_pipeline 

# open the environment
conda activate sequencing_pipeline

# install packages
mamba install -c bioconda freebayes fastqc bcftools samtools bwa samblaster vcflib
mamba install -c bioconda biopet-vcffilter
mamba install -c bioconda trim-galore

# set working directory at the base folder
wd=~/PRJEB47945/clone

# set other directories
raw_files=$wd
trimmed_files=$wd/trimmed_files
fastqc_output=$wd/fastqc_output

# check
cd $wd

# run trim galore on each file ####
mkdir -p $trimmed_files

# trim illumina reads
for fwd_file in $raw_files/*R1.fastq.gz; do
        
        rev_file=${fwd_file%.R1.fastq.gz}.R2.fastq.gz
        
        trim_galore --paired --quality 10 --output_dir $trimmed_files $fwd_file $rev_file

done

# run fastqc on all files
mkdir -p $fastqc_output
fastqc -o $fastqc_output $trimmed_files/*.fastq.gz

#-------------------------#
# map to reference genome #
#-------------------------#

# make directories for bam files
bam_files=$wd/processed/bams
mkdir bam_files

ref_pa14=~/PRJEB47945/ref/PA14_reference_genome.fna

# prepare the reference for mapping
samtools faidx $ref_pa14
bwa index $ref_pa14

# map reads to reference genome
for file in $trimmed_files/*R1_val_1.fastq.gz

    do
        echo $file
        
        # assign fwd and rev files
        file_fwd=$file
        file_rev=${file%_R1_val_1.fastq.gz}_R2_val_2.fastq.gz

        #stub
        stub=$(basename ${file%.R1_val_1.fastq.gz})

        # do the mapping for fwd and reverse pairs
        # use samblaster to mark duplicates
        
        bwa mem -t 12 $ref_pa14 $file_fwd $file_rev | samblaster | samtools view -S -b -o $wd/processed/bams/${stub}_bwa.bam

done

# make directories for mapped and unmapped files
mkdir -p $bam_files/mapped
mkdir -p $bam_files/unmapped
mkdir -p $bam_files/stats

for bam_file in $wd/processed/bams/*.bam

    do
        echo $bam_file
        stub=${bam_file%.bam}
        stub2=$(basename ${bam_file%.bam})

        # 1. split bam into mapped and unmapped reads
        samtools view -b -F4 $bam_file > $bam_files/mapped/${stub2}_refpa14_mapped.bam
        samtools view -f4 $bam_file > $bam_files/unmapped/${stub2}_refpa14_unmapped.bam

        # 2. sort mapped file by position in genome and not by order of mapped sequence
        samtools sort -o $bam_files/mapped/${stub2}_refpa14_sorted.bam $bam_files/mapped/${stub2}_refpa14_mapped.bam

        # 3. index the sorted bam file
        samtools index $bam_files/mapped/${stub2}_refpa14_sorted.bam

        # 3. remove intermediate (unsorted mapped bam file)
        rm $bam_files/mapped/${stub2}_refpa14_mapped.bam

# 4. extract stats of file
        # sorted and mapped bam file
        samtools flagstat $bam_files/mapped/${stub2}_sorted.bam > $wd/processed/bams/stats/${stub2}_mapped_stats.txt
        
        # raw bam file
        samtools flagstat $bam_file > $wd/processed/bams/stats/${stub2}_raw_stats.txt

done

#--------------------------------------------------#
# run freebayes to call variants for each bam file #
#--------------------------------------------------#

# make directory for output
mkdir -p $wd/processed/vcf_output
mkdir -p $wd/processed/vcf_output/freebayes

for bam_file in $wd/processed/bams/mapped/*sorted.bam

    do
        echo $bam_file
        stub=${bam_file%.bam}
        stub2=$(basename ${bam_file%.bam})
        
        # skip if file already exists
        if [ -e $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf ]
        then
                continue
        fi

        # call SNPs using FreeBayes, set ploidy to 1
        freebayes -f $ref_pa14  -p 1 $bam_file > $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf
        # keep only snps with a QUAL score > 20, based on recommendation of Erik Garrison
        #https://www.biostars.org/p/71299/
        vcffilter -f "QUAL > 20" $wd/processed/vcf_output/freebayes/${stub2}_freebayes.vcf > $wd/processed/vcf_output/freebayes/${stub2}_freebayes_filter.vcf
done
