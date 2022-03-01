# Parallel evolution of _Pseudomonas aeruginosa_ phage resistance and virulence loss in response to phage treatment _in vivo_ and _in vitro_

Coding and data repository for research article: "Parallel evolution of _Pseudomonas aeruginosa_ phage resistance and virulence loss in response to phage treatment _in vivo_ and _in vitro_". 

- Accepted in eLife: [DOI: 10.7554/eLife.73679](https://elifesciences.org/articles/73679) 
- Preprint available on biorxiv: [DOI: 10.1101/2021.09.06.459069](https://www.biorxiv.org/content/10.1101/2021.09.06.459069v2)
- Sequencing data available on ENA: [Accession code: PRJEB47945](https://www.ebi.ac.uk/ena/browser/view/PRJEB47945)

## Outline

This repository contains the final datasets and R scripts to run all the analyses and create the figures present in the main text of the above-mentioned paper. 

## Feedback

Please report any problems or bugs in the code in the Issues tab of this GitHub repository. Alternatively, please email mcastledine96@gmail.com

## Licensing

This code is licensed under GPL-3.

## Running the scripts and analyses

- The project can be cloned or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. R projects automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in **scripts**.

### Phenotype scripts and analyses

- these can be found in `scripts/phenotype`
- 1. **time_series_invitro.R** and **time_series_invivo.R** are scripts used for analysis of bacteria - phage coevolution and levels of bacterial resistance to phage. To generate Figure 1, run both scripts following in-script prompts within the same RStudio session. Data files used for this analysis include **invitro_timeseries.csv**, **invivo_timeseries.csv** and **anc_resistance.csv** respectively.
- 2. **virulence_models.R** includes analysis of the *Galleria* virulence assays for both the in vivo and in vitro clones. Data for this file is found in **data/phenotype/virulence_data.csv**. This script recreates **Figure 2** of the manuscript.
- 3. **growth_analysis.R** includes analysis of in vivo and in vitro data, growth rates of clonal populations for phenotypic analysis. Data for this file is found in **invitro_growth_rates.csv** and **invivo_growth_rates.csv**
- 4. **biofilm_analysis.R** includes analysis of in vivo and in vitro data, biofilm data of clonal populations for phenotypic analysis. Data for this file is found in **biofilm_data.csv**

### Sequencing scripts and analysis

- these can be found within `scripts/sequencing`
- 1. **download_from_ena.sh** is a bash script to facilitate the download of the raw whole genome sequencing files (fastq.gz) from the European Nucleotide Archive. They are publicly available under the accession number [PRJEB47945](https://www.ebi.ac.uk/ena/browser/view/PRJEB47945).
- 2. **invitro_pool_mapping_pipeline.sh** and **invivo_clone_mapping_pipeline.sh** are bash scripts to run the bioinformatic analysis of the raw sequencing files. They filter the raw reads, map to the reference genome, and call variants using freebayes. The resulting .vcf files are saved and available in `data/sequencing/vcf`.
- 3. **invitro_process_vcf.R** processes the pool seq vcf files and produces **ancestors.rds** and **invitro_pool.rds** present in `data/sequencing` that are used in subsequent analyses.
- 4. **invitro_analysis_genomics.R** runs extra processing on the pool seq data and runs the genomics analysis on pool seq data. It produces **Figure 5** of the manuscript.
- 5. **gene_change_invivo_invitro.R** processes the invivo clone sequencing data and produces **Figure 4** of the manuscript.

# Data files

- Currently, the meaning of column names and variables are described in the relevant R script. Please see above to find where this information can be found. If anything is missing or unclear, please feel free to get in contact to the above email address and we will aim to respond promptly. 
