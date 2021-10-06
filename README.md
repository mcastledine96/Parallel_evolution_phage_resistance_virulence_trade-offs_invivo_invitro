# Parallel_evolution_phage_resistance_virulence_trade-offs_invivo_invitro
Coding and data repository for manuscript: "Parallel evolution of phage resistance - virulence trade - offs during in vitro and nasal Pseudomonas aeruginosa phage treatment". Currently available on biorxiv: https://doi.org/10.1101/2021.09.06.459069


# Outline
This repository contains the final datasets, R scripts containing the analyses, and figures of the above-mentioned paper. It can recreate the analyses and figures in the main text.

# Feedback
Please report any problems or bugs in the code in the Issues tab of the GitHub repository. Alternatively, please email mcastledine96@gmail.com
Licensing
This code is licensed under GPL-3.

# Running the scripts and analyses
- The project can be cloned or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. R projects automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in scripts/.
- time_series_invitro.R and time_series_invivo.R are scripts used for analysis of bacteria - phage coevolution and levels of bacterial resistance to phage. To generate Figure 1, run both scripts following in-script prompts within the same RStudio session. Datafiles used for this analysis include In_vitro_timeseries.csv and In_vivo_timeseries.csv respectively.
- virulence_models.R includes analysis of in vivo and in vitro data, virulence of clonal populations for phenotypic analysis. Data for this file is found in virulence_data.csv 
- growth_analysis.R includes analysis of in vivo and in vitro data, growth rates of clonal populations for phenotypic analysis. Data for this file is found in invitro_growth_rates.csv and invivo_growth_rates.csv
-  biofilm_analysis.R includes analysis of in vivo and in vitro data, biofilm data of clonal populations for phenotypic analysis. Data for this file is found in biofilm_data.csv

# Data files
- Currenly, the meaning of column names and variables are described in the relevant R script. Please see above to find where this information can be found. If anything is missing or unclear, please feel free to get in contact to the above email address and we will aim to respond promptly. 
