# Qiime2_Mariana_forearc_16S_Amplicons
Updated 9.17.25

How to process 16S rRNA amplicon sequences using Qiime2. These samples were very low-biomass, from pH 12.5 fluid. Contaminant removal follows methods in Sheik et al., 2018.

45 mL of serpentinized fluid incubated with <sup>13</sup>CH<sub>4</sub> gas to track methanotrophy was extracted using phenol:chloroform methods described in Zhou et al., 1996 and Crump et al., 2003 and modified by Dr. Caroline Fortunato and Dr. Amy Smith. DNA was amplified using two EMP primer sets (below). This pipeline includes code for both of them. This code in modified from Dr. Sarah Hu's github: https://github.com/shu251/qiime2-2024/tree/main?tab=readme-ov-file.

## Primer sets: 
- 515F (Parada et al., 2016), 806R (Apprill et al., 2015)
- 515F (Parada et al., 2016), 926R (Quince et al., 2011)

## Adapters: 
- Forward: CTGTCTCTTATACACATCT
- Reverse: CTGTCTCTTATACACATCT

## Step 1: Install Qiime2 (v2024.10.1) on Poseidon using a singularity containiner: 
```
#enter scavenger 
srun -p scavenger --time=04:00:00 --ntasks-per-node 1 --mem=5gb --pty bash

#install using a singularity container
module load singularity/3.7
singularity pull \
  docker://quay.io/qiime2/amplicon:2024.2

#test it to see if it works
singularity exec amplicon_2024.2.sif qiime --help
#docker image file should be around 2 gb.
```

## Step 2: Git clone Dr. Sarah Hu's Qiime2 repo to get the code for the manifest file
```
#navigate to the directory where your sequences are stored. 
#git clone
git clone git@github.com:shu251/qiime2-2024.git

#If it won’t let you git clone, download the GitHub repo to your local computer and scp it into the directory where your sequences are stored. 
scp -r /Users/sabrinaelkassas/Downloads/qiime2-2024-main selkassas@poseidon.whoi.edu:/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences 
```


## Step 3: Modify the create-manifest R script from Dr. Sarah Hu's Github (printed below with my edits)
### 515F/806R
```
#my file paths: /vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R
#my file name structure (showing two, just in case):
1-dil1to10-J21475-10C-515F-806R_S13_L001_R1_001.fastq.gz
1-dil1to10-J21475-10C-515F-806R_S13_L001_R2_001.fastq.g

GNU nano 2.3.1                                       File: create-manifest.R                                                                                      

# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import

# Change working directory to sequences
setwd("/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R/")

# Required libraries
library(tidyverse)

# Import all files in current directory
paths <- as.data.frame(list.files(pattern = ".fastq.gz", full.names = TRUE))

colnames(paths)[1]<-"FASTQ" #this makes a df with the first column name “FASTQ” that will list all of your fastq file names

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files #lists fastq file paths

# Extract sample ID
paths_run <- paths %>%
## Use for sampleid_L001_Rx_001.fastq.gz
        mutate(SAMPLEID = str_extract(FASTQ, "^.+(?=_L001_R[12]_001.fastq.gz)")) #This command creates a new column called SAMPLEID by extracting and replacing part of the filename stored in the FASTQ column.

## Use for _R1.fastq.gz
#	mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_R((\\d)).fastq.gz","\\1")) %>%
#	separate(SAMPLEID, c("SAMPLEID", "else"), sep = "_L001_")
#
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID) #This is a base R function that searches for a pattern in a string and replaces all occurrences with something else.
paths_run$FASTQ <- gsub("./","", paths_run$FASTQ)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$FASTQ, sep="/") #make a new column with the full path

forward <- paths_run %>%
        filter(grepl("R1_001.fastq.gz", FASTQ)) %>% 
        select(SAMPLEID, `forward-absolute-filepath` = FULL_PATH)

# Write manifest file
manifest <- paths_run %>%
        filter(grepl("R2_001.fastq.gz", FASTQ)) %>%
        select(SAMPLEID, `reverse-absolute-filepath` = FULL_PATH) %>%
        right_join(forward) %>%
        select('sample-id' = SAMPLEID, `forward-absolute-filepath`, `reverse-absolute-filepath`) #join forward and reverse file paths

# Write output as a manifest file
write.table(manifest, file = "/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R/manifest_515F_806R.tsv", quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t") #write manifest file
```

### 515F/926R
```
#my file paths: /vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R
#my file name structure (showing two, just in case):
1-dil1to10-J21475-10C-515F-806R_S13_L001_R1_001.fastq.gz
1-dil1to10-J21475-10C-515F-806R_S13_L001_R2_001.fastq.g

GNU nano 2.3.1                                       File: create-manifest.R                                                                                      

# Read in SRA sample IDs and list of raw fastq files to generate a manifest file ready for qiime2 import

# Change working directory to sequences
setwd("/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R/")

# Required libraries
library(tidyverse)

# Import all files in current directory
paths <- as.data.frame(list.files(pattern = ".fastq.gz", full.names = TRUE))

colnames(paths)[1]<-"FASTQ" #this makes a df with the first column name “FASTQ” that will list all of your fastq file names

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files #lists fastq file paths

# Extract sample ID
paths_run <- paths %>%
## Use for sampleid_L001_Rx_001.fastq.gz
        mutate(SAMPLEID = str_extract(FASTQ, "^.+(?=_L001_R[12]_001.fastq.gz)")) #This command creates a new column called SAMPLEID by extracting and replacing part of the filename stored in the FASTQ column.

## Use for _R1.fastq.gz
#	mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_R((\\d)).fastq.gz","\\1")) %>%
#	separate(SAMPLEID, c("SAMPLEID", "else"), sep = "_L001_")
#
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID) #This is a base R function that searches for a pattern in a string and replaces all occurrences with something else.
paths_run$FASTQ <- gsub("./","", paths_run$FASTQ)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$FASTQ, sep="/") #make a new column with the full path

forward <- paths_run %>%
        filter(grepl("R1_001.fastq.gz", FASTQ)) %>% 
        select(SAMPLEID, `forward-absolute-filepath` = FULL_PATH)

# Write manifest file
manifest <- paths_run %>%
        filter(grepl("R2_001.fastq.gz", FASTQ)) %>%
        select(SAMPLEID, `reverse-absolute-filepath` = FULL_PATH) %>%
        right_join(forward) %>%
        select('sample-id' = SAMPLEID, `forward-absolute-filepath`, `reverse-absolute-filepath`) #join forward and reverse file paths

# Write output as a manifest file
write.table(manifest, file = "/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R/manifest_515F_926R.tsv", quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t") #write manifest file
```

## Step 4: Run manifest file
### 515F/806R
```
#enter scavenger: srun -p scavenger --time=04:00:00 --ntasks-per-node 1 --mem=10gb --pty bash

#Start R conda environment (you will have to make your own R environment using conda if you don't have one already), see docs here on how to install one: https://anaconda.org/conda-forge/r-base
conda activate R_environment

#run the file: 
Rscript manifest_515F_806R.R

#check your manifest file
head(manifest)
nrow(manifest)

#fix manifest file to work with the singularity docker: 
sed 's|/proj/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R/|/data/|g' manifest_515F_806R.tsv > manifest_515F_806R_fixed.tsv

#this fix ensures that singularity can bind your WORKDIR with data/ in the next step. Essentially this means that it will change the absolute file path of your WORKDIR into a relative file path (data/)
```

### 515F/926R
```
#enter scavenger: srun -p scavenger --time=04:00:00 --ntasks-per-node 1 --mem=10gb --pty bash

#Start R conda environment (you will have to make your own R environment using conda if you don't have one already), see docs here on how to install one: https://anaconda.org/conda-forge/r-base
conda activate R_environment

#run the file: 
Rscript manifest_515F_926R.R

#check your manifest file
head(manifest)
nrow(manifest)

#fix manifest file to work with the singularity docker: 
sed 's|/proj/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R/|/data/|g' manifest_515F_926R.tsv > manifest_515F_9266R_fixed.tsv

#this fix ensures that singularity can bind your WORKDIR with data/ in the next step. Essentially this means that it will change the absolute file path of your WORKDIR into a relative file path (data/)
```






# References: 
1. Sheik CS, Reese BK, Twing KI, Sylvan JB, Grim SL, Schrenk MO, Sogin ML, Colwell FS. Identification and removal of contaminant sequences from ribosomal gene databases: lessons from the census of deep life. Frontiers in microbiology. 2018 Apr 30;9:840.
2. Zhou J, Bruns MA, Tiedje JM. DNA recovery from soils of diverse composition. Applied and environmental microbiology. 1996 Feb;62(2):316-22.
3. Crump BC, Kling GW, Bahr M, Hobbie JE. Bacterioplankton community shifts in an arctic lake correlate with seasonal changes in organic matter source. Applied and Environmental Microbiology. 2003 Apr;69(4):2253-68.
4. Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome. 2018 Dec 17;6(1):226.
5. Parada, A. E., Needham, D. M., & Fuhrman, J. A. (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. Environmental Microbiology, 18(5), 1403–1414. http://doi.org/10.1111/1462-2920.13023
6. Apprill, A., McNally, S., Parsons, R., & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. Aquatic Microbial Ecology, 75(2), 129–137. http://doi.org/10.3354/ame01753
7. Quince, C., Lanzen, A., Davenport, R.J., & Turnbaugh, P.J. (2011) Removing noise from pyrosequenced amplicons. BMC Bioinformatics 12: 38. https://doi.org/10.1186/1471-2105-12-38
