# Qiime2_Mariana_forearc_16S_Amplicons
Updated 9.17.25

How to process 16S rRNA amplicon sequences using Qiime2. These samples were very low-biomass, from pH 12.5 fluid. Contaminant removal follows methods in Sheik et al., 2018.

45 mL of serpentinized fluid incubated with <sup>13</sup>CH<sub>4</sub> gas to track methanotrophy was extracted using phenol:chloroform methods described in Zhou et al., 1996 and Crump et al., 2003 and modified by Dr. Caroline Fortunato and Dr. Amy Smith. DNA was amplified using two EMP primer sets (below). This pipeline includes code for both of them. This code in modified from Dr. Sarah Hu's github: https://github.com/shu251/qiime2-2024/tree/main?tab=readme-ov-file.

## Primer sets: 
- 515F (Parada et al., 2016), 806R (Apprill et al., 2015)
- 515F (Parada et al., 2016), 926R (Quince et al., 2011)

## Adapters: Nextera kit 
- Forward: CTGTCTCTTATACACATCT
- Reverse: CTGTCTCTTATACACATCT

## Sequencing
Illumina Miseq using an Illumina PE250 Nano (500 cycles) reagent kit.

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
This creates a file with the sampleID's in one column and their paths. 
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

## Step 5: create slurm script to import your sequences as a Qiime2 artifact
### 515F/806R
```
nano qiime2_import_806.sh

#write slurm script (cp into nano doc)
#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=qiime2_import806R                 # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=4                            # Number of CPU cores per task
#SBATCH --mem=10gb                                   # Job memory request
#SBATCH --time=24:00:00								               # Time limit hrs:min:sec
#SBATCH --output=qiime2_import806R.log               # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module (adjust if needed on your system)
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R"

# QIIME2 import (manifest-based paired-end FASTQ)
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /data/manifest_515F_806R_fixed.tsv \
  --output-path /data/paired-end-input_806R.qza \
  --input-format PairedEndFastqManifestPhred33V2

# QIIME2 visualization of imported data
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime demux summarize \
  --i-data /data/paired-end-input_806R.qza \
  --o-visualization /data/paired-end-input_806R.qzv
```

### 515F/926R
```
#create slurm script to import your sequences
nano qiime2_import_926.sh

#write slurm script (cp into nano doc)
#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=qiime2_import926R                 # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=4                            # Number of CPU cores per task
#SBATCH --mem=10gb                                   # Job memory request
#SBATCH --time=24:00:00								               # Time limit hrs:min:sec
#SBATCH --output=qiime2_import926R.log               # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module (adjust if needed on your system)
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R"

# QIIME2 import (manifest-based paired-end FASTQ)
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /data/manifest_515F_926R_fixed.tsv \
  --output-path /data/paired-end-input_926R.qza \
  --input-format PairedEndFastqManifestPhred33V2

# QIIME2 visualization of imported data
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime demux summarize \
  --i-data /data/paired-end-input_926R.qza \
  --o-visualization /data/paired-end-input_926R.qzv
```

## Step 6: Trim primer and adapter sequences using Cutadapt in Qiime2
### 515F/806R
```
#create slurm script to import your sequences
#script name: qiime_primer_806R.slurm
#write slurm script (cp into nano doc)

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=qiime2_primer_806R                # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=4                            # Number of CPU cores per task
#SBATCH --mem=10gb                                   # Job memory request
#SBATCH --time=24:00:00								               # Time limit hrs:min:sec
#SBATCH --output=qiime2_primer_806R.log              # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module (adjust if needed on your system)
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R"

# Remove primer sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /data/paired-end-input_806R.qza \
  --p-cores $SLURM_CPUS_PER_TASK \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-error-rate 0.1 \
  --p-overlap 3 \
   --p-adapter-f  \
   --p-adapter-r CTGTCTCTTATACACATCT \
  --p-match-adapter-wildcards \
  --o-trimmed-sequences /data/paired-end-input-trimmed_806R.qza

# Summarize trimmed sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime demux summarize \
  --i-data /data/paired-end-input-trimmed_806R.qza \
  --o-visualization /data/paired-end-input-trimmed_806R.qzv
```

### 515F/926R
```
#create slurm script to import your sequences
#script name: qiime_primer_926R.slurm
#write slurm script (cp into nano doc)

#!/bin/bash
#SBATCH --partition=compute                          # Queue selection
#SBATCH --job-name=qiime2_primer_926R                # Job name
#SBATCH --mail-type=ALL                              # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu               # Where to send mail
#SBATCH --ntasks=1                                   # Run a single task
#SBATCH --cpus-per-task=4                            # Number of CPU cores per task
#SBATCH --mem=10gb                                   # Job memory request
#SBATCH --time=24:00:00								               # Time limit hrs:min:sec
#SBATCH --output=qiime2_primer_926R.log              # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module (adjust if needed on your system)
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R"

# Remove primer sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /data/paired-end-input_926R.qza \
  --p-cores $SLURM_CPUS_PER_TASK \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r CCGYCAATTYMTTTRAGTTT \
  --p-error-rate 0.1 \
  --p-overlap 3 \
   --p-adapter-f CTGTCTCTTATACACATCT \
   --p-adapter-r CTGTCTCTTATACACATCT \
  --p-match-adapter-wildcards \
  --o-trimmed-sequences /data/paired-end-input-trimmed_926R.qza

# Summarize trimmed sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime demux summarize \
  --i-data /data/paired-end-input-trimmed_926R.qza \
  --o-visualization /data/paired-end-input-trimmed_926R.qzv
```

## Step 7: View trimmed files on Qiime2View and decided where to truncate sequences
```
Qiime2View: https://view.qiime2.org/

Forum for truncation info (I found this very helpful!)
https://forum.qiime2.org/t/deciding-where-to-trim-in-dada2-on-paired-end-reads/18762
```

## Step 8: Remove primer and adapter sequences and run DADA2 to call ASV's
### 515F/806R
```
#slurm script name: dada2_denoise_806R.sh

#!/bin/bash
#SBATCH --partition=compute                                                                  # Queue selection
#SBATCH --job-name=qiime2_dada2_denoise_806R                                                 # Job name
#SBATCH --mail-type=ALL                                                                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu                                                       # Where to send mail
#SBATCH --ntasks=1                                                                           # Run a single task
#SBATCH --cpus-per-task=4                                                                    # Number of CPU cores per task
#SBATCH --mem=25gb                                                                           # Job memory request
#SBATCH --time=04:00:00                                                                      # Time limit hrs:min:sec
#SBATCH --output=qiime2_dada2_denoise_806R.log                                               # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R"
export MPLCONFIGDIR=/tmp/mplconfig && mkdir -p $MPLCONFIGDIR

# Run DADA2 denoising
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /data/paired-end-input-trimmed_806R.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-min-overlap 12 \
 --p-pooling-method independent \
  --p-n-reads-learn 1000000 \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-chimera-method pooled \
  --o-table /data/paired-end-input-asv-table_806R.qza \
  --o-representative-sequences /data/paired-end-input-asvs-ref-seqs_806R.qza \
  --o-denoising-stats /data/paired-end-input-dada2-stats_806R.qza

# Export ASV table
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools export \
  --input-path /data/paired-end-input-asv-table_806R.qza \
  --output-path /data/

# Convert BIOM to TSV
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
biom convert \
  -i /data/feature-table.biom \
  -o /data/samples-asv-table_806R.tsv \
  --to-tsv

# Tabulate DADA2 stats
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime metadata tabulate \
  --m-input-file /data/paired-end-input-dada2-stats_806R.qza \
  --o-visualization /data/paired-end-input-dada2-stats_806R.qzv
```

### 515F/926R
```
#!/bin/bash
#SBATCH --partition=compute                                                                  # Queue selection
#SBATCH --job-name=qiime2_no_trunc_926R                                                      # Job name
#SBATCH --mail-type=ALL                                                                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu                                                       # Where to send mail
#SBATCH --ntasks=1                                                                           # Run a single task
#SBATCH --cpus-per-task=4                                                                    # Number of CPU cores per task
#SBATCH --mem=25gb                                                                           # Job memory request
#SBATCH --time=04:00:00                                                                      # Time limit hrs:min:sec
#SBATCH --output=qiime2_no_trunc_926R.log                                                    # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module 
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R"
export MPLCONFIGDIR=/tmp/mplconfig && mkdir -p $MPLCONFIGDIR

# Remove primer sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /data/paired-end-input_926R.qza \
  --p-cores $SLURM_CPUS_PER_TASK \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r CCGYCAATTYMTTTRAGTTT \
  --p-adapter-f CTGTCTCTTATACACATCT \
  --p-adapter-r CTGTCTCTTATACACATCT \
  --p-error-rate 0.1 \
  --p-overlap 10 \
  --p-match-adapter-wildcards True \
  --o-trimmed-sequences /data/paired-end-input-trimmed_926R.qza

# Summarize trimmed sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime demux summarize \
  --i-data /data/paired-end-input-trimmed_926R.qza \
  --o-visualization /data/paired-end-input-trimmed_926R.qzv

# Run DADA2 denoising
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /data/paired-end-input-trimmed_926R.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-min-overlap 12 \
 --p-pooling-method independent \
  --p-n-reads-learn 1000000 \
  --p-n-threads $SLURM_CPUS_PER_TASK \
  --p-chimera-method pooled \
  --o-table /data/paired-end-input-asv-table_926R.qza \
  --o-representative-sequences /data/paired-end-input-asvs-ref-seqs_926R.qza \
  --o-denoising-stats /data/paired-end-input-dada2-stats_926R.qza

# Export ASV table
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools export \
  --input-path /data/paired-end-input-asv-table_926R.qza \
  --output-path /data/

# Convert BIOM to TSV
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
biom convert \
  -i /data/feature-table.biom \
  -o /data/samples-asv-table_926R.tsv \
  --to-tsv

# Tabulate DADA2 stats
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime metadata tabulate \
  --m-input-file /data/paired-end-input-dada2-stats_926R.qza \
  --o-visualization /data/paired-end-input-dada2-stats_926R.qzv
```
### Output explanation: 
	•	samples-asv-table.tsv: ASV table that is tab-delimited. It should show your samples as columns, ASVs (or FeatureIDs) as rows, and numbers will represent the number of sequences assigned to the respective ASV and sample.
	•	paired-end-input-asvs-ref-seqs.qza: Reference sequences for the output ASVs. Use this below to assign taxonomy.
	•	paired-end-input-dada2-stats.qzv: Output visualization of the DADA2 stats. Import this to the QIIME2 viewer to learn more about how your sequences did.


## Step 9: Format Silva 16S rRNA db v138.2 for use in QIIME2 using the Rescript plugin. Same formatting for both primer sets. 
```
#script name: qiime2_silva.slurm

#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=qiime2_silva_806R
#SBATCH --mail-type=ALL
#SBATCH --mail-user=selkassas@whoi.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --output=qiime2_silva_806R.log

export OMP_NUM_THREADS=4

module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R/"
export MPLCONFIGDIR=/tmp/mplconfig && mkdir -p $MPLCONFIGDIR

# Get SILVA database
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime rescript get-silva-data \
  --p-version '138.1’ \
  --p-target 'SSURef_NR99' \
  --o-silva-sequences /data/silva-138.1-ssu-nr99-rna-seqs_806R.qza \
  --o-silva-taxonomy /data/silva-138.1-ssu-nr99-tax_806R.qza

# Reverse-transcribe RNA to DNA
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime rescript reverse-transcribe \
  --i-rna-sequences /data/silva-138.1-ssu-nr99-rna-seqs_806R.qza \
  --o-dna-sequences /data/silva-138.1-ssu-nr99-seqs_806R.qza

# Cull low-quality sequences
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime rescript cull-seqs \
  --i-sequences /data/silva-138.1-ssu-nr99-seqs_806R.qza \
  --o-clean-sequences /data/silva-138.1-ssu-nr99-seqs-cleaned_806R.qza

# Filter sequences by taxon-specific length thresholds
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences /data/silva-138.1-ssu-nr99-seqs-cleaned_806R.qza \
  --i-taxonomy /data/silva-138.1-ssu-nr99-tax_806R.qza \
  --p-labels Archaea Bacteria Eukaryota \
  --p-min-lens 900 1200 1400 \
  --o-filtered-seqs /data/silva-138.1-ssu-nr99-seqs-filt_806R.qza \
  --o-discarded-seqs /data/silva-138.1-ssu-nr99-seqs-discard_806R.qza

#submit to slurm
sbatch qiime2_silva.slurm

```

## Step 10: Classify ASV's with VSearch and compute the core microbiome
This function, qiime feature-table core-features, allows you to see if there are taxa common to all samples. Especially if you have samples from different sources, this can tell you if there is a contaminant common to all sequences. You will then be able to remove this contaminant before moving forward with the pipeline. 
Use this file and the taxonomy bar plot to scope out potential contaminants. 
### 515F/806R
```
#slurm script: qiim2_classify_806.sh

#!/bin/bash
#SBATCH --partition=compute                                                                  # Queue selection
#SBATCH --job-name=qiime2_classify_806R                                                      # Job name
#SBATCH --mail-type=ALL                                                                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu                                                       # Where to send mail
#SBATCH --ntasks=1                                                                           # Run a single task
#SBATCH --cpus-per-task=4                                                                    # Number of CPU cores per task
#SBATCH --mem=25gb                                                                           # Job memory request
#SBATCH --time=04:00:00                                                                      # Time limit hrs:min:sec
#SBATCH --output=qiime2_classify_806R.log                                                    # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_806R"
export MPLCONFIGDIR=/tmp/mplconfig && mkdir -p $MPLCONFIGDIR

# Taxonomic classification using VSEARCH
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-classifier classify-consensus-vsearch \
  --i-query /data/paired-end-input-asvs-ref-seqs_806R.qza \
  --i-reference-reads /data/silva-138.1-ssu-nr99-seqs-filt_806R.qza \
  --i-reference-taxonomy /data/silva-138.1-ssu-nr99-tax_806R.qza \
  --o-classification /data/classification_806R.qza \
  --o-search-results /data/search_results_806R.qza \
  --p-threads $SLURM_CPUS_PER_TASK \
  --p-maxaccepts 10 \
  --p-perc-identity 0.90 \
  --p-min-consensus 0.70 \
  --p-query-cov 0.90

# Export classification results to TSV
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools export \
  --input-path /data/classification_806R.qza \
  --output-path /data/

# Create a .qzv visualization of the taxonomy table
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime metadata tabulate \
  --m-input-file /data/classification_806R.qza \
  --o-visualization /data/classification_806R.qzv

# Taxonomy barplots
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime taxa barplot \
  --i-table /data/paired-end-input-asv-table_806R.qza \
  --i-taxonomy /data/classification_806R.qza \
  --o-visualization /data/taxonomy-barplot_806R.qzv

# choose a single threshold (e.g., 80% of samples):
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-table core-features \
  --i-table /data/paired-end-input-asv-table_806R.qza \
  --p-min-fraction 0.8 \
  --p-max-fraction 0.8 \
  --o-visualization /data/core-features_asv_80_806R.qzv

# collapse to genus (level 6: k=1 … g=6, s=7)
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime taxa collapse \
  --i-table /data/paired-end-input-asv-table_806R.qza \
  --i-taxonomy /data/classification_806R.qza \
  --p-level 6 \
  --o-collapsed-table /data/table_genus_806R.qza

# core genera at 70–100% with 7 steps (0.70, 0.75, …, 1.00):
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-table core-features \
  --i-table /data/table_genus_806R.qza \
  --p-min-fraction 0.70 \
  --p-max-fraction 1.0 \
  --p-steps 7 \
  --o-visualization /data/core-features_genus_70_100_806R.qzv
```

### 515F/926R:
```
#slurm script: qiim2_classify_926.sh

#!/bin/bash
#SBATCH --partition=compute                                                                  # Queue selection
#SBATCH --job-name=qiime2_classify_926R                                                      # Job name
#SBATCH --mail-type=ALL                                                                      # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu                                                       # Where to send mail
#SBATCH --ntasks=1                                                                           # Run a single task
#SBATCH --cpus-per-task=4                                                                    # Number of CPU cores per task
#SBATCH --mem=25gb                                                                           # Job memory request
#SBATCH --time=04:00:00                                                                      # Time limit hrs:min:sec
#SBATCH --output=qiime2_classify_926R.log                                                    # Job log name
export OMP_NUM_THREADS=4

# Load Singularity module
module load singularity/3.7

CONTAINER=/vortexfs1/home/selkassas/amplicon_2024.2.sif
WORKDIR="/vortexfs1/omics/huber/selkassas/Mariana_amplicon_sequences_2/data_515F_926R"
export MPLCONFIGDIR=/tmp/mplconfig && mkdir -p $MPLCONFIGDIR

# Taxonomic classification using VSEARCH
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-classifier classify-consensus-vsearch \
  --i-query /data/paired-end-input-asvs-ref-seqs_926R.qza \
  --i-reference-reads /data/silva-138.1-ssu-nr99-seqs-filt_926R.qza \
  --i-reference-taxonomy /data/silva-138.1-ssu-nr99-tax_926R.qza \
  --o-classification /data/classification_926R.qza \
  --o-search-results /data/search_results_926R.qza \
  --p-threads $SLURM_CPUS_PER_TASK \
  --p-maxaccepts 10 \
  --p-perc-identity 0.90 \
  --p-min-consensus 0.70 \
  --p-query-cov 0.90

# Export classification results to TSV
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime tools export \
  --input-path /data/classification_926R.qza \
  --output-path /data/

# Create a .qzv visualization of the taxonomy table
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime metadata tabulate \
  --m-input-file /data/classification_9266R.qza \
  --o-visualization /data/classification_926R.qzv

# Taxonomy barplots
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime taxa barplot \
  --i-table /data/paired-end-input-asv-table_926R.qza \
  --i-taxonomy /data/classification_926R.qza \
  --o-visualization /data/taxonomy-barplot_926R.qzv

# choose a single threshold (e.g., 80% of samples):
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-table core-features \
  --i-table /data/paired-end-input-asv-table_926R.qza \
  --p-min-fraction 0.8 \
  --p-max-fraction 0.8 \
  --o-visualization /data/core-features_asv_80_926R.qzv

# collapse to genus (level 6: k=1 … g=6, s=7)
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime taxa collapse \
  --i-table /data/paired-end-input-asv-table_926R.qza \
  --i-taxonomy /data/classification_926R.qza \
  --p-level 6 \
  --o-collapsed-table /data/table_genus_926R.qza

# core genera at 70–100% with 7 steps (0.70, 0.75, …, 1.00):
singularity exec --bind ${WORKDIR}:/data ${CONTAINER} \
qiime feature-table core-features \
  --i-table /data/table_genus_926R.qza \
  --p-min-fraction 0.70 \
  --p-max-fraction 1.0 \
  --p-steps 7 \
  --o-visualization /data/core-features_genus_70_100_926R.qzv
```

## Step 11: Create metadata file for your samples, add DADA2 stats to it
It has to be a very specific format, see link below for guide. 
```
https://docs.qiime2.org/2024.10/tutorials/metadata/

# To grab the DADA2 stats in table form
Qiime2 View —> paired-end-input-dada2-stats_926R.qzv —> download metadata table

# Make sure the first line "#q2:…" is deleted because this will be imported as a Qiime2 object, and Qiime2 won’t be 	able to do it correctly if that line is still in there. 
```

## Step 12: 
View the classification.qzv on Qiime2 View. download it as a .tsv. This file has the ASV ID’s matched with the taxonomy. 

# References: 
1. Sheik CS, Reese BK, Twing KI, Sylvan JB, Grim SL, Schrenk MO, Sogin ML, Colwell FS. Identification and removal of contaminant sequences from ribosomal gene databases: lessons from the census of deep life. Frontiers in microbiology. 2018 Apr 30;9:840.
2. Zhou J, Bruns MA, Tiedje JM. DNA recovery from soils of diverse composition. Applied and environmental microbiology. 1996 Feb;62(2):316-22.
3. Crump BC, Kling GW, Bahr M, Hobbie JE. Bacterioplankton community shifts in an arctic lake correlate with seasonal changes in organic matter source. Applied and Environmental Microbiology. 2003 Apr;69(4):2253-68.
4. Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome. 2018 Dec 17;6(1):226.
5. Parada, A. E., Needham, D. M., & Fuhrman, J. A. (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. Environmental Microbiology, 18(5), 1403–1414. http://doi.org/10.1111/1462-2920.13023
6. Apprill, A., McNally, S., Parsons, R., & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. Aquatic Microbial Ecology, 75(2), 129–137. http://doi.org/10.3354/ame01753
7. Quince, C., Lanzen, A., Davenport, R.J., & Turnbaugh, P.J. (2011) Removing noise from pyrosequenced amplicons. BMC Bioinformatics 12: 38. https://doi.org/10.1186/1471-2105-12-38
