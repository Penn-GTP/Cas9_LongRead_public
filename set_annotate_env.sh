#!/bin/bash

## Load required modules on HPC
module load samtools/1.11 # SAMtools
module load R/4.2 # R

# set envs
export PATH=$PATH:/project/gtplab/apps/bin:/project/gtplab/apps/subread-2.0.3-source/bin
