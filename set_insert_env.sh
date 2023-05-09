#!/bin/bash

## Load required modules on HPC
module load samtools/1.11 # SAMtools
module load java/openjdk-1.8.0 # java for running picard

# set envs
PATH=$PATH:/project/gtplab/apps/bin # gtplab apps
PATH=$PATH:/project/gtplab/apps/minimap2 # minimap2

export PATH
