# Cas9_LongRead pipeline manual

This pipeline guideline was initiated by: Qi Zheng (<zhengqi@pennmedicine.upenn.edu>).

This manual is for the new Cas9_LongRead pipeline developed based on previous unpublished procedures.

We thank Lili Wang (<liliwang@upenn.edu>), Claude Warzecha (<claude@upenn.edu>), Xiang Huang, and Jenny Greig (<greg@upenn.edu>) for their help.

Cas9 LongRead sequencing assay is a method to based on CRISP/Cas9 crRNA enrichment of targeted genomic DNAs and Long-read sequencing.

---

## Citations

---

## Dependency and pre-requirement
You need working Perl and Java distributions and available in your *PATH* to run this pipeline.
Both Perl and Java are generally available under most model Operation Systems,
including Windows, Mac OS X, and Unix/Linux, and is likely pre-install on most Unix/Linux based systems.

Besides, this pipeline also depends the following programs/tools pre-installed and available in your *PATH*.

If you are running this pipeline on PMACS HPC, all the following dependencies are pre-installed and pre-configured,
so you can ignore this step and jump to the next step.

1. **Samtools** - basic tool for manipulating standard alignment files in SAM/BAM format.
You can download and install **Samtools** by following the instructions at <https://www.htslib.org/>.

2. **Minimap2** - Default Long-read read aligner used by this pipeline.
You can download and install **Minimap2** by following the instructions at <https://github.com/lh3/minimap2>.

3. **Picard tools** - Picard is a set of Java based command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
This pipeline ships a JAR copy of **Picard tools** with it, but you can replace it with your own copy, you can achieve this by
simply download and replace the latest `picard.jar` file from https://broadinstitute.github.io/picard/ in this pipeline.

4. **R** - R is a free software environment for statistical computing and graphics.
You can download and install **R** by following the instructions at <https://www.r-project.org/>.

5. **Subread** - Subread package: high-performance read alignment, quantification and mutation discovery.
You can download and install **Subread** package by following the instructions at <https://subread.sourceforge.net/>.

 
## Getting Started

This pipeline's scripts only "prepare" the running commands (cmds) and output them into the shell/bash files.
After running these "prepare" scripts, you will need to run them either by directly execute the output files,
or "submitting" them via an LSF system such as `bsub` that is available on PMACS/HPC.
The "prepare" scripts will also check whether the output files already exist,
and if yes, it will not override the existing files and instead generate "commented-out" commands that you can feel free to enable/disable them manually.

**Note:** If you are using this pipeine in an HPC environment, be sure to do any of these steps in an interactive node,
initiated by typing `bsub -Is "bash"`.

### Create the experimental design config file

You need to create a special tab-delimited config file to specify all the experimental designs for a given project (**PROJECT_ID**).
You should do this by copying the template file `Template_Cas9_LongRead_experimental_design.conf` in this pipeline into your working directory,
and fill/edit it with a text-editor or Excel/Libre-Office,
and save it in tab-delimited format to your working directory (assuming it is saved as `PROJECT_ID_EXP_DESIGN.conf`)
 
Each project-specific config file contains global options and per-sample local options,
both are explained in the "comment lines" at the beginning of the template,
while the global options are in all upper-cases, and the values are given in **GLOBAL\_OPT=GLOBAL\_VAL** format,
and local options are given in corresponding table cells.

**Note:** If you are setting any of the working dir options below to values other than the current dir (./),
you will have to create or link the directories manually before running this pipeline.

#### Global options

- **NUM\_PROC**: # of processors (cores/threads) to use [default 8]
- **BASE\_DIR**: base directory for generating major result files (sequences, alignments, tables, figures, etc.) [default ./]
- **WORK\_DIR**: working directory for all intermediate files [default WORK/]
- **SCRIPT\_DIR**: path to this pipeline feel free to link it to your working directory [default scripts/]
- **VEC\_DIR**: path to the AAV vector/plasmid map files in GenBank format, also used to build vector databases [default AAV_vec/]
- **NGS\_ALIGNER**: NGS aligner name to use, only supports 'minimap2' [default minimap2]
- **FEATURE\_TAG**: tag used to search for functional donor features in customized GFF3 annotation files [default label]
- **MIN\_COVER\_RATIO**: minimum cover-ratio (aligned bases / feature-length) required for defining a functional insert [default 0.9]
- **ARM\_KEY**: keyword/regex to search the feature labels to determine whether it is an ARM to determine HDR vs. NHEJ inserts [default (?i:ARM)|shHDR]
- **ARM\_MIN\_RATIO**: minimum cover-ratio for 5'/3' HDR at end of long-reads to determine HDR reads [default 0.1]

#### Per-sample options
- **sample\_name**: unique sample name
- **longread\_tech**: LongRead technology, must be one of 'ont' (Oxford Nanopore Technology), 'pb' (PacBio CLR) or 'hifi' (PacBio HiFi) [default ont]
- **read\_fastq**: per-sample contatenated LongRead FASTQ read file, supports .gz or .bz2 compressed files
- **ref\_db**: path to ref genome database
- **ref\_gff**: path to ref genome GFF3/GTF annotation file(s), multiple files can be separated by comma
- **enrich\_bed**: path to the BED file of the known Cas9/crRNA recognition site on the ref genome
- **target_bed**: path to the BED file of the known ARCUS nuclease target site on the ref genome
- **nuclease\_gb**: path to the nuclease vector sequence in GenBank format, leave blank if not used
- **donor\_gb**: path to the donor vector sequence in GenBank format
- **trans\_gb**: path to the trans-plasmid GenBank file used in vector production, optional
- **helper\_gb**: path to the helper-plasmid GenBank file used in vector production, optional
- **ref2\_db**: path to additional reference genome/database, ignored if empty
- **ref2\_gff**: path to ref2 genome GFF3/GTF annotation file(s), ignored if empty, multiple files can be separated by comma
- **vec2\_db**: path to additional virus/vector database/seq, ignored if empty
- **vec2\_gff**: path to vec2 genome GFF3/GTF annotation file(s), ignored if empty, multiple files can be separated by comma
- **ref\_map\_opts**: additional ref mapping options you want to invoke the NGS aligner, i.e. '-N 50 -p 0.1' to increase the sensitivity for finding secondary alignments
- **min\_mapQ**: minimum mapping quality (mapQ) required for the ref mapped reads (recommend 30)
- **min\_insert**: minimum insert size required at the known ARCUS insert site for the long reads  [default 20]
- **max\_dist**: maximum distance to known ARCUS target site for merging multiple-inserts on the same long read; suggest to use the homologue arm length [default 520]
- **vec\_map\_opts**: additional nuclease/donor vec mapping options you want to invoke the NGS aligner, i.e. '-p 0.1'
- **ref2\_map\_opts**: any additonal ref mapping options you want to invoke the NGS aligner, i.e. '-p 0.1'
- **vec2\_map\_opts**: any additonal options you want to use to invoke the NGS aligner, i.e. '-p 0.1'
- **ref\_anno\_opts**: additional options for generating ref-based insert annotations, i.e. '--display Name'
- **vec\_anno\_opts**: additional options for generating vec-based insert annotations, i.e. '--display label'
- **ref2\_anno\_opts**: additional options for generating ref2-based insert annotations, i.e. '--display Name'
- **vec2\_anno\_opts**: additional options for generating vec2-based insert annotations, i.e. '--display product'
- **nuclease\_count\_opts**: additional nuclease vector mapped featureCount options you want to invoke the featureCount program, i.e. '-g Name'
- **donor\_count\_opts**: additional donor vector mapped featureCount options you want to invoke the featureCount program, i.e. '-g Name'
- **functional\_feature**: vector feature(s) (given in the label) required to define a functional donor target/off longread, seperated by '|', i.e. 'transgene'

---

## Step REF -- Build Reference Database files for the project
- Prepare REF cmds by running:
`SCRIPT_DIR/prepare_ref_cmd.pl PROJECT_ID_EXP_DESIGN.conf ref.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `ref.sh`
   
- Run REF cmds in `ref.sh` by running:
   - On a Linux cluster: `bsub -J REF -o ref.log ./ref.sh`
   - On a regular Linux server: `./ref.sh > ref.log 2>&1`
   
- Output for each PROJECT:
   - REF database files in `VEC_DIR`, including a sequence file (FASTA format), a region file (BED format), a feature annotaiton file (GFF3 format) for nuclease, donor, trans-plasmid, and helper-plasmid, if provided in config file
   - Combined vectors' sequence file (FASTA format) and feature annotation file (GFF3 format)

## Step MAP -- Map long-reads to host/ref genome
- Prepare MAP cmds by running:
`SCRIPT_DIR/prepare_map_cmd.pl PROJECT_ID_EXP_DESIGN.conf map.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `map.sh`
   
- Run MAP cmds in `map.sh` by running:
   - On a Linux cluster: `bsub -J MAP -o map.log -n 24 -M 25600 ./map.sh`
   - On a regular Linux server: `./map.sh > map.log 2>&1`
   
- Output for each **SAMPLE**:
   - Host reference (ref) mapped alignment: `BASE_DIR/SAMPLE_ref_map_filtered_sorted.bam`
   - Ref mapped alignment in Cas9 enrichment region: `BASE_DIR/SAMPLE_ref_map_enrich_sorted.bam`
   - Ref mapped alignment covering nuclease-vector target region: `BASE_DIR/SAMPLE_ref_map_target_sorted.bam`
          
## Step INSERT -- Extract ref inserts and map to vec genomes 
- Prepare INSERT cmds by running:
`SCRIPT_DIR/prepare_insert_cmd.pl PROJECT_ID_EXP_DESIGN.conf insert.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `insert.sh`
   
- Run INSERT cmds in `insert.sh` by running:
   - On a Linux cluster: `bsub -J FILTER -o filter.log -n 24 -M 25600 ./filter.sh`
   - On a regular Linux server: `./filter.sh > filter.log 2>&1`
   
- Output for each **SAMPLE**:
   - On-target  insert sequences and info extracted from long-reads (merged if close enough): `BASE_DIR/SAMPLE_target_insert_seq.fastq`, `BASE_DIR/SAMPLE_target_insert_seq.fasta`, `BASE_DIR/SAMPLE_target_insert_info.tsv`
   - Off-target insert sequences and info extracted from long-reads (merged if close enough): `BASE_DIR/SAMPLE_off_insert_seq.fastq`, `BASE_DIR/SAMPLE_off_insert_seq.fasta`, `BASE_DIR/SAMPLE_off_insert_info.tsv`
   - Vec mapped alignment of on-target inserts: `BASE_DIR/SAMPLE_target_insert_vec_filtered_sorted.bam`   
   - Ref2 mapped alignment of on-target inserts: `BASE_DIR/SAMPLE_target_insert_ref2_filtered_sorted.bam`, if ref2 provided   
   - Vec2 mapped alignment of on-target inserts: `BASE_DIR/SAMPLE_target_insert_vec2_filtered_sorted.bam`, if vec2 provided
   - Vec mapped alignment of off-target inserts: `BASE_DIR/SAMPLE_off_insert_vec_filtered_sorted.bam`   
   - Ref2 mapped alignment of off-target inserts: `BASE_DIR/SAMPLE_off_insert_ref2_filtered_sorted.bam`, if ref2 provided   
   - Vec2 mapped alignment of off-target inserts: `BASE_DIR/SAMPLE_off_insert_vec2_filtered_sorted.bam`, if vec2 provided
   - Vec mapped sequences of on-target inserts: `BASE_DIR/SAMPLE_target_insert_vec_map_seq.fasta`   
   - Vec mapped sequences of off-target inserts: `BASE_DIR/SAMPLE_off_insert_vec_map_seq.fasta`   
   - Ref mapped alignment of vec-mappable on-target inserts: `BASE_DIR/SAMPLE_target_insert_ref_map_sorted.bam`   
   - Ref mapped alignment of vec-mappable off-target inserts: `BASE_DIR/SAMPLE_off_insert_ref_map_sorted.bam`   
   - Ref mapped, vec unmapped (genomic) alignment of on-target inserts: `BASE_DIR/SAMPLE_target_genomic_ref_map_sorted.bam`   
   - Ref mapped, vec unmapped (genomic) sequence of on-target inserts: `BASE_DIR/SAMPLE_target_genomic_ref_map_seq.fasta`

## Step ANNOTATE -- annotate called inserts
- Prepare ANNO cmds by running:
`SCRIPT_DIR/prepare_annotate_cmd.pl PROJECT_ID_EXP_DESIGN.conf annotate.sh`
   - **Input**: `PROJECT_ID_EXP_DESIGN.conf`
   - **Output**: `annotate.sh`
   
- Run ANNOTATE cmds in `annotate.sh` by running:
   - On a Linux cluster: `bsub -J ANNOTATE -o annotate.log -n 4 -M 48000 ./annotate.sh`
   - On a regular Linux server: `./annotate.sh > annotate.log 2>&1`
   
- Output for each **SAMPLE**:
   - size distribution of on-target inserts: `BASE_DIR/SAMPLE_target_insert_size_distrib.pdf`
   - size distribution of off-target inserts: `BASE_DIR/SAMPLE_off_insert_size_distrib.pdf`
   - per long-read on-target insert feature annotation in standard GFF3 format (readable for IGV) for nuclease|donor|trans|helper|combined vectors/plasmids: `BASE_DIR/SAMPLE_target_insert_(nulcease|donor|trans|helper)_vec_anno.gff3`
   - per long-read off-target insert feature annotation in standard GFF3 format (readable for IGV) for nuclease|donor|trans|helper|combined vectors/plasmids: `BASE_DIR/SAMPLE_off_insert_(nulcease|donor|trans|helper)_vec_anno.gff3`
   - per long-read on-target insert feature annotation in standard GFF3 format (readable for IGV) for ref2: `BASE_DIR/SAMPLE_target_insert_ref2_anno.gff3`, if provided
   - per long-read off-target insert feature annotation in standard GFF3 format (readable for IGV) for ref2: `BASE_DIR/SAMPLE_off_insert_ref2_anno.gff3`, if provided
   - per long-read on-target insert feature annotation in standard GFF3 format (readable for IGV) for vec2: `BASE_DIR/SAMPLE_target_insert_vec2_anno.gff3`, if provided
   - per long-read off-target insert feature annotation in standard GFF3 format (readable for IGV) for vec2: `BASE_DIR/SAMPLE_off_insert_vec2_anno.gff3`, if provided
   - feature-count table of nuclease-vector mapped long-reads for all on-target|off-target inserts: `BASE_DIR/SAMPLE_insert_nuclease_vec_feature_count.tsv`
   - feature-count table of donor-vector mapped long-reads for all on-target|off-target inserts: `BASE_DIR/SAMPLE_insert_donor_vec_feature_count.tsv`
   - per long-read on-target insert summary in TSV format2: `BASE_DIR/SAMPLE_target_insert_vec_summ.tsv`
   - per long-read off-target insert summary in TSV format2: `BASE_DIR/SAMPLE_off_insert_vec_summ.tsv`

- Output for the entire experiment/run:
   - PROJECT statistics summary in TSV format: `PROJECT_ID_EXP_DESIGN_sample_stats.tsv`, including the following tab-delimited fields for each SAMPLE
      - **sample\_name**: Sample name
      - **total\_read**: Total reads (pairs)
      - **ref\_mapped**: Host-mapped reads
      - **ref\_enrich**: Host-mapped Cas9 enriched region mapped reads
      - **ref\_target**: Host-mapped ARCUS target mapped reads
      - **target\_insert**: On-target inserts
      - **target\_insert\_complete**: On-target complete inserts
      - **target\_insert\_incomplete**: On-target incomplete inserts
      - **target\_insert\_vec_mapped**: On-target inserts vector mapped
      - **target\_insert\_nuclease\_mapped**: On-target inserts nuclease-vector mapped
      - **target\_insert\_donor\_mapped**: On-target inserts donor-vector mapped
      - **target\_insert\_trans\_mapped**: On-target inserts trans-plasmid mapped
      - **target\_insert\_helper\_mapped**: On-target inserts helper-plasmid mapped
      - **target\_insert\_ref2\_mapped**: On-target inserts host2 mapped
      - **target\_insert\_vec2\_mapped**: On-target inserts vector2 mapped
      - **target\_insert\_functional\_count**: On-target inserts functional count
      - **target\_insert\_functional\_clone**: On-target inserts functional clone
      - **target\_insert\_functional\_freq**: On-target inserts functional clone frequency
      - **target\_insert\_type\_freq**: On-target inserts type frequency
      - **target\_insert\_functional\_type\_freq**: On-target inserts functional type frequency
      - **off\_insert**: Off-target inserts
      - **off\_insert\_complete**: Off-target complete inserts
      - **off\_insert\_incomplete**: Off-target incomplete inserts
      - **off\_insert\_vec\_mapped**: Off-target inserts vector mapped
      - **off\_insert\_nuclease\_mapped**: Off-target inserts nuclease-vector mapped
      - **off\_insert\_donor\_mapped**: Off-target inserts donor-vector mapped
      - **off\_insert\_trans\_mapped**: Off-target inserts trans-plasmid mapped
      - **off\_insert\_helper\_mapped**: Off-target inserts helper-plasmid mapped
      - **off\_insert\_ref2\_mapped**: Off-target inserts host2 mapped
      - **off\_insert\_vec2\_mapped**: Off-target inserts vector2 mapped
      - **off\_insert\_functional\_count**: Off-target inserts functional count
      - **off\_insert\_functional\_clone**: Off-target inserts functional clone
      - **off\_insert\_functional\_freq**: Off-target inserts functional clone frequency
      - **off\_insert\_type\_freq**: Off-target inserts type frequency
      - **off\_insert\_functional\_type\_freq**: Off-target inserts functional type frequency
      
## Notes:
- Some common options for the LSF `bsub` is explained below, try `man bsub` for more options and details:
   - `-J`: Job name (for display only)
   - `-o`: redirect job logs (stdout and stderr) to the given file 
   - `-n`: # of cores/CPUs requests for this job
   - `-M`: memory required for this job

---
