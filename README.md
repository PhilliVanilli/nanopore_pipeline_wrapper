# nanopore_pipeline_wrapper
a wrapper for a nanopore pipeline based on the artic-ebov pipeline

## Requirements
This pipeline wrapper requires python v3.6 or higher (Pyhton 3.7 is preferred)

Several dependencies are noted in the requirements.yml file and are installed automatically when creating the conda environment.
external dependencies jvarkit and porechop are bundled with the repo but have specific install steps highlighted below.
this version of porechop is needed as it has been expanded to allow the 24 combination native borcode set to be used.

This pipeline also requires a "sample_names.csv" to indicate which barcodes correspond to which samples for the demultiplexing step.
A template of this file is included in the repo.

This pipeline assumes that the primer scheme bed files are of a specific format:
 tap seperated ".bed" file 
 
with the column headings: "genome", "start", "end", "Primer_ID" and "number"
start is the start position of the primer, relative to the reference genome, using a zero based index (pos 1 = 0)
Primer names need to include "LEFT" or "RIGHT", using "_" as a delimeter to refer to Fwd and Rev primers

If you need to use a primer scheme that is not included here, please create an issue on this github page and I will add it to the repo for you

## Step 1
Download and install the 64-bit Python 3.7 version of Miniconda/Anaconda

## Step 2 clone the pipeline repo
`git clone --recursive https://github.com/ColinAnthony/nanopore_pipeline_wrapper.git`

 change into the repo directory
 
 `cd nanopore_pipeline_wrapper`

## Step 3 create the conda environment
`conda env create -f requirements.yml`

# Step 3.1 activate the conda env
activate the environment

`source activate nanop`

## Step 3.2
install jvarkit (http://lindenb.github.io/jvarkit/SAM4WebLogo.html, https://github.com/lindenb/jvarkit.git)

This tool converts a bam/sam file to a multi sequence alignment

To install run:

cd into the jvarkit folder in this repo

` cd jvarkit`

and run

`./gradlew sam4weblogo`

If this fails

install java jkd 11

`sudo apt-get install openjdk-11-jdk`

and repeat the step above

change directory back to the main repo

`cd ..`

##  Step 3.3
install porechop 

change into the Porechop directory

`cd Porechop`

install porechop

`python3 setup.py install`

# Running the pipeline:

## activate the environment

`source activate nanop`

## make QCplots or raw data
`python step_1_plot_sequencing_quality.py -in <path to where the ouput folder will be made> -s <the path and name of the sequencing_summary.txt file>`

## process the basecalled fastq files to sample consensus sequences

`python step_2_process_basecalled_reads.py -in <path to fastq/fast5 folders> -s < path and name of sample_names.csv file>`
