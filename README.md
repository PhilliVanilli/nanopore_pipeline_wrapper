# nanopore_pipeline_wrapper
a wrapper for a nanopore pipeline based on the artic-ebov pipeline

## Step 1
Download and install the 64-bit Python 3.7 version of Miniconda/Anaconda

## Step 2 clone the pipeline repo
`git clone https://github.com/ColinAnthony/nanopore_pipeline_wrapper.git`

 change into the repo directory
 
 `cd nanopore_pipeline_wrapper`

## Step 3 create the conda environment
`conda env create -f environment.yml`

# Step 4 activate the conda env
`source activate nanop`

# Running the pipeline:

## make QCplots
`python step_1_plot_sequencing_quality.py -in <path to where the ouput folder will be made> -s <the path and name of the sequencing_summary.txt file>`

## process the basecalled fastq files to sample consensus sequences

`python step_2_process_basecalled_reads.py -in <path to fastq/fast5 folders> -s < path and name of sample_names.csv file>`