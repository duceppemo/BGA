## Info
This pipeline is designed to assembly bacterial genomes from long reads. Although it could technically work for both Nanopore and Pacbio,
it was design with Nanopore support in mind. For example, only Nanopore long read pre-processing (trimming and filtering) is supported (optional).
If Illumina paired-end sequencing data is available, this pipeline will polish the long read assembly using the short reads (optional). Short read assembly only is not supported.

This pipeline takes fastq files as input and contains up to 5 steps:
1- Long read adapter removal using Porchop (optional).
2- Long read quality filtering using Filtlong (optional).
3- Long read assembly using Flye or Shasta.
4- Short read trimming using Fastp (optional).
5- Long read assembly polishing using short reads (optional).

If you need to perform basecalling to produce fastq files from Nanopore fast5 files, please refer this sister pipeline [here](https://github.com/duceppemo/basecall_nanopore).

## Important notes
* The default assembler Flye expect high quality Nanopore data. Please use fastq basecalled with the most recent Guppy version using the Super Accuracy model for better results.
* To be able to perform the hybrid assembly, long reads and short reads of the same sample must have the same name. Everything before the first "_" will be considered the sample name. If not "_" are present, then it's everything before the first ".".
* The Illumina paired-end file names must contain "_R1" (forward reads) or "_R2" (reverse reads) to be recognized by the pipeline.

## Installation
```bash
# Create virtual environment
mamba create -n BGA -c conda-forge -c bioconda -c plotly -y python=3.10.11 nextpolish=1.4.1 bwa=0.7.17 samtools=1.17 \
    porechop=0.2.4 filtlong=0.2.1 minimap2=2.26 flye=2.9.2 shasta=0.11.1 qualimap=2.2.2d bbmap=39.01 bandage=0.8.1 \
    fastp=0.22.0 ntedit=1.3.5 polypolish=0.5.0 pandas=1.5.3 seqtk=1.4 quast=5.2.0 medaka=1.8.0 mummer4=4.0.0rc1 \
    gnuplot=5.4.5 plotly=5.15.0

# Activate virtual environment
conda activate BGA

# Install tool to you favorite location
cd ~/prog
git clone https://github.com/duceppemo/BGA

# Test tool
cd BGA
python bga.py -h
```
## Usage
```commandline
usage: python bga.py [-h] -l /path/to/nanopore_folder [-s /path/to/illumina_folder] -o /path/to/output_folder/ [-a {flye,shasta}]
                     [--min-size 3000] [--size 5000000] [--trim-long] [--filter-long] [--polish] [--trim-short] [-t 16] [-p 2]
                     [-m 57] [-v]

Bacterial Genome assembly.

options:
  -h, --help            show this help message and exit
  -l /path/to/nanopore_folder, --long-reads /path/to/nanopore_folder
                        Folder that contains the fastq files from Nanopore. Mandatory.
  -s /path/to/illumina_folder, --short-reads /path/to/illumina_folder
                        Folder that contains the fastq files from Illumina. Illumina file names must match Nanopore files
                        names (everything before the first underscore ("_"). Optional.
  -o /path/to/output_folder/, --output /path/to/output_folder/
                        Folder to hold the result files. Mandatory.
  -a {flye,shasta}, --assembler {flye,shasta}
                        Assembly method. Default "flye". Optional.
  --min-size 3000       Minimum read size for Shasta assembler or minimum read overlap for Flye. Default 3000 for Shasta and
                        auto for Flye. Optional.
  --size 5000000        Override automatically detected reference size for Flye. If entered manually, it will also be used
                        to during the long read filtering step (to retain the 100x to reads, if coverage allows. Optional.
  --trim-long           Trim long reads with Porechop prior assembly. Default is False.
  --filter-long         Filter long reads with Filtlong prior assembly. Drop bottom 5%. Default is False.
  --polish              Polish long read assembly with Illumina paired-end data. Default is False.
  --trim-short          Trim paired-end Illumina reads with FastP prior polishing. Default is False.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional.
  -p 2, --parallel 2    Number of samples to process in parallel. Default is 2. Optional.
  -m 57, --memory 57    Memory in GB. Default is 85% of total memory (57)
  -v, --version         show program's version number and exit
```
## Examples
Here we perform a hybrid assembly using both Nanopore and Illumina:  
```bash
python bga.py \
    -l /data/nanopore \
    -s /data/illumina \
    -o /analyses/hybrid_assemblies \
    --trim-long \
    --filter-long \
    --trim-short \
    --polish \
    -t 16 \
    -p 3
```
## TODO
* Maybe add support for single-end Illumina data?
* Process multiple sample in parallel for the polishing step.
* Maybe add more options to control the long and short read pre-processing steps.
