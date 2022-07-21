
= Workflows for making VEIMEs

This repo contains a collection of makefiles that will take concatemeric long
reads, generate polished monomers, cluster the results alongside erference
sequences, and build gene trees.

NOTE: currently, it is assumed that you have already used the
concatemer_finding workflow to isolate concatemeric reads. The starting point
is:

 * a fastq file of repetetive long reads
 * a tsv file indexed on read ids with a "repeat_size" column

You will also need to download the VOG HMM db and eggnog DB for the last two
steps (clustering and tree building).

== Installation

First, clone this repository locally.

Second, install the necessary software with conda:

    conda env create -p ./conda.env -f ./conda.yaml

Note: replace "conda" with "mamba" for faster installation.

At this point, you can test the installation with the reduced set of data in ./test/data
without downlading any database.

Finally, to run on real data, you'll need to download the VOG and EggNOG
databases:

### TODO ###

== Test Run
The following commands will run through the workflows. Note, that for the
snakemake commandas, you can alther the number of threads used with the "-j"
option.

First, activate the conda environment:

    conda activate ./conda.env

Next, extract monomers:

    snakemake -s paper/Snakefile.frags -j 10 -p

Then, self-polish the monomers:

    snakemake -s paper/Snakefile.veime -j 10 -p

Polish again, but with short reads:

    snakemake -s paper/Snakefile.racon -j 10 -p

Cluster polished monomers alongside known satellites:

    snakemake -s paper/Snakefile.module -j 10 -p

Build a tree of tyrosine integrases:

    snakemake -s paper/Snakefile.iqtree -j 10 -p

