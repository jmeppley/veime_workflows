
# Workflows for making VEIMEs

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

## Installation

Note: This workflow wil not work on windows or macos. The required dependencies are only available for linux at the moment.

First, clone this repository locally.

Second, install the necessary software with conda:

    conda env create -p ./conda.env -f ./conda.yaml

Note: replace "conda" with "mamba" (if you have mamba installed) for faster installation.

At this point, you can test the installation with the reduced set of data in ./test/data
without downlading any database.

## Test Run
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

Optionally, cluster polished monomers alongside known satellites:

    snakemake -s paper/Snakefile.module -j 10 -p

Optionally, build a tree of tyrosine integrases:

    snakemake -s paper/Snakefile.iqtree -j 10 -p

Note, iqtree will fail if you request more threads (via -j) than are available on your
system.

## Download Databases

Finally, to run on real data, you'll need to download the VOG and EggNOG
databases:

### VOGDB
  
To get the latest VOG definitions, create a local folder to store them in,
[download](https://vogdb.org/download) and upack the archive, and concatenate all the hmm files into one DB:

   mkdir -p /local/path/to/VOGS
   wget -c http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
   tar -zxvf vog.hmm.tar.gz
   cat VOG[0-9]*.hmm > VOGS.hmm

You will also want to download the [VOG
descriptions](http://fileshare.csb.univie.ac.at/vog/latest/vog.annotations.tsv.gz) file. This contains the
functional annotations for each VOG family.

### EggNOG

Eggnog mapper, installed in the conda environment, comes with a script to
download the database for you. All you have to do is tell snakemake where to
store the data. The first time, snakemake will download the data. On subsequent
runs, if you specify the same db location, the previously downlaoded data will
be reused.

WARNING: The current EggNOG database will take up about 50GB of space and needs
an extra 20GB or so to download and uncompress files. 

If you do not speciy an eggnog DB location, only HMM annotations are used.

If you do specify an eggnog datanase, you must tell snakemake to use conda. The 
eggnog environment is incompatible with other items in the main environment.

### Specifying databse locations

When running the Snakemake.module or snakemake.iqtree workflows, specify
database locations as follows:

    snakemake -j 10 -p -s paper/Snakefile.module \
        --config hmm_db=/local/path/to/VOGS/VOGS.hmm \
                 eggnog_data_dir=/local/path/to/eggnog_data \
        --use-conda --conda-frontend=mamba

## Configuring the Workflows
As with the database locations (see above), configuration parameters can be
specified on the command line with the "--config key=value" pattern. User
modifiable params are found at the top of each Snakefile. Other than the
database locations, the most useful ones are:

 * working_dir: (all workflows) this is where the output will be written
 * reads_fasta: (Snakefile.frags) the initial set of concatemeric reads
 * cmer_data: (Snakefile.frags) table with repeat sizes
 * short_reads: (Snakefile.racon) fastq with short reads for polishing
 * hmm_db: (Snakefile.module, Snakefile.iqtree) hmm file with VOG models
 * eggnog_data_dir: (Snakefile.module) location to find (or put) eggnog
 * eggnog_tmp_dir: (Snakefile.module) fast scratch location for tmp files
 * tree_annot: (Snakefile.iqtree) gene or genes (HMM ids) to build trees on

## Troubleshooting

### MacOS is not supported

Sorry. As of now, some key dependencies are unavailable for MacOS. 

### Snakemake errors: Check the log file first

In general, if there is an error, snakemake will report something like this:

    Error in rule align_tree_dir_vog:
    jobid: 4
    output: test/run/tree_VOG00035/VOG00035.aln
    shell:
        hmmalign -o test/run/tree_VOG00035/VOG00035.aln --amino --informat FASTA --trim  test/run/VOG_tmp/VOG00035.hmm test/run/tree_VOG00035/VOG00035.faa             > test/run/tree_VOG00035/VOG00035.aln.log 2>&1
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Most rules are configured to capture all output to rule-specific log files. The first thing to do is to look at the log file. In the case above:

    $ cat test/run/tree_VOG00035/VOG00035.aln.log
    
    	Error: Sequence file test/run/tree_VOG00035/VOG00035.faa is empty or misformatted

In this case, there were no genes found with VOG00035. 

### Snakefile.iqtree: empty faa file
If you encounter the error above, where the faa file for tree bulding is empty, check that the VOG is in your annotation table:

    $ grep -c VOG00035 test/run/genomes.genes.tsv
    0

This can happen if you use default gene (VOG00035) for tree building or updated your VOG HMMs. The VOG ids can change between versions, and genes that were present before may have new names.

