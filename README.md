# CanVar
Framework based upon ExAC to provide study specific variant frequencies. Dependencies and installation instructions are almost identical to ExAC as described in https://github.com/konradjk/exac_browser/blob/master/README.md

Installation
=======

### Getting the code

Create a directory to put all this stuff in. This will serve as the parent directory of the actual canvar_browser repository 

    mkdir canvar
    cd canvar

First (as this can run in parallel), get the datasets that the browser uses and put them into an 'exac_data' directory. All the files required can be found at:

    wget https://personal.broadinstitute.org/konradk/exac_browser/omim_info.txt.gz
    wget https://personal.broadinstitute.org/konradk/exac_browser/canonical_transcripts.txt.gz
    wget https://personal.broadinstitute.org/konradk/exac_browser/dbNSFP2.6_gene.gz
    wget https://personal.broadinstitute.org/konradk/exac_browser/gencode.gtf.gz
    wget https://personal.broadinstitute.org/konradk/exac_browser/omim_info.txt.gz
For the dbsnp file:

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz 
    gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
    tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgzsudo 
    zcat b142_SNPChrPosOnRef_105.bcp.gz | awk '$3 != ""' | perl -pi -e 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp142.txt.bgz
    tabix -s 2 -b 3 -e 3 dbsnp142.txt.bgz

Now clone the repo: 

    git clone https://github.com/konradjk/exac_browser.git

### Dependencies

You need Python and pip installed

Follow these instructions to get Python and Homebrew installed on your Mac:
http://docs.python-guide.org/en/latest/starting/install/osx/

For linux systems 

http://docs.python-guide.org/en/latest/starting/install/linux/


Install MongoDB on Mac:

    brew install mongodb
    # or
    sudo port install mongodb

Alternatively, if you are on an linux system:

https://docs.mongodb.com/v3.0/tutorial/install-mongodb-on-linux/

Create a directory to hold your mongo database files: 

    mkdir database

In a separate tab, start the mongo database server:

    mongod --dbpath database

This local server needs to be running at all times when you are working on the site.
You could do this in the background if you want or set up some startup service, or within a screeen session
if you don't want to keep an active terminal window open. 

Finally, you may want to keep the system in a virtualenv, to install:

For Mac:

    sudo port install py27-virtualenv # Or whatever version
For linux:

    sudo apt-get install virtualenv
    or sudo pip install virtualenv

You can then create a python virtual environment where the browser will live:

    mkdir canvar_env
    virtualenv canvar_env
    source canvar_env/bin/activate

Install the python requirements:

    pip install -r requirements.txt

Some packages will require Python headers (python-dev on some systems).

### Setup

The code provided here is primarily used to produce the sites VCF file and coverage file required by the framework. Other scripts have been modified to remove hardcoded references to ExAC and the populations contained within it, replacing them with arbitrary groups such as disease state. 

VCF -> sites vcf:

    python vcf_to_site_canvar.py --infile /path/to/file.vcf --phenotypes /path/to/phenotypes.csv --pops populations.csv

phenotype.csv contains phenotype data for each sample in the vcf file, it must have the columns:

    sample phenotype sex

where sample is the sample ID as stated in the header line of the VCF, phenotype is the disease (or other) state you want to define as a population and sex is specified as M or F. 

populations.csv contains all the "populations" present. These correspond to the different phenotypes. e.g. the Colorectal cancer population, it must have the columns:

    code phenotype

where code is the abbreviated form of the phenotype that will appear in the sites.vcf e.g. CRC Colorectal

To generate coverage file

    java -jar GenomeAnalysisTK.jar -T DepthOfCoverage -R human_g1k_v37.fasta -I bams.list  -L capture.bed --omitIntervalStatistics --omitLocusTable --omitPerSampleStats -nt n -o list.coverage

bams.list contains the paths to all the BAM files you want to calc coverage for. 

    perl prepare_coverage_for_tabix.pl list.coverage | bgzip -c > list.coverage.gz
    tabix -s 1 -b 2 -e 2 list.coverage.gz
    python average_coverage_calculate.py coverage_file_list capture.bed

in this example coverage_file_list would contain the path to list.coverage.gz. If you split your GATK DepthOfCoverage command in to different sample sets then each tabixed coverage file goes on a different line. A .cov file is produced with the average coverage for each position across all samples.

    bgzip -c coverage_file_list.cov > all_coverage.txt.gz
    tabix -s 1 -b 2 -e 2 all_coverage.txt.gz

all_coverage.txt.gz is then used in the canvar.py script as the coverage file. It can also be split by chromosome and added as per the original exac.py code.

Now the database is loaded from the flat files.
This is a single command, but it can take a while (can take advantage of parallel loads by modifying LOAD\_DB\_PARALLEL\_PROCESSES in exac.py):

    python manage.py load_db

You won't have to run this often - most changes won't require rebuilding the database.
That said, this is (and will remain) idempotent,
so you can run it again at any time if you think something might be wrong - it will reload the database from scratch.
You can also reload parts of the database using any of the following commands:

    python manage.py load_variants_file
    python manage.py load_dbsnp_file
    python manage.py load_base_coverage
    python manage.py load_gene_models

Then run:

    python manage.py precalculate_metrics

Then, you need to create a cache for autocomplete and large gene purposes:

    python manage.py create_cache

All genes contained within the genes_to_cache.txt file will be pre-cached. You can add as many Transcript IDs as you want, up to every transcript within the database in order to fully cache the data.

### Running the site

Note that if you are revisiting the site after a break, make sure your virtualenv is `activate`'d.

You can run the development server with:

    python canvar.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T


For testing, you can open up an interactive shell with:

    python manage.py shell



