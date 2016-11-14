# CanVar
Framework based upon ExAC to provide study specific variant frequencies. Dependencies and installation instructions are the same as ExAC and are described in https://github.com/konradjk/exac_browser/blob/master/README.md

The code provided here is primarily used to produce the sites VCF file and coverage file required by the framework. Other scripts have been modified to remove hardcoded references to ExAC and the populations contained within it, replacing them with arbitrary groups such as disease state. 

VCF -> sites vcf

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

all_coverage.txt.gz is then used in the exac.py script as the coverage file. It can also be split by chromosome and added as per the original exac.py code.





