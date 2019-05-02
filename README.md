Prerequisite:
1.Python 2.7 with the following packages
	a.Pandas
	b.Matplotlib
2.Shapeit (only if you want to use our package for phasing your variants). This can be downloaded here: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download
3.Make sure to change the path to different tools in the config file (config_file/config.txt)

Included in the package:
1.Samtools 1.9
2.1000 genome SNPs 

This package consists of three scripts which should be run sequentially:

Script 1:
Input: 
This script takes a valid bulk vcf as input and does the following:
1.Subsets the variants in the vcf to only heterozygous SNP present in 1000 genome.
2.Takes the heterozygous 1000 genome SNPs and runs a phasing algorithm called Shapeit on it.

Output: 
The output consists of two vcfs:
1.VCF with heterozygous snps from 1000 genome.
2.Phased VCF with heterozygous snps from 1000 genome (this is the vcf you will be using in the second and third scripts).

Usage:
python First_process_vcf.py -h
usage: First_process_vcf.py [-h] -V VCF_FILE -O OUTPUT_DIR -S SAMPLE_NAME
                            [-l CHROMOSOME_NUMBER] [-f]

Preprocessing of vcf file

optional arguments:
  -h, --help            show this help message and exit
  -V VCF_FILE, --VCF_file VCF_FILE
                        Path to VCF file
  -O OUTPUT_DIR, --Output_dir OUTPUT_DIR
                        Path to directory where results will be written
  -S SAMPLE_NAME, --Sample_name SAMPLE_NAME
                        Name of the sample
  -l CHROMOSOME_NUMBER, --chromosome_number CHROMOSOME_NUMBER
                        chromosome number
  -f, --fast_option     run all chromosome in parallel.Uses more memory

Example:
python First_process_vcf.py -V test.vcf -O test_out_directory -S test_sample_name -f

Notes:
1.It is recommended that the vcf is run through GATK VQSR before this step. It is not necessary.
2.By using the "-l"option you can run each chromosome separately
3.By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.
4.This script generates files per chromosome and you will have to merge them to one file before proceeding to the next step.

Script 2:
Input: 
This script takes a valid vcf (either from script 1 or otherwise) and a valid bam as input.
1.It uses the variants from the vcf and generates a pileup using the bam files. By default it filters any read with less than 20 mapping quality and any base with less than 20 base quality by default. But this can be changed using the parameters

Output: 
It outputs a file with allele frequency of each base at each position in the vcf provided. This file will be used in the third script to calculate haplotype allele frequency.

Usage:
python Second_calculate_allele_frequency.py  -h
usage: Second_calculate_allele_frequency.py [-h] -B BAM_FILE -V VCF_FILE -O
                                            OUTPUT_DIR -S SAMPLE_NAME
                                            [-q MIN_BASE_QUALITY]
                                            [-Q MIN_MAPPING_QUALITY]
                                            [-l CHROMOSOME_NUMBER] [-a] [-n]

Preprocessing of vcf file

optional arguments:
  -h, --help            show this help message and exit
  -B BAM_FILE, --BAM_file BAM_FILE
                        Path to bam file
  -V VCF_FILE, --VCF_file VCF_FILE
                        Path to VCF file
  -O OUTPUT_DIR, --Output_dir OUTPUT_DIR
                        Path to directory where results will be written
  -S SAMPLE_NAME, --Sample_name SAMPLE_NAME
                        Name of the sample
  -q MIN_BASE_QUALITY, --min_base_quality MIN_BASE_QUALITY
                        Mininum base quality
  -Q MIN_MAPPING_QUALITY, --min_mapping_quality MIN_MAPPING_QUALITY
                        Mininum mapping quality
  -l CHROMOSOME_NUMBER, --chromosome_number CHROMOSOME_NUMBER
                        chromosome number
  -f, --all_chromosome_in_parallel
                        run all chromosome in parallel
  -n, --no_chr_splitting
                        Doesn't split by chromosome

Example:

python allele_frequency.py -B test.bam -V test.1000G_het_snps.vcf -O test_directory -S test_sample_name -f

Notes:
1.You can skip the first step and start with the second step if you already have a phase vcf, but keep the following in mind:
	a.Make sure to subset the variants to only 1000 genome SNPS
	b.Make sure to include only heterozygous SNPs
2.By using the "-l"  option you can run each chromosome separately
3.By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.
4.This script generates also gives you the option to not split by chromosome and run everything as a whole with the “-n” option. This can be very time consuming if your data is too big.
5.If you are not using the "-n" option, you will have to merge them to one file before proceeding to the next step.

Script 3:
Input: 
This script takes a the following files:
1.Output file generate by the script 2 (AF_FILE)
2.Phased vcf either generated by script 1 (GERM_HAP), or provided by the uses. If the phased vcf is generated outside of this package, please make sure:
	a.The vcf had heterozygous 1000 genome SNPs only
	b.The genotype (GT field) is in the right phased format

Output: 
The output consists of an allele frequency plot which shows the quality of amplification for the single cell samples.

python Third_generate_allele_frequency_plot.py -h
usage: Third_generate_allele_frequency_plot.py [-h] -a AF_FILE -g GERM_HAP -o
                                               OUTPUT_DIR -S SAMPLE_NAME
                                               [-n SNPS]

optional arguments:
  -h, --help            show this help message and exit
  -a AF_FILE, --AF_file AF_FILE
                        path sample allele freq file
  -g GERM_HAP, --Germ_hap GERM_HAP
                        Path to germline hap file
  -o OUTPUT_DIR, --Output_dir OUTPUT_DIR
                        Path to output file
  -S SAMPLE_NAME, --Sample_name SAMPLE_NAME
                        Name of the sample
  -n SNPS, --Snps SNPS  Number_of_snps

Example:
python Third_generate_allele_frequency_plot.py -a test_AF.txt -g Test_germ_hap.vcf -o test_directory/ -S test_sample_name

Notes:
1.The "-n" option is for the number of SNP unit used. The default it 100. This means that the allele frequency is calculated over 100 SNPs and represented as a unit in the plot. This number is based on a whole genome sequenced sample with 5 million reads. But please change this unit appropriately for your samples to assure proper comparisons across samples with different read counts.

