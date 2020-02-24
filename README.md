# **Prerequisite:**
1. Python 3.6 with the following packages
	1. Pandas
	1. Matplotlib
	1. Scipy
	1. Numpy
1. Shapeit (only if you want to use our package for phasing your variants). This can be downloaded here: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download. Make sure to change the path to different tools in the config file (config_file/config.txt)
1. Samtools
# **Configuration setup:** 
1. Reference
1. Dbsnp
1. Shapeit

# **This package consists of three scripts which should be run sequentially:**

## **Script 1:**
### Input:
This script takes a valid bulk vcf as input and does the following:
1. Subsets the variants in the vcf to only heterozygous SNP present in 1000 genome.
1. Takes the heterozygous 1000 genome SNPs and runs a phasing algorithm called Shapeit on it.

### Output:
The output consists of two vcfs:
1. VCF with heterozygous snps from 1000 genome.
1. Phased VCF with heterozygous snps from 1000 genome (this is the vcf you will be using in the second and third scripts).
### Usage:
```javascript
python First_process_vcf.py -h
usage: First_process_vcf.py [-h] -v VCF_FILE -o OUTPUT_DIR -s SAMPLE_NAME
                            [-l CHROMOSOME_NUMBER] [-f]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF_FILE, --VCF_file VCF_FILE
                        Path to VCF file
  -o OUTPUT_DIR, --Output_dir OUTPUT_DIR
                        Path to directory where results will be written
  -s SAMPLE_NAME, --Sample_name SAMPLE_NAME
                        Name of the sample
  -l CHROMOSOME_NUMBER, --chromosome_number CHROMOSOME_NUMBER
                        chromosome number
  -f, --fast_option     run all chromosome in parallel.Uses more memory
```
Example:

```python First_process_vcf.py -v test.vcf -o test_out_directory -s test_sample_name -f```

#### Notes:
1. It is recommended that the vcf is run through GATK VQSR before this step. It is not necessary.
1. By using the "-l"option you can run each chromosome separately.
1. By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.

## **Script 2:**
### Input: 
This script takes a valid vcf (either from script 1 or otherwise) and a valid bam as input.
1.It uses the variants from the vcf and generates a pileup using the bam files. By default it filters any read with less than 20 mapping quality and any base with less than 20 base quality by default. But this can be changed using the parameters

### Output: 
It outputs a file with allele frequency of each base at each position in the vcf provided. This file will be used in the third script to calculate haplotype allele frequency.This script also provide a recommendation for the number of SNPs to use as a SNP unit in the next step. 
### Usage:
```javascript
python Second_calculate_allele_frequency.py  -h
usage: Second_calculate_allele_frequency.py [-h] -b BAM_FILE -v VCF_FILE -o
                                            OUTPUT_DIR -s SAMPLE_NAME
                                            [-q MIN_BASE_QUALITY]
                                            [-Q MIN_MAPPING_QUALITY]
                                            [-l CHROMOSOME_NUMBER] [-a] [-n]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_FILE, --BAM_file BAM_FILE
                        Path to bam file
  -v VCF_FILE, --VCF_file VCF_FILE
                        Path to VCF file
  -o OUTPUT_DIR, --Output_dir OUTPUT_DIR
                        Path to directory where results will be written
  -s SAMPLE_NAME, --Sample_name SAMPLE_NAME
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
```
Example:

```python Second_calculate_allele_frequency.py -b test.bam -v test.1000G_het_snps.vcf -o test_directory -s test_sample_name -f```

#### Notes:
1. You can skip the first step and start with the second step if you already have a phase vcf, but keep the following in mind:
	a.Make sure to subset the variants to only 1000 genome SNPS
	b.Make sure to include only heterozygous SNPs
1. By using the "-l"  option you can run each chromosome separately. You will have to merge the per chromosome files prior to the third step.
1. By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.
1. This script generates also gives you the option to not split by chromosome and run everything as a whole with the “-n” option. This can be very time consuming if your data is too big.

## **Script 3:**
### Input: 
This script takes a the following files:
1. Output file generate by the script 2 (AF_FILE)
1. Phased vcf either generated by script 1 (GERM_HAP), or provided by the uses. If the phased vcf is generated outside of this package, please make sure:
	1. The vcf had heterozygous 1000 genome SNPs only
	1. The genotype (GT field) is in the right phased format

### Output: 
The output consists of an allele frequency plot which shows the quality of amplification for the single cell samples. It also generates a file with standard deviation and allele dropout values
### Usage
```javascript
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
```
Example:
```python Third_generate_allele_frequency_plot.py -a test_AF.txt -g Test_germ_hap.vcf -o test_directory/ -S test_sample_name```

#### Notes:
1. The "-n" option is for the number of SNP unit used. The default it 100. This means that the allele frequency is calculated over 100 SNPs and represented as a unit in the plot. This number is based on a whole genome sequenced sample with 5 million reads. But please change this unit appropriately for your samples to assure proper comparisons across samples with different read counts. The second script gives a recommendation on the number of SNPs to use as a SNP unit.

