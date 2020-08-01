## **Prerequisite:**
1. Python 3.6 with the following packages
	1. Pandas
	1. Matplotlib
	1. Scipy
	1. Numpy
1. Other dependencies:
    1. Samtools
    1. Shapeit
## **Setup:** 
### Download
  ```
  git clone https://github.com/abyzovlab/Scellector.git    
  ```
### Configuration setup
The method required the following tools and references. Please follow the steps to download the dependencies. 
1. Shapeit (download and update "SHAPEIT_PATH"  in config_file/config.txt):

   ```
   wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz
   
   tar zvfx shapeit.v2.r904.glibcv2.12.linux.tar.gz
   ```

1. Shapeit reference (download and update "SHAPIT_REF"  in config_file/config.txt):
    
   ```
   wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz

   tar -xvzf 1000GP_Phase3.tgz 
   ```
    
1. Samtools (download and update "SAMTOOLS"  in config_file/config.txt):

    ```
   wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
    
   tar -xvf samtools-1.9.tar.bz2
    
   cd samtools-1.9/
    
    /.configure
    
    make
   ```

1. Reference (download and update "REFERENCE"  in config_file/config.txt):

    ```
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
    
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
    
    gunzip human_g1k_v37_decoy.fasta.gz
    
    gunzip human_g1k_v37_decoy.fasta.fai.gz
   ```

1. 1000 genome SNPs (download and update "G1K_SNP"  in config_file/config.txt):

    ```
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
   ```

## This package consists of three scripts which can be run sequencially:
### Script-1:
#### Input:
This script takes a valid bulk vcf as input and does the following:
1. Subsets the variants in the vcf to only heterozygous SNP present in 1000 genome.
1. Takes the heterozygous 1000 genome SNPs and runs a phasing algorithm called Shapeit on it.

#### Output:
The output consists of a vcf:
1. Phased VCF with heterozygous snps from 1000 genome (this is the vcf you will be using in the second and third scripts).
#### Usage:
```
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

```
python First_process_vcf.py -v example_data/example_input/B01_bulk_downsampled.vcf -o example_data/example_output/script_1/ -s B01_cell_1 -f
```

##### Notes:
1. It is recommended that the vcf is run through GATK VQSR before this step. It is not necessary.
1. By using the "-l" option you can run each chromosome separately.
1. By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.

### Script-2:
#### Input: 
This script takes a valid vcf (either from script 1 or otherwise) and a valid bam as input.
1.It uses the variants from the vcf and generates a pileup using the bam files. By default it filters any read with less than 20 mapping quality and any base with less than 20 base quality. But this can be changed using the parameters

#### Output: 
It outputs a file with allele frequency of each base at each position in the vcf provided. This file will be used in the third script to calculate haplotype allele frequency.This script also provide a recommendation for the number of SNPs to use as a SNP unit in the next step. 
#### Usage:
```
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

```
python Second_calculate_allele_frequency.py -b example_data/example_input/B01_cell_1_downsampled.bam -v example_data/example_output/script_1/B01_cell_1.vcf -o example_data/example_output/script_2/ -s B01_cell_1 -f
```

##### Notes:
1. You can skip the first step and start with the second step if you already have a phase vcf, but keep the following in mind:
	1. Make sure to subset the variants to only 1000 genome SNPS
	1. Make sure to include only heterozygous SNPs
1. By using the "-l"  option you can run each chromosome separately. You will have to merge the per chromosome files prior to the third step.
1. By using "-f" option you will be running all chromosome in parallel. This is memory intensive but can be really fast depending on your computing capacity.
1. This script generates also gives you the option to not split by chromosome and run everything as a whole with the “-n” option. This can be very time consuming if your data is too big.

### Script-3:
#### Input: 
This script takes a the following files:
1. Output file generate by the script 2 (AF_FILE)
1. Phased vcf either generated by script 1 (GERM_HAP), or provided by the uses. If the phased vcf is generated outside of this package, please make sure:
    1. The vcf had heterozygous 1000 genome SNPs only
    1. The genotype (GT field) is in the right phased format

#### Output: 
The output consists of an allele frequency plot which shows the quality of amplification for the single cell samples. It also generates a file with standard deviation and allele dropout values
#### Usage
```
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
```
python Third_generate_allele_frequency_plot.py -a example_data/example_output/script_2/B01_cell_1.AF.txt -g example_data/example_output/script_1/B01_cell_1.vcf -o example_data/example_output/script_3/ -s B01_cell_1 -n 200
```

##### Notes:
1. The "-n" option is for the number of SNP unit used. The default it 100. This means that the allele frequency is calculated over 100 SNPs and represented as a unit in the plot. This number is based on a whole genome sequenced sample with 5 million reads. But please change this unit appropriately for your samples to assure proper comparisons across samples with different read counts. The second script gives a recommendation on the number of SNPs to use as a SNP unit.

### Example input and output files:

```
# Script-1 : Takes a vcf, preprocesses it and generates phased data using shapeit
python First_process_vcf.py -v example_data/example_input/B01_bulk_downsampled.vcf -o example_data/example_output/script_1/ -s B01_cell_1 -f

# Script-2 : Generates allele frequency of variants using Bam file and vcf files
python Second_calculate_allele_frequency.py -b example_data/example_input/B01_cell_1_downsampled.bam -v example_data/example_output/script_1/B01_cell_1.vcf -o example_data/example_output/script_2/ -s B01_cell_1 -f

# Script-3 : Combines the data from phasing file and allele frequency file from first and second script to generate VAF plot
python Third_generate_allele_frequency_plot.py -a example_data/example_output/script_2/B01_cell_1.AF.txt -g example_data/example_output/script_1/B01_cell_1.vcf -o example_data/example_output/script_3/ -s B01_cell_1 -n 200
```
1. Input examples can be found in "example_data/example_input":
    1. B01_cell_1_downsampled.bam : This is a downsampled single cell bam file that can be used to test the tool.
    1. B01_bulk_downsampled.vcf : This is a downsampled bulk vcf file that can be used to test the tool.
    
1. Output of the tool can be found in "example_data/example_output":
    1. Script-1: results from running First_process_vcf.py
    1. Script-2: results from running Second_calculate_allele_frequency.py
    1. Script-3: results from running Third_generate_allele_frequency_plot.py
   
##### Notes:
 Please keep in mind that these examples are only to test the execution of the tool and is not a real example. You can find real examples in our paper. 

