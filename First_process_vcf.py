"""Takes a vcf, preprocesses it and generates phased data using shapeit"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Development"

import os
import argparse
import sys
import Util
import gzip
import threading
import time


def argument_parse(script_path):
    '''Parses the command line arguments'''
    parser = argparse.ArgumentParser(description='Preprocessing of vcf file')
    parser.add_argument("-V", "--VCF_file", help="Path to VCF file", required=True, type=Util.FileValidator)
    parser.add_argument("-O", "--Output_dir", help="Path to directory where results will be written", required=True)
    parser.add_argument("-S", "--Sample_name", help="Name of the sample", required=True)
    parser.add_argument("-l", "--chromosome_number", help="chromosome number ")
    parser.add_argument("-f", "--fast_option", help="run all chromosome in parallel.Uses more memory ",
                        action='store_true')
    return parser


def shapeit(vcf_file, output_dir, sample_name, chromosome_number):
    '''running shapeit on vcf file'''
    MAP = os.path.join(config["SHAPIT_REF"], "genetic_map_chr" + chromosome_number + "_combined_b37.txt")
    HAP = os.path.join(config["SHAPIT_REF"], "1000GP_Phase3_chr" + chromosome_number + ".hap.gz")
    LEG = os.path.join(config["SHAPIT_REF"], "1000GP_Phase3_chr" + chromosome_number + ".minimal.legend.gz")
    SAM = os.path.join(config["SHAPIT_REF"], "1000GP_Phase3.sample")
    CHK = os.path.join(output_dir, sample_name + "." + chromosome_number + ".alignments")
    OUT = output_dir + "/" + sample_name + ".chr" + chromosome_number
    OUT_VCF = output_dir + "/" + sample_name + ".chr" + chromosome_number + ".vcf"
    OUT_VCF_fh = open(OUT_VCF, 'w')
    shapeit = os.path.join(config["SHAPEIT_PATH"], "shapeit")
    vcf = vcf_file
    cmd = ' '.join([shapeit, "check",
                    "--input-vcf", vcf,
                    "--input-map", MAP,
                    "--input-ref", HAP, LEG, SAM,
                    "--output-log", CHK,
                    "--output-max", OUT])
    print(cmd)
    os.system(cmd)
    if os.path.exists(CHK + ".snp.strand.exclude") and os.path.getsize(CHK + ".snp.strand.exclude") > 0:
        cmd = ' '.join([shapeit, "--thread 4",
                        "--input-vcf", vcf,
                        "--input-map", MAP,
                        "--exclude-snp", CHK + ".snp.strand.exclude",
                        "--input-ref", HAP, LEG, SAM,
                        "--output-log", OUT,
                        "--output-max", OUT])
    else:
        cmd = ' '.join([shapeit, "--thread 4",
                        "--input-vcf", vcf,
                        "--input-map", MAP,
                        "--input-ref", HAP, LEG, SAM,
                        "--output-log", OUT,
                        "--output-max", OUT])
    print(cmd)
    os.system(cmd)
    OUT_VCF_fh.write("##fileformat=VCFv4.1\n")
    OUT_VCF_fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")
    for i in open(OUT + ".haps"):
        line = i.strip().split(" ")
        line_out = "\t".join([line[0], line[2], ".", line[3], line[4], ".", ".", ".", "GT", line[5] + "|" + line[6]])
        OUT_VCF_fh.write(line_out + "\n")


def read_1000_genome(vcf, chromosome_number=None):
    '''reads the 1000 genome file and returns a dict of positions'''
    if chromosome_number == None:
        print("---------- reading all SNPs from 1000 genome file----------")
    else:
        print("---------- reading chr" + chromosome_number + " SNPs from 1000 genome file----------")
    g1k_snp = config["G1K_SNP"]
    if g1k_snp.endswith(".gz"):
        g1k_snp_fh = gzip.open(g1k_snp)
    else:
        g1k_snp_fh = open(g1k_snp)
    g1k_chr_status = Util.check_vcf_with_chr_or_not(g1k_snp)
    vcf_chr_status = Util.check_vcf_with_chr_or_not(vcf)
    chromosome_number = chromosome_number
    if g1k_chr_status != "with_chr" and vcf_chr_status == "with_chr" and chromosome_number != None:
        g1k_chromosome_number = "chr" + chromosome_number
    else:
        g1k_chromosome_number = chromosome_number
    g1k_snp_dict = {}
    for i in g1k_snp_fh:
        line = i.strip().split("\t")
        if i.startswith("#"):
            continue
        if vcf_chr_status != g1k_chr_status and g1k_chr_status == "without_chr":
            line[0] = "chr" + line[0]
        if chromosome_number != None and line[0] != g1k_chromosome_number:
            continue
        ref = line[2]
        alt = line[3]
        if len(ref) > 1 or len(alt) > 1:
            continue
        snp_id = line[2]
        pos = "_".join([line[0], line[1], line[2], line[3]])
        g1k_snp_dict[pos] = snp_id
    g1k_snp_fh.close()
    return g1k_snp_dict


def extract_heterozygous_snps(vcf, vcf_out, chromosome_number, g1k_snp_dict):
    '''extracts the SNPs present in 1000G and keeps only mono allelic heterozygous snps'''
    g1k_snp = config["G1K_SNP"]
    vcf_out_fh = open(vcf_out, 'w')

    if g1k_snp.endswith(".gz"):
        g1k_snp_fh = gzip.open(g1k_snp)
    else:
        g1k_snp_fh = open(g1k_snp)
    print("filtering vcf file......")
    if vcf.endswith(".gz"):
        vcf_fh = gzip.open(vcf)
    else:
        vcf_fh = open(vcf)

    vcf_chr_status = Util.check_vcf_with_chr_or_not(vcf)
    if vcf_chr_status == "without_chr":
        vcf_chromosome_number = chromosome_number
    else:
        vcf_chromosome_number = "chr" + chromosome_number
    print(vcf_chromosome_number)
    pos_dict = {}
    for i in vcf_fh:
        line = i.strip().split("\t")
        if i.startswith("#"):
            vcf_out_fh.write(i)
            continue
        if line[0] != vcf_chromosome_number:
            continue
        filter_col = line[6]
        GT = ""
        for n, j in enumerate(line[8].split(":")):
            if j == "GT":
                GT = line[9].split(":")[n]
        pos = "_".join([line[0], line[1], line[3], line[4]])
        if pos in pos_dict:
            continue
        else:
            pos_dict[pos] = 1
        if pos not in g1k_snp_dict or filter_col != "PASS" or GT != "0/1":
            continue
        new_i = "\t".join([line[0], line[1], line[2], line[3], line[4], ".", ".", ".", "GT", GT])
        vcf_out_fh.write(new_i + "\n")
    vcf_out_fh.close()


def running_preprocessing(vcf, output_dir, sample_name, chromosome_number, g1k_snp_dict):
    '''given a chromosome number this runs pre processing step'''
    print("-------------- processing for chr" + chromosome_number + "-------------------------")
    vcf_snp_het = os.path.join(output_dir, sample_name + ".chr" + chromosome_number + ".het_snps.vcf")
    print("1. Extracting heterozygous allelic SNPS present in 1000 Genomes")
    extract_heterozygous_snps(vcf, vcf_snp_het, chromosome_number, g1k_snp_dict)
    print("---------------Done extracting SNPS-----------------------")
    print("2. Generating phased VCF")
    shapeit(vcf_snp_het, output_dir, sample_name, chromosome_number)
    print("---------------Done phasing vcf---------------------------")


def main():
    script_name = os.path.basename(sys.argv[0])
    script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    parser = argument_parse(script_name)
    try:
        arg = parser.parse_args()
    except:
        if "JOB_ID" not in os.environ:
            print(sys.exc_info()[1])
            exit()
        else:
            print(sys.exc_info()[1])
            exit(100)
    global config
    config_file = os.path.join(script_path, "config_file/config.txt")
    try:
        config = Util.ParseConfig(config_file)
    except:
        if "JOB_ID" not in os.environ:
            print(sys.exc_info()[1])
            exit()
        else:
            print(sys.exc_info()[1])
            exit()
    ### assigning values to variable
    output_dir = arg.Output_dir
    sample_name = arg.Sample_name
    vcf = arg.VCF_file
    all_chr_at_one = arg.fast_option
    chromosome_number = arg.chromosome_number
    thread_list = []
    list_of_chromosomes = config["CHROMOSOMES"].split(":")
    if chromosome_number == None:
        g1k_snp_dict = read_1000_genome(vcf, None)
        for chromosome_number in list_of_chromosomes:
            if all_chr_at_one == False:
                running_preprocessing(vcf, output_dir, sample_name, chromosome_number, g1k_snp_dict)
            else:
                print("submitting chr" + chromosome_number + " in parallel")
                t = threading.Thread(target=running_preprocessing,
                                     args=(vcf, output_dir, sample_name, chromosome_number, g1k_snp_dict))
                time.sleep(5)
                t.start()
                time.sleep(5)
                thread_list.append(t)
        time.sleep(60)
        for thread in thread_list:
            thread.join()

    else:
        g1k_snp_dict = read_1000_genome(vcf, chromosome_number)
        running_preprocessing(vcf, output_dir, sample_name, chromosome_number, g1k_snp_dict)


if __name__ == "__main__":
    Util.OsEnviron(os.environ)
    main()
