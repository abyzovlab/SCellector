"""Generates allele frequency of variants using Bam file and vcf files"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Development"

import sys
import threading
import time
import os
import argparse
import Util
import gzip
import re

from subprocess import PIPE, Popen


def cmdline(command):
    """Executes the command and returns the result"""
    process = Popen(args=command, stdout=PIPE, shell=True)
    return process.communicate()[0]


def argument_parse(script_path):
    """Parses the command line arguments"""
    parser = argparse.ArgumentParser(description='Preprocessing of vcf file')
    parser.add_argument("-B", "--BAM_file", help="Path to bam file", required=True)
    parser.add_argument("-V", "--VCF_file", help="Path to VCF file", required=True, type=Util.FileValidator)
    parser.add_argument("-O", "--Output_dir", help="Path to directory where results will be written", required=True)
    parser.add_argument("-S", "--Sample_name", help="Name of the sample", required=True)
    parser.add_argument("-q", "--min_base_quality", help="Minimum base quality", default="20")
    parser.add_argument("-Q", "--min_mapping_quality", help="Minimum mapping quality", default="20")
    parser.add_argument("-l", "--chromosome_number", help="chromosome number ")
    parser.add_argument("-f", "--all_chromosome_in_parallel", help="run all chromosome in parallel ",
                        action='store_true')
    parser.add_argument("-n", "--no_chr_splitting", help="Doesn't split by chromosome", action='store_true')
    return parser


def rmv_indels(info):
    """
    Takes a string of bases with INDELs and returns the strings after removing INDELs
    Can be used to detect the Indels in future
    """
    if "+" in info or "-" in info:
        tmp_info = ""
        j = 0
        while j < len(info):
            if info[j] == "+" or info[j] == "-":
                if j + 3 < len(info) and info[j + 3].isdigit() and info[j + 2].isdigit():  # new adds
                    indel = int(info[j + 1] + info[j + 2] + info[j + 3]) + 2  # new adds
                elif j + 2 < len(info) and info[j + 2].isdigit():  # new add (only the if)
                    indel = int(info[j + 1] + info[j + 2]) + 1
                else:
                    indel = info[j + 1]
                skip = int(indel) + 2
                j = j + skip
                try:
                    tmp_info = tmp_info + info[j]
                except IndexError:
                    continue
            else:
                tmp_info = tmp_info + info[j]
            j = j + 1
        info = tmp_info
    return info

def snp_unit_calculation(bam):
    cmd = " ".join([samtools, "view", bam,
                    "| head -100 | cut -f10 | awk '{print length($1)}'|sort |uniq -c | sort -rnk 1,1| awk '{print $2}'| head -1"])
    read_length = float(cmdline(cmd))
    cmd = " ".join([samtools, "idxstats", bam,
                    "| awk '{sum=sum+$3;print sum}'| tail -1"])
    total_reads = float(cmdline(cmd))
    bases_in_human_genome = float("3088286401")
    coverage = float((total_reads * read_length) / bases_in_human_genome)
    snp_unit = int(30 / coverage)
    return snp_unit

def base_quality_check(info, base_info, min_base_quality):
    """
    Does a base quality check.
    Removes bases lower than a base quality score provided by the user
    """
    base = ""
    qual = ""
    if len(info) == 0:
        return "", ""
    if not len(info) == len(base_info):
        print("Number of read bases doesn't match number of base quality score")
        exit()
    for i in range(0, len(info)):
        if ord(base_info[i]) < int(min_base_quality):
            continue
        if info[i] == "*":
            continue
        base += info[i]
        qual += base_info[i]
    return base, qual


def create_pileup_and_parse(bam_file, vcf, output_dir, sample_name, chromosome_number, minimum_base_quality,
                            minimum_mapping_quality):
    '''creates a pileup file, parses it and stores it in a file'''
    global samtools
    samtools = config["SAMTOOLS"]
    reference = config["REFERENCE"]
    pileup_file = os.path.join(output_dir, sample_name + ".pileup")
    output_file = os.path.join(output_dir, sample_name + ".AF.txt")
    pileup_file_pos = os.path.join(output_dir, sample_name + ".pileup.pos")
    cmd = " ".join(["cat", vcf, "|grep -v \"#\" |awk '{print $1\"\\t\"$2}' >", pileup_file_pos])
    print(cmd)
    os.system(cmd)
    cmd = " ".join([samtools, "mpileup",
                    "-q", minimum_base_quality,
                    "-Q", minimum_mapping_quality,
                    "-d 10000",
                    bam_file,
                    "-f", reference,
                    "-l", pileup_file_pos,
                    "-o", pileup_file])
    ### with chromosome option
    if chromosome_number != None:
        output_file = os.path.join(output_dir, sample_name + ".chr" + chromosome_number + ".AF.txt")
        pileup_file = os.path.join(output_dir, sample_name + ".chr" + chromosome_number + ".pileup")
        pileup_file_pos = os.path.join(output_dir, sample_name + ".chr" + chromosome_number + ".pileup.pos")
        if Util.check_vcf_with_chr_or_not(vcf) == "with_chr":
            chromosome_number = "chr" + chromosome_number
        cmd = " ".join(
            ["cat", vcf, "|grep -v \"#\" |grep -w \"^" + chromosome_number + "\"", "|awk '{print $1\"\\t\"$2}' >",
             pileup_file_pos])
        print(cmd)
        os.system(cmd)
        cmd = " ".join([samtools, "mpileup",
                        "-q", minimum_base_quality,
                        "-Q", minimum_mapping_quality,
                        "-d 10000",
                        bam_file,
                        "-f", reference,
                        "-l", pileup_file_pos,
                        "-o", pileup_file])
    print(cmd)
    os.system(cmd)
    vcf_dict = {}
    if vcf.endswith(".gz"):
        vcf_fh = gzip.open(vcf)
    else:
        vcf_fh = open(vcf)
    for i in vcf_fh:
        if i.startswith("#"):
            continue
        line = i.strip().split("\t")
        chr_pos = line[0] + "_" + line[1]
        ref = line[3]
        alt = line[4]
        vcf_dict[chr_pos] = [ref, alt]
    vcf_fh.close()
    output_file_fh = open(output_file, "w")
    output_file_fh.write(
        "#Chr\tPos\tRef\tAlt\tDepth\tA\ta\tT\tt\tG\tg\tC\tc\tref_depth\talt_depth\tref_freq\talt_freq\n")

    ### parsing pileup files
    output_dict = {}
    for i in open(pileup_file):
        line = i.strip().split("\t")
        chr_pos = line[0] + "_" + line[1]
        alt = vcf_dict[chr_pos][1]
        line1 = line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + alt
        tot_depth = line[3]
        line2 = ""
        rd = 3
        ref_depth = 0
        alt_depth = 0
        info_raw = rmv_indels(re.sub(r'\^.', "", i.split("\t")[rd + 1]).replace("$", ""))
        info, base_qual = base_quality_check(info_raw, i.split("\t")[rd + 2].strip(), minimum_base_quality)
        ref = line[2]
        line2 = line1 + "\t" + i.split("\t")[rd].strip()
        for k in "ATGC":
            if k == ref:
                line2 = line2 + "\t" + str(info.count(".")) + "\t" + str(info.count(","))
                ref_depth = info.count(".") + info.count(",")
                continue
            if k == alt:
                alt_depth = info.count(k) + info.count(k.lower())
            line2 = line2 + "\t" + str(info.count(k)) + "\t" + str(info.count(k.lower()))
        if tot_depth == "0":
            ref_freq = "0.0"
            alt_freq = "0.0"
        else:
            ref_freq = float(float(ref_depth) / float(tot_depth))
            alt_freq = float(float(alt_depth) / float(tot_depth))
        line2 = line2 + "\t" + str(ref_depth) + "\t" + str(alt_depth) + "\t" + str(ref_freq) + "\t" + str(alt_freq)
        output_dict[chr_pos] = line2
    for i in open(pileup_file_pos):
        line = i.strip().split("\t")
        chrm = line[0]
        pos = line[1]
        chr_pos = chrm + "_" + pos
        ref = vcf_dict[chr_pos][0]
        alt = vcf_dict[chr_pos][1]
        if chr_pos in output_dict:
            output_file_fh.write(output_dict[chr_pos] + "\n")
        else:
            line_out = "\t".join(
                [chrm, pos, ref, alt, "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0.0", "0.0"])
            output_file_fh.write(line_out + "\n")

    output_file_fh.close()


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
            Util.SendEmail(script_name, sys.exc_info()[1], os.environ['SGE_STDERR_PATH'])
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
            Util.SendEmail(script_name, sys.exc_info()[1], os.environ['SGE_STDERR_PATH'])
            exit(100)
    # assigning values to variable
    output_dir = arg.Output_dir
    sample_name = arg.Sample_name
    vcf = arg.VCF_file
    all_chr_at_one = arg.all_chromosome_in_parallel
    chromosome_number = arg.chromosome_number
    bam_file = arg.BAM_file
    minimum_base_quality = arg.min_base_quality
    minimum_mapping_quality = arg.min_mapping_quality
    all_chr_at_one = arg.all_chromosome_in_parallel
    no_split = arg.no_chr_splitting
    list_of_chromosomes = config["CHROMOSOMES"].split(":")
    # creating pileup files
    thread_list = []
    if chromosome_number == None and no_split == False:
        for chromosome_number in list_of_chromosomes:
            if all_chr_at_one == False:
                print("Running chr" + chromosome_number)
                create_pileup_and_parse(bam_file, vcf, output_dir, sample_name, chromosome_number, minimum_base_quality,
                                        minimum_mapping_quality)
            else:
                print("Running chr" + chromosome_number + " in parallel-------------")
                t = threading.Thread(target=create_pileup_and_parse, args=(
                bam_file, vcf, output_dir, sample_name, chromosome_number, minimum_base_quality,
                minimum_mapping_quality))
                thread_list.append(t)
                time.sleep(5)
                t.start()
        time.sleep(60)
        for thread in thread_list:
            thread.join()
    else:
        create_pileup_and_parse(bam_file, vcf, output_dir, sample_name, chromosome_number, minimum_base_quality,
                                minimum_mapping_quality)
    snp_unit = snp_unit_calculation(bam_file)
    print("Recommended SNP unit to use for whole genome human sample=", snp_unit)


if __name__ == "__main__":
    main()