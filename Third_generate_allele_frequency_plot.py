"""Combines the data from phasing file and allele frequency file from first and second script to generate allele frequency plot"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu,viveksarangi@gmail.com"

import os
import argparse
import sys
import Util
import numpy as np
import matplotlib
from scipy.stats import norm

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def argument_parse():
    '''Parses the command line arguments'''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-a", "--AF_file", help="path sample allele freq file", required=True, type=Util.FileValidator)
    parser.add_argument("-g", "--Germ_hap", help="Path to germline hap file", required=True, type=Util.FileValidator)
    parser.add_argument("-o", "--Output_dir", help="Path to output file", required=True)
    parser.add_argument("-S", "--Sample_name", help="Name of the sample", required=True)
    parser.add_argument("-n", "--Snps", help="Number_of_snps", default="100")
    return parser


def plot_data(output, sample_name):
    """plots allele frequency plot"""
    plot_out = output.replace(".txt", ".png")
    af_numbers = []
    head = {}
    for i in open(output):
        line = i.strip().split("\t")
        if i.startswith("#"):
            for n, j in enumerate(line):
                head[j] = n
            continue
        af_numbers.append(float(line[head["h1_af"]]))
    sigma = np.std(af_numbers)
    mu = np.mean(af_numbers)
    n, bins, patches = plt.hist(af_numbers, bins=31, range=(-0.01, 1.01), edgecolor='steelblue', color='steelblue',
                                density='true')
    y = norm.pdf(bins, mu, sigma)
    if sigma <= 0.25:
        plt.plot(bins, y, 'g-', linewidth=1.5)
    elif sigma <= 0.33:
        plt.plot(bins, y, 'y-', linewidth=1.5)
    else:
        plt.plot(bins, y, 'r-', linewidth=1.5)
    plt.tick_params(labelsize=10, left='off', top='off', right='off', bottom='off')  # removing tick marks
    sample = sample_name.split("_")[0].split(".")[0]
    plt.title(sample, fontsize=20, loc='left')  # sample name
    plt.title("std=" + str(sigma)[:4], fontsize=20, loc='right')  # std valu
    # lines for creating percent on y axis and limiting it to 30%
    sum_n = sum(n)
    plt.xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    plt.yticks([sum_n * 0.0, sum_n * 0.05, sum_n * 0.10, sum_n * 0.15, sum_n * 0.20, sum_n * 0.25, sum_n * 0.30])
    formatter = mticker.FuncFormatter(lambda n, y: str((n / sum_n))[:4])
    n_for_y = sum_n * 0.3
    plt.ylim(0, n_for_y)  ## setting the upper lim of y axis
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.xlabel("Allele Frequency")
    plt.ylabel("Percent of SNP units")
    plt.tight_layout()  ## make every subplot fit properly
    plt.grid(True, color='steelblue', linestyle="-", alpha=0.5)  ## adds grid
    plt.savefig(plot_out)


def main():
    parser = argument_parse()
    arg = parser.parse_args()
    script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    output_dir = arg.Output_dir
    sample_name = arg.Sample_name
    output = os.path.join(output_dir, sample_name + ".hap_af.txt")
    af_file = arg.AF_file
    germ_hap = arg.Germ_hap
    number_of_snps = int(arg.Snps)
    hap_dict = {}
    for i in open(germ_hap):
        if i.startswith("#"):
            continue
        line = i.strip().split("\t")
        chr_pos = line[0] + "_" + line[1]
        format_col = line[8]
        genotype_col = line[9]
        for n, j in enumerate(format_col.split(":")):
            if j == "GT":
                haplotype = genotype_col.split(":")[n].split("|")[0]
        hap_dict[chr_pos] = haplotype
    snp_increament = 0
    n_h1 = 0
    n_h2 = 0
    depth = 0
    output_fh = open(output, 'w')
    output_fh.write("#chr\tpos\tn_snps\th1_af\tn_h1\tnh2\n")
    for i in open(af_file):
        line = i.strip().split("\t")
        if i.startswith("#"):
            continue
        if snp_increament == 0:
            chr_first = line[0]
            pos_first = line[1]
        if snp_increament == number_of_snps or chr_first != line[0]:
            if depth == 0:
                h1_af = 0.0
            else:
                h1_af = float(float(n_h1) / float(depth))
            out = "\t".join(
                [chr_first, pos_first + "-" + pos_end, str(snp_increament), str(h1_af), str(n_h1), str(n_h2)])
            output_fh.write(out + "\n")
            snp_increament = 0
            n_h1 = 0
            n_h2 = 0
            depth = 0
            chr_first = line[0]
            pos_first = line[1]

        snp_increament += 1

        chr_pos = line[0] + "_" + line[1]
        ref_depth = line[13]
        alt_depth = line[14]
        depth += int(line[4])
        if hap_dict[chr_pos] == "0":
            n_h1 += int(ref_depth)
            n_h2 += int(alt_depth)
        else:
            n_h1 += int(alt_depth)
            n_h2 += int(ref_depth)
        pos_end = line[1]
    output_fh.close()

    plot_data(output, sample_name)


if __name__ == "__main__":
    Util.OsEnviron(os.environ)
    main()
