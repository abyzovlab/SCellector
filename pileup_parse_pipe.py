"""Parses the pileup file as STDIN input,to give counts of bases for each position provided"""
__author__ = "Vivekananda Sarangi"
__email__ = "sarangi.vivekananda@mayo.edu"
__status__ = "Validation"

import sys
import re
import math
import os

from subprocess import PIPE, Popen


def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]

try:
	sample_name=sys.argv[1]
	min_base_quality=int(sys.argv[2])
	sys.stdin
except:
	print "Usage -> samtools mpileup bam_file | python pileup_parse_pipe.py <sample_name> <Minimum_base_quality> <position_bed_file>"
	exit()

def average(s): return sum(s)*1.0/ len(s)

def std(n):
	'''
	Finds the standard deviation of a base quality string
    	Takes a string as input and returns average and standard deviation as output
	Makes call to another function calles average()
	'''
	qu=[]
	for i in n:
		qu.append(ord(i))
        avg = average(qu)
	variance = map(lambda x: (x - avg)**2,qu)
	standard_deviation = math.sqrt(average(variance))
	return avg,standard_deviation

def rmv_indels(info):
	'''
	Takes a string of bases with INDELs and returns the strings after removing INDELs
	Can be used to detect the Indels in future 
	'''
	if "+" in info or "-" in info:
        	tmp_info=""
		j=0
                while j<len(info):
                	if info[j]=="+" or  info[j]=="-":
				if j+3 < len(info) and info[j+3].isdigit() and info[j+2].isdigit():  # new adds
					indel=int(info[j+1]+info[j+2]+info[j+3])+2  #new adds
				elif j+2 < len(info) and info[j+2].isdigit():# new add (only the if)
                                        indel=int(info[j+1]+info[j+2])+1
				else:
					indel=info[j+1]
                        	skip=int(indel) + 2
                                j=j+skip
                                try:
                                        tmp_info=tmp_info+info[j]
                                except IndexError:
                                       continue
			else:
				tmp_info=tmp_info+info[j]
			j=j+1
                info=tmp_info
	return info

def base_quality_check(info,base_info,min_base_quality):
	'''
	Does a base quality check.
	Removes bases lower than a base quality score provided by the user
	Future developement include taking only phred score as input
	'''
	base=""
	qual=""
	if len(info)==0:
		return "",""
	if not len(info)==len(base_info):
		print "Number of read bases doesn't match number of base quality score"
		exit()
	for i in range(0,len(info)):
		if ord(base_info[i]) < min_base_quality:
			continue
		if info[i]=="*":
			continue
		base+=info[i]
		qual+=base_info[i]
	return base,qual

def main():
	sample=''
	position=[]
	for i in  sys.stdin:
		line=i.split("\t")[0]+"\t"+i.split("\t")[1]+"\t"+i.split("\t")[2]
		position.append(i.split("\t")[0]+"_"+i.split("\t")[1])
		line2=""
		rd=3
		info_raw=rmv_indels(re.sub(r'\^.', "", i.split("\t")[rd+1]).replace("$",""))
		info,base_qual=base_quality_check(info_raw,i.split("\t")[rd+2].strip(),min_base_quality)
		if i.split("\t")[rd]=="0" or len(info)==0:
			mean=0.0
			stdv=0.0
		else:
			mean,stdv=std(base_qual)
		ref=i.split("\t")[2]
		line2=line+"\t"+i.split("\t")[rd].strip()
		for k in "ATGC":
			if k==ref:
				line2=line2+"\t"+str(info.count("."))+"\t"+str(info.count(","))
				continue
			line2=line2+"\t"+str(info.count(k))+"\t"+str(info.count(k.lower()))
		line2=line2+"\t"+str(info.count("<")+info.count(">"))
		sample=line2.split("\t")[0]
		print line2.strip()
	pos=sys.argv[3]
	if pos not in position:
		cmd=" ".join(["/projects/bsi/bictools/apps/alignment/samtools/0.1.19/samtools faidx /data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa",pos.split("_")[0]+":"+pos.split("_")[1].strip()+"-"+pos.split("_")[1].strip(),"|grep -v \">\""])
		ref_base=cmdline(cmd).strip()
		line2="\t".join([pos.split("_")[0],pos.split("_")[1].strip(),ref_base,"0","0","0","0","0","0","0","0","0","0"])
		print line2
if __name__ == "__main__":
        main()
