#! /usr/bin/env python

####Usage statement: This script will subset a vcf file returning only the contigs that
#### match a list of specified contigs
#### written by 18AdvBioinf
#### vcfsubsetter.py subsetlist.txt VCFIN.vcf

import sys

subsetlist=open(sys.argv[1], 'r')
vcfin=open(sys.argv[2], 'r')
vcfout=open('%s_subset.vcf'%(sys.argv[2][:-4]), 'w')

subsetdict={}
for line in subsetlist:
	line=line.rstrip()
	subsetdict[line]=''

linecount=0
for line in vcfin:
	linecount+=1
	line=line.rstrip()
	if line.startswith('#'):
		vcfout.write('%s\n'%(line))
	else:
		cols=line.split('\t')
		try:
		    subsetdict[cols[0]]
		    vcfout.write('%s\n'%(line))
		except KeyError:
			continue

subsetlist.close()
vcfin.close()
vcfout.close()
