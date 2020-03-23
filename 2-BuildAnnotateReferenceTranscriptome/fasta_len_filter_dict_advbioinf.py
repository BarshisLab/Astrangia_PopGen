#! /usr/bin/env python

import sys
import os

#This script takes any number of input fasta files and outputs sequences that are >= a user specified length threshold
#####USAGE fasta_len_filter.py lengththresh(i.e.200) anynumberofinputfiles.fasta

def read_fasta(file):
	fin = open(file, 'r')
	count=0
	
	contigs={}
	seq=''
	for line in fin:
		line=line.strip()
		if line and line[0] == '>':                #indicates the name of the sequence
			count+=1
			if count>1:
				contigs[name]+=seq
			name=line[1:]
			contigs[name]=''
			seq=''
		else:
			seq +=line
	contigs[name]+=seq
	fin.close()
	return contigs
	
fastas=sys.argv[2:]

for file in fastas:
	contigs=read_fasta(file)
	fin=open(file,'r')
	lengthresh = int(sys.argv[1])
	ofilestring='%s_over%d%s'%(file[:-6],lengthresh,'.fasta')
	outfile= open(ofilestring, 'w')
	overthresh=1
	print 'Number of total seqs for %s: %d' %(file,len(contigs.keys()))
	for line in fin:
		line=line.strip()
		if line and line[0] == '>':
			name=line[1:]
			if len(contigs[name]) >= lengthresh:
				overthresh+=1
				outfile.write('>%s\n%s\n'%(name,contigs[name]))
		else:
			continue
	fin.close()
	outfile.close()
	print 'Number of seqs over %d for %s: %d' %(lengthresh, file, overthresh)
