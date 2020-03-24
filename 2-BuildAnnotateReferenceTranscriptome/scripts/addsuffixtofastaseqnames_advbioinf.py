#! /usr/bin/env python
####usage addsuffixtofastaseqnames.py suffixtoadd infile.fasta
import sys

suffix = sys.argv[1]
fin = open(sys.argv[2], 'r')
outfile = open('%s_suffixed.fasta'%(sys.argv[2][:-6]), 'w')
count=0

for line in fin:
	line=line.rstrip()
	if line and line[0] == '>':                #indicates the name of the sequence
		count+=1
		outfile.write('%s_%s %s\n'%(line.split(' ')[0],suffix,line.split(' ')[1]))
	else:
		outfile.write('%s\n'%(line))
outfile.close()
