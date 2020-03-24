#!/usr/bin/env python

import sys
#Opens an infile specified by the user. Should be a list of the blast hits, each name on a new line
HITS = open(sys.argv[1], 'r')

#Opens a fasta file with the fasta sequences to extract from
DBSEQS = open(sys.argv[2], 'r')

#Opens an output text file as specified by user
OUT = open(sys.argv[3], 'w')

status=0
matches={}

linenum=0
for match in HITS:
	linenum+=1
	match=match.rstrip().strip(' ')
	cols=match.split(' ')
	if linenum>=1:
		matches[cols[1]]=''
	
#print matches

for line in DBSEQS:

	line=line.rstrip()
	if line[0] == '>':
	
#	For returning sequences IN the HitsList
#		if matches.has_key(line[1:]):
		if matches.has_key(line[1:].split(' ')[0]):
 			OUT.write(str(line) + '\n')
 			status = 1
 		else: status=0

#   For returning sequences NOT in the HitsList
#  		if matches.has_key(line[1:]):
# 			status=0
#  		else:
#  			OUT.write(str(line) + '\n')
#  			status = 1
 	else:
 		if status ==1:
 			OUT.write(str(line) + '\n')
