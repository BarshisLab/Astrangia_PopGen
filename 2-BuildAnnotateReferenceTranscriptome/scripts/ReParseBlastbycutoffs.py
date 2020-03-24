#!/usr/bin/env python

######################################################################
# This script takes a parsed blast output file which has a single line for each hit
# and returns the best top hit meeting either an evalue threshold or a
# length and percent identity threshold
# Usage is:
# outsuffix=sys.argv[1]  # string to be added to the output file
# parsetype=sys.argv[2] either 'evalue' or lengthidentity depending on which threshold
# if evalue:
# 	evalue = float(sys.argv[3]) # for parsing based on evalue threshold
# 	Blasts = sys.argv[4:]             # Any number of blast files for the same query seqs
# if lengthidentity:
# 	percentmatch=float(sys.argv[3])		  # the threshold of percent match
# 	lengthhit=float(sys.argv[4])  #  the threshold for number of nucleotides in overlap
# 	Blasts = sys.argv[5:]             # Any number of blast files for the same query seqs
######################################################################

import sys
filecount=0

outsuffix=sys.argv[1]  # string to be added to the output file

parsetype=sys.argv[2]  # either 'evalue' or 'lengthidentity' depending on which threshold

if parsetype=='evalue':
	evalue = float(sys.argv[3]) # for parsing based on evalue threshold
# match = sys.argv[2] # for counting the occurence of a string
	Blasts = sys.argv[4:]             # Any number of blast files for the same query seqs

if parsetype=='lengthidentity':
	percentmatch=float(sys.argv[3])		  # the threshold of percent match
	lengthhit=float(sys.argv[4])  #  the threshold for number of nucleotides in overlap
	Blasts = sys.argv[5:]             # Any number of blast files for the same query seqs

all=[]
contig=''

for file in Blasts:
	BLAST=open(file, 'r')
	linenum=0
	goodhits=0
	uniques={}
	filecount += 1
	OUT=open(file[:-4]+'_'+outsuffix, 'w')
	for line in BLAST:
		line=line.rstrip()
		linenum+=1
		if linenum == 1:
		
			OUT.write(line)	#Used to copy header to new files
	
		if linenum>1:	#Just used to skip past the header row

			cols=line.split('\t')
			
			if cols[0] == contig:
				continue
			else: 
#				if match in line:
				if parsetype=='evalue':
					if float(cols[12]) <= evalue: #for parsing based on evalue threshold
						OUT.write('\n'+line)    #Will be for outputting to a new parsed file
						goodhits = goodhits + 1
						contig = cols[0]
						if uniques.has_key(cols[2]) :
							continue
						else:
							uniques[cols[2]]=''
				if parsetype=='lengthidentity':
					if float(cols[14]) >= percentmatch and float(cols[4]) >= lengthhit: #for parsing based on threshold value
#					name=cols[2].split(" ")
#					for item in name:
#						if 'Match_Acc=' in item:
#							accession = item[10:]
						OUT.write('\n'+line)    #Will be for outputting to a new parsed file
						goodhits = goodhits + 1
						contig = cols[0]
						if uniques.has_key(cols[2]) :
							continue
						else:
							uniques[cols[2]]=''
	OUT.close()

	print file + "\tNumber of Good Hits:", "\t", goodhits
	print file + "\tNumber of unique matches:\t" + str(len(uniques.keys()))

