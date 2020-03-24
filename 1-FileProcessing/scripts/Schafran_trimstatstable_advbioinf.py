#! /usr/bin/python
# Written by Peter Schafran 1-7 Oct 2015
#
# Usage: Schafran_trimstatstable.py [-c, -v, -h ] trimclipstats.txt outputfilename.ext

import sys,re
#verbose initial state off
verbose=0
#check for -v parameter, if yes, turn on
if '-v' in sys.argv:
	verbose=1
#display help message with -h parameter
if '-h' in sys.argv:
	print '''Written by Peter Schafran pscha005@odu.edu 5-Oct-2015

This script takes a stats output file from fastx_clipper and converts it into a table.

Usage: Schafran_trimstatstable.py [-c, -v, -h] inputfile.txt outputfile.txt

Options (-c and -v must be listed separately to run together):
-h	Display this help message
-c	Use comma delimiter instead of tabs
-v	Verbose mode (print steps to stdout)
'''
else:
	#Get input file from 2nd to last position of sys.argv
	infile = open(sys.argv[-2], 'r')

	#Create empty dictionary to store all stats
	statdict={}

	#Loop to pull filenames and put in dictionary
	for line in infile:
		if 'fastq' in line:
			#split lines on various characters to isolate just the sample ID
			splitline = line.strip('\n').split(':')
			splitline2 = splitline[1].split('.')
			splitline3 = splitline2[0].split('_clipped')
			filename = splitline3[0].strip( )
			#store sample IDs in dictionary 
			statdict[filename]=[]
			if verbose ==1:
				print 'Adding %s to dictionary...' %(filename)
	#script for some reason won't enter 'for' loop without redefining infile
	infile = open(sys.argv[-2], 'r')
	#loop through each line in stats file, compare to keys in dictionary to match 
	for line in infile:
		linecount=0
		#set currentkey variable so each sampleID is retained through all the lines with its stats
		for key in statdict.keys():
			if key in line:
				currentkey=key
				linecount+=1
		#use regex to pull adapter seq or numbers from each line
		splitline = re.findall('[ACTG\d]+',line)
		#use line counter to avoid storing lines with filename 
		if linecount==0:
			#sort list of data found in re.findall from longest to shortest. Longest item is one we want to keep
			data=sorted(splitline,key=len, reverse=True)
			#data in first position of sorted list appended to dictionary
			statdict[currentkey].append(data[0])
			if verbose==1:
					print "Appending %s to dictionary key %s..." %(data[0],currentkey)
	infile.close()
	#open output file with name from last position of sys.argv
	outfile=open(sys.argv[-1],'w')
	#check for comma-delimiting option
	if '-c' in sys.argv:
		if verbose==1:
			print "Writing output file..."
		#write header line
		outfile.write('OriginalFileName,ClippingAdapter,MinAdapterLength,AdapterInput,AdapterOutput,AdapterTooShort,AdapterOnly,MinQualThresh,MinQualLength,QualInput,QualOutput,Qualdiscarded,Quality_cut-off,Minpercentage,QualFilterInput,QualFilterOutput,QualFilterdiscarded\n')
		#loop through keys in statdict, sorted alphabetically
		for key in sorted(statdict.keys()):
			outfile.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(key, statdict[key][0],statdict[key][1],statdict[key][2],statdict[key][3],statdict[key][4],statdict[key][5],statdict[key][6],statdict[key][7],statdict[key][8],statdict[key][9],statdict[key][10],statdict[key][11],statdict[key][12],statdict[key][13],statdict[key][14],statdict[key][15]))
	#default mode is tab-delimited
	else:
		if verbose==1:
			print "Writing output file..."
		#write header line
		outfile.write('OriginalFileName	ClippingAdapter	MinAdapterLength	AdapterInput	AdapterOutput	AdapterTooShort	AdapterOnly	MinQualThresh	MinQualLength	QualInput	QualOutput	Qualdiscarded	Quality_cut-off	Minpercentage	QualFilterInput	QualFilterOutput	QualFilterdiscarded\n')
		#loop through keys in statdict, sorted alphabetically
		for key in sorted(statdict.keys()):
			#write data line to file from dict entries
			outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(key, statdict[key][0],statdict[key][1],statdict[key][2],statdict[key][3],statdict[key][4],statdict[key][5],statdict[key][6],statdict[key][7],statdict[key][8],statdict[key][9],statdict[key][10],statdict[key][11],statdict[key][12],statdict[key][13],statdict[key][14],statdict[key][15]))
	outfile.close()
