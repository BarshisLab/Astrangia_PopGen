#! /usr/bin/env python

# convert_bayenv_to_baypass.py
# 7 April 2016 by C Tepolt
# 
# Reads in a file in BayEnv format,
# and creates a new file in BayPass format.

Usage = "USAGE: convert_bayenv_to_baypass.py BAYENV_FILE"

import numpy
import sys

if (len(sys.argv) < 2) or (len(sys.argv) > 2):
	print "Please provide one input file after the script name."
	print Usage
	exit()
else:
	BayEnvFile = sys.argv[1]

BayEnvFileName = BayEnvFile.strip('.txt')
OutFileName = BayEnvFileName + "-to-baypass.txt"

InputFile = open(BayEnvFile, 'r')
OutFile = open(OutFileName, 'w')

LineCounter = 1
LocusCount = 0

for Line in InputFile:
	Line = Line.strip('\n')
	Elements = Line.split("\t")

	if not LineCounter %2 == 0:
		AlleleOne = Elements
	else:
		AlleleTwo = Elements
		Interleaved = [x for t in zip(AlleleOne, AlleleTwo) for x in t]
		NewLine = '\t'.join(Interleaved)
		OutFile.write(NewLine + '\n')
		LocusCount += 1

	LineCounter += 1

print "Converted %d loci from BayEnv to BayPass format. Share and enjoy." % (LocusCount)	

# Clean up after your damn self.

InputFile.close()
OutFile.close()
