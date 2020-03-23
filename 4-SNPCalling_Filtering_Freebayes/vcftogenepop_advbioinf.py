#!/usr/bin/env python

#########################################################################################################################################
#                                                                                                                                       #
#   Written by Jason Ladner, modified by Dan Barshis                                                                                    #
#   questions?: jtladner@gmail.com                                                                                                      #
#   Usage:                                                                                                                              #
#   python vcftogenepop_advbioinf.py input_vcf_file inputPopulationsFile                                                                #
#                                                                                                                                       #
#   input_vcf_file = the .vcf file that you want to extract genotype information from                                                   #
#   inputPopulationsFile = this should contain two columns with no header, tab delimited text                                           #
#                          the first column should be the individual name                                                               #
#                          and the second column the population to which that individual belongs                                        #
#                                                                                                                                       #
#########################################################################################################################################


'''---------------------------> VCF Format <-----------------------------
	cols[0] = contig name
	cols[1] = SNP position
	cols[2] = ID???? - empty in my vcf files
	cols[3] = Reference base
	cols[4] = Alternative base
	cols[5] = SNP quality (Phred score)
	cols[6] = filter flags
	cols[7] = SNP Info
	cols[8] = genotype format (e.g., GT:AD:DP:GQ:PL)             ### This is different for some!!!!!!!!   It can even be different for different SNPs in the same .vcf
	cols[9:] = individual genotypes'''

import sys

#-------------------------------------------------------------#
#  This function takes a file and creates a dictionary with:  #
#    -one entry for each line                                 #
#    -keys will be first column of each line                  #
#    -values will be the complete line                        #
#-------------------------------------------------------------#
def file_dict(file):
	fin = open(file, 'r')       
	d={}

	for line in fin:
		line=line.rstrip()
		cols=line.split()
		d[cols[0]]=line
	fin.close()
	return d


#-------------------------------------------------------------#
#  This function takes the genotype description for a SNP     #
#   (e.g., GT:AD:DP:GQ:PL)                                    #
#  and returns the position of the genotype quality (GQ)      #
#  with the first position being 0                            #
#-------------------------------------------------------------#
def get_GQ_pos(genotype_description):            
	return genotype_description.split(':').index('GQ')

#-------------------------------------------------------------#
#  This function writes genotype data to an outfile           #
#   with 2 columns per SNP (one for each allele)              #
#                                                             #
#  This is very similar to the input required by smartpca     #
#-------------------------------------------------------------#
def write_snps_as_cols(outfile_name, headers, contigs, pos, ref, alt, all_gens, pops):
	#Opens an output text file as specified by user
	OUT = open(outfile_name, 'w')

	OUT.write('AllSNPs\n')
	OUT.write('\t')
	for index in range(len(contigs))[:-1]:
		OUT.write('%s_%s,\t' % (contigs[index],pos[index]))          #Can be easily changed so that the contig name is a combination of the contig_pos
	OUT.write('%s_%s' % (contigs[index+1],pos[index+1]))
	for pop in set(pops.values()):
		OUT.write('\nPop')
		indivcount=0
		for indiv in headers[4:]:
			if pops[indiv] == pop:
				OUT.write('\n%s,' % (indiv))
				for index in range(len(contigs)):
					OUT.write('\t')
					for each in all_gens[indivcount][index]:
						if each=='.':
							OUT.write('%s' %('00'))
						if each=='0':
							OUT.write('%s' %('01'))
						if each=='1':
							OUT.write('%s' %('02'))
			indivcount+=1

	OUT.close()

def extract_genos_from_vcf(vcf_file, outfile_name, popsfile):
	popsfilein=open(popsfile, 'r')
	pops={}
	for popline in popsfilein:
		cols=popline.rstrip().split('\t')
		pops[cols[0]]=cols[1]
	popsfilein.close()
	#Opens an infile specified by the user. Should be a .vcf file 
	fin = open(sys.argv[1], 'r')
	linecount=0
	contigs=[]
	pos=[]
	ref=[]
	alt=[]
	all_gens=[]
	locus=[]
	headers=[]

	#Reads through the vcf file line by line
	for line in fin:
		cols=line.rstrip().split('\t')
		gens=[]                           #Will hold individual level genotype data, resets list from previous SNP
		gencount=0                        #Will be used below to keep track of the individuals when appending new SNPs
	
		if cols[0][0] == '#' and cols[0][1] != '#':       #To pull out headers

			headers=cols[0:2] + cols[3:5] + cols[9:]      #Creates list of headers of interest
			print "Indivs with genotypes in vcf file: %s" % ("\t".join(headers[4:]))                             #Prints the IDs of the individuals with genotypes in the vcf

		if cols[0][0] != '#':                             #Specifies only SNP lines
			linecount+=1                                  #Keeps track of the numbers of SNPs
			
			contigs.append(cols[0])   #Creates a list of all the contig names
			pos.append(cols[1])       #Creates a list of all the SNP base positions
			ref.append(cols[3])       #Creates a list of all the reference bases
			alt.append(cols[4])       #Creates a list of all the alternative bases

			genotype_description=cols[8]
#			GQ_position = get_GQ_pos(genotype_description)
			
			gens=cols[9:]                  #Creates a list of all the genotype info
			
			if linecount==1:                      #will only execute these commands once, when reading the first line of the file
				for gen in gens:
					all_gens.append([])           #Creates seperate list entry for each individual, will hold genotypes at all SNPs
					locus.append([])              #Creates seperate list entry for each individual, will hold genoypes at a specific locus 
			
			for gen in gens:                    #Steps through the genotype info for each individual
				parts=gen.split(':')            #Breaks apart the different sections of the genotype info
				genotype=parts[0]               #Most probable genotype for the individual, with the two alleles seperated by '/'
				genotypes=genotype.split('/')
	
				locus[gencount]=genotypes           #Replace genotype data from last locus with genotype from this locus
				gencount+=1                         #Keeping track of the individual whose data is being examined
		
			locuscount=0
			for each in locus:
				all_gens[locuscount].append(each)                  #Appends genotype data to the all_gens list
				locuscount+=1                                      #keeps track of the number of loci with no missing data

	#--------------->'cols-style output is what you need to go into smartpca, but it still needs to be *****tweaked to include population*****
	#This will have the contig name, pos, ref and alt alleles listed as rows at the top and then each individual and all of its genotyoes listed in the rows that follow. There will be two columns of data for each SNP
	print len(headers), len(contigs), len(pos), len(ref), len(alt), len(all_gens)
	write_snps_as_cols(outfile_name, headers, contigs, pos,  ref, alt, all_gens, pops)


###------->>>

extract_genos_from_vcf(sys.argv[1], '%s_genepop.gen'%(sys.argv[1][:-4]), sys.argv[2])
