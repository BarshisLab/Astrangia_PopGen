#! /usr/bin/env python

###Usage####
#Authors Dan Barshis
import sys, os

#argv[1]=Inputfasta name
#argv[2]=Number of seq per splitfile
#argv[3]=BlastDatabase path, names, output formats (only 5 or 7), and numthreads

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta(file):
	fin = open(file, 'r')
	count=0
	
	names=[]
	seqs=[]
	seq=''
	for line in fin:
		line=line.strip()
		if line and line[0] == '>':                #indicates the name of the sequence
			count+=1
			names.append(line[1:])
			if count>1:
				seqs.append(seq)
			seq=''
		else: seq +=line
	seqs.append(seq)
	
	return names, seqs

def get_fasta(maxperfile, names, seqs):
	indices=[]
	splitfilelist=[]
	second=maxperfile
	first=0
	total_digits = len(str(len(names)))
	while True:
		indices.append([first, second])
		first+=maxperfile
		second+=maxperfile
		if second > len(names):
			second=len(names)
			indices.append([first, second])
			break
	count=1
	
	for i in indices:
		sub_names=names[i[0]:i[1]]
		sub_seqs=seqs[i[0]:i[1]]
		last = count+(int(maxperfile)-1)
		if last > len(names): last=len(names)
		o=open("./splitfastas/%s_%s%d-%s%d.fasta" % (sys.argv[1].split('.')[0], (total_digits-len(str(count)))*'0',count,(total_digits-len(str(last)))*'0' ,last), 'w')
		splitfilelist.append("%s_%s%d-%s%d.fasta" % (sys.argv[1].split('.')[0], (total_digits-len(str(count)))*'0',count,(total_digits-len(str(last)))*'0' ,last))
		count+=int(maxperfile)
		for j in range(len(sub_names)):
			o.write('>' + sub_names[j] + '\n')
			o.write(sub_seqs[j] + '\n')
		o.close()
	return splitfilelist

os.popen("mkdir ./splitfastas")
names, seqs=read_fasta(sys.argv[1])
splitfilelist=get_fasta(int(sys.argv[2]), names, seqs)

# now making individual qsub scripts for each blast job
dbparams=open(sys.argv[3], 'r')
os.popen("mkdir ./submissionscripts ./blastoutputs")
linecount=0
qsubslist=[]
for line in dbparams:
	cols=line.rstrip().split('\t')
	linecount+=1
	if linecount>1:
		dbpath=cols[0]
		db=cols[1]
		if cols[2]=='7':
			suffix='.txt'
		if cols[2]=='5':
			suffix='.xml'
		os.popen("mkdir	./blastoutputs/%s" %(db))
		outfmt=cols[2]
		numthreads=cols[3]
		for fasta in splitfilelist:
			header='#SBATCH -o blastx2%s.txt\n#SBATCH -n %s\n#SBATCH --mail-user=haich001@odu.edu\n#SBATCH --mail-type=END\n#SBATCH --job-name=splitblast_hea\nmodule load blast/2.6.0\n' %(db, numthreads)
			outfile=open('./submissionscripts/%s_blastsub_%s.sh' %(fasta[:-6],db),'w')
			qsubslist.append('./submissionscripts/%s_blastsub_%s.sh' %(fasta[:-6],db))
			outfile.write('%s\n\n' %('#!/bin/bash -l'))
#For blastn
#			outfile.write('%sblastn -db %s%s -query ./splitfastas/%s -out ./blastoutputs/%s/%s_blastx2%s%s -evalue 0.05 -outfmt %s -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -word_size 11 -max_target_seqs 20 -num_threads %s\n' %(header,dbpath,db,fasta,db,fasta[:-6],db,suffix, outfmt, numthreads))
#For blastx
			outfile.write('%sblastx -db %s%s -query ./splitfastas/%s -out ./blastoutputs/%s/%s_blastx2%s%s -evalue 0.05 -outfmt %s -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads %s\n' %(header,dbpath, db,fasta,db,fasta[:-6],db,suffix, outfmt, numthreads))
			outfile.close()
subcount=0
batchdir=1
for sub in qsubslist:
	subcount+=1
	if subcount==1:
		os.popen("mkdir ./submissionscripts/qsubsbatch_%d"%(batchdir))
		os.popen("mv %s ./submissionscripts/qsubsbatch_%d/"%(sub,batchdir))
	if subcount>1 and subcount<20:
		os.popen("mv %s ./submissionscripts/qsubsbatch_%d/"%(sub,batchdir))
	if subcount==20:
		os.popen("mv %s ./submissionscripts/qsubsbatch_%d/"%(sub,batchdir))
		batchdir+=1
		subcount=0
