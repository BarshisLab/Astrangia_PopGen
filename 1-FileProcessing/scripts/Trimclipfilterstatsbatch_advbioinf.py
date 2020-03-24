#!/usr/bin/env python
# Written by Dan Barshis

import sys, os

########### Usage #############
# Trimclipfilterstatsbatch.py barcodefile.txt anynumberoffastqfiles
# This should be run from within the folder with all of your original .fastqstrim
# Will quality trim multiple SINGLE-END fastq files
# Things to customize for your particular platform:
# qualoffset = 33 Quality score offset
# -t option in qualitytrim (the lower threshold quality score for trimming)
# -l option in qualitytrim and adapterclip (the length threshold for throwing out short reads)
# -q -p options in quality filter (the low quality, -q, and percentage -p options for filtration)

qualoffset=33
adapters=open(sys.argv[1], 'r')
adaptdict={}
for line in adapters:
	line=line.rstrip()
	cols=line.split('\t')
	adaptdict[cols[0]]=cols[1]

def adapterclip_batch_single(fastqs, outfile, qualoffset):
#    adapter='ATCACCGACTGCCCATAGAGAGG'   #this is the reverse complement of the trP1 Ion Torrent adapter
    for fastq in fastqs:
		fastq_prefix=fastq[:-6]
		o=open(outfile,'a')
		o.write('Adapter Clipping stats for: %s\n' %(fastq))
		o.close()
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastx_clipper -a %s -l 20 -Q %d -n -v -i %s -o %s >> %s" % (adaptdict[fastq], qualoffset, fastq, fastq_prefix+'_clipped.fastq', outfile)
		os.popen(cmd)
		print 'finished clipping: %s' %(fastq)

def qualitytrim_batch_single(fastqs, outfile, qualoffset):

	for fastq in fastqs:
		o=open(outfile,'a')
		o.write('Quality Trimming stats for: %s\n' %(fastq + '_clipped.fastq'))
		o.close()
		fastq_prefix=fastq[:-6]
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastq_quality_trimmer -v -Q %d -t 20 -l 20 -i %s -o %s >> %s" % (qualoffset, fastq_prefix + '_clipped.fastq', fastq_prefix + '_clippedtrimmed_nofilter.fastq', outfile)
		os.popen(cmd)
#		os.popen('rm %s' %(fastq_prefix+'_clipped.fastq'))
		print 'finished trimming: %s' %(fastq+ '_clipped.fastq')

# Will filter low quality reads from multiple SINGLE-END fastq files 
def qualfilter_batch_single(fastqs, outfile, qualoffset):
    for fastq in fastqs:
		fastq_prefix=fastq[:-6]
		o=open(outfile,'a')
		o.write('Quality Filtering stats for: %s\n' %(fastq_prefix+'_clippedtrimmed_nofilter.fastq'))
		o.close()
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastq_quality_filter -q 20 -p 90 -Q %d -v -i %s -o %s >> %s" % (qualoffset, fastq_prefix+'_clippedtrimmed_nofilter.fastq', fastq_prefix+'_clippedtrimmedfilterd.fastq', outfile)
		os.popen(cmd)
#		os.popen('rm %s' %(fastq_prefix+'_qualtrimmed.fastq'))
		print 'finished filtering: %s' %(fastq_prefix+'_clippedtrimmed_nofilter.fastq')

def qualstats_batch_single(fastqs):
    for fastq in fastqs:
		fastq_prefix=fastq[:-6]
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastx_quality_stats -Q %d -i %s -o %s " % (qualoffset, fastq_prefix+'_clippedtrimmedfilterd.fastq', fastq_prefix+'_clipped_trimmed_stats.txt')
		os.popen(cmd)
		print 'finished stats: %s' %(fastq_prefix+'_clippedtrimmedfilterd.fastq')

def makegraphs_batch_single(fastqs):
	for fastq in fastqs:
		fastq_prefix=fastq[:-6]
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastx_nucleotide_distribution_graph.sh -i %s -o %s -t %s" % (fastq_prefix+'_clipped_trimmed_stats.txt', fastq_prefix+'_nucdist.png', fastq_prefix)
		os.popen(cmd)
		cmd="/cm/shared/courses/dbarshis/15AdvBioinf/scripts/fastx_toolkit/fastq_quality_boxplot_graph.sh -i %s -o %s -t %s" % (fastq_prefix+'_clipped_trimmed_stats.txt', fastq_prefix+'_qualbox.png', fastq_prefix)
		os.popen(cmd)
		print 'finished graphs: %s' %(fastq_prefix)


adapterclip_batch_single(sys.argv[2:], 'trimclipstats.txt', qualoffset)
qualitytrim_batch_single(sys.argv[2:], 'trimclipstats.txt', qualoffset)
qualfilter_batch_single(sys.argv[2:], 'trimclipstats.txt', qualoffset)
qualstats_batch_single(sys.argv[2:])
makegraphs_batch_single(sys.argv[2:])
os.popen('mkdir ./QCFastqs')
os.popen('mkdir ./QCFastqs/nofilter')
os.popen('mv *_clippedtrimmed_nofilter.fastq ./QCFastqs/nofilter')
os.popen('mkdir ./QCFastqs/filtered')
os.popen('mv *_clippedtrimmedfilterd.fastq ./QCFastqs/filtered')
os.popen('mkdir ./filteringstats')
os.popen('mv *.png ./filteringstats')
os.popen('mv *stats.txt ./filteringstats')
os.popen('mkdir ./originalfastqs')
os.popen('mv %s ./originalfastqs'%(' '.join(sys.argv[1:])))
