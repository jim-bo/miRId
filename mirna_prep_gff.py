#!/usr/bin/python
'''
makes a gff files of regions of the genome hit by
a sam entry.

usage: reference.fasta output.gff [-window integer] sam_file(s)

'''

# imports.
import sys
import os
import numpy as np
from bitarray import bitarray
import biobuffers
import logging

# logging.
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

# parameters.
fasta_file = sys.argv[1]
gff_file = sys.argv[2]

if sys.argv[3] == "-window":
	window_size = int(sys.argv[4])
	start = 5
else:
	window_size = 0
	start = 3
	
sam_files = []
for i in range(start, len(sys.argv)):
	sam_files.append(sys.argv[i])
	
############## functions ##############
	
def create_fasta(ref_file):
	''' reads contig sizes '''
	
	# create dictionary.
	ref_fasta = dict()
	
	# count number of entries and get their size.
	for entry in biobuffers.fasta_buffer(ref_file):

		# make entry.
		ref_fasta[entry['name']] = entry['seq']

	# return array.
	return ref_fasta

	
def translate_file(ref_fasta, sam_file, writer, window_size):
	''' write out hit pairs.'''
	
	# make entry.
	gff_entry = writer.setup()
	
	# loop over entries.
	cnt = 0
	for entry in biobuffers.sam_buffer(sam_file):
			
		# debug.
		if cnt % 100000 == 0:
			logging.debug(".")
		cnt += 1
			
		# pull info.
		qname = entry['QNAME']
		rname = entry['RNAME']
		mirna = entry['SEQ']
			
		# check existance.
		if rname not in ref_fasta:
			continue
			
		# determine window.
		start = entry['POS'] - window_size
		stop = entry['POS'] + len(entry['QUAL']) + window_size
			
		# fix window bounds.
		if start < 0:
			start = 0
		if stop > len(ref_fasta[rname]):
			stop = len(ref_fasta[rname])
			
		# pull seqs.
		flank = ref_fasta[rname][start:stop]
		
		# populate entry.
		gff_entry['seqname'] = rname
		gff_entry['source'] = "SAM"
		gff_entry['feature'] = "hairpin?"
		gff_entry['start'] = start
		gff_entry['end'] = stop
		gff_entry['score'] = "."
		gff_entry['strand'] = "."
		gff_entry['frame'] = "."
		gff_entry['group'] = "%s;%s;%s" % (qname, mirna, flank)
		
		# write entry.
		writer.write(gff_entry)
		

############## scripts	 ##############


# read sizes of reference.
logging.info("loading reference")
ref_fasta = create_fasta(fasta_file)

# create a gff writer.
writer = biobuffers.Writer(gff_file, biobuffers.gff_dt, biobuffers.gff_fmt)

# loop over sam files.
for sam_file in sam_files:
	
	# annotate from this file.
	logging.info("annotating from %s" % sam_file)
	translate_file(ref_fasta, sam_file, writer, window_size)
	
# close output file.
writer.finalize()
