#!/usr/bin/python
'''
Created on Aug 1, 2012

filters gff styled entry for complexity.

@author: eljimbo
'''
import os
import sys
import string
import subprocess
import mmap
import tempfile
import nwalign as nw

if __name__ == '__main__':
    pass

# parameters.
gff_file = sys.argv[1]
out_file = sys.argv[2]



######## functions ###########

def token_gff(line):
	''' target start stop tokenizer for GFF'''
	
	# sanity check.
	if len(line) == 0 or line[0] == "#": 
		return None
	
	# tokenize by tab.
	tmp = line.strip().split("\t")	

	if len(tmp) < 8: 
		return None
	
	# grab info.
	target = tmp[0]
	start = int(tmp[3])
	stop = int(tmp[4])
	
	# return info.
	return [target, start, stop]
	
def token_general(line):
	''' general tokenizer '''
	
	# sanity check.
	if len(line) == 0: 
		return None
	
	# tokenize by tab.
	tmp = line.strip().split("\t")	

	# return it.
	return tmp
	
def tab_gen(map_file, line_lookup, func_token):
	''' generator for tab deliminated files '''
	
	# loop over line.
	line_idx = 0
	char_idx = 0
	ecnt = 0
	line = map_file.readline()
	while line != None:
		
		# tokenize.
		res = func_token(line)
		
		# set line lookup.
		line_lookup.append(char_idx)
		
		# error check.
		if res != None:
			
			# yield it.
			res.append(line_idx)
			yield res
			ecnt = 0
			
		else:
			
			# count empties.
			ecnt += 1
			
			# sanity.
			if ecnt > 100:
				break
		
		# increment indicies.
		line_idx += 1
		char_idx += len(line)
		
		# increment line.
		line = map_file.readline()

def kmer_complexity(a, k, x):
	''' ensures a has x k-mers'''
	
	# build kmer dictionary.
	kmers = set()
	for i in range(k, len(a) + 1):
		kmers.add(a[i-k:i])
		
	# return decision.
	if len(kmers) >= x:
		return True
	else:
		return False
	
	

######## script ###########

# open gff input file.
fin_t = open(gff_file, "rb")
map_t = mmap.mmap(fin_t.fileno(), 0, access=mmap.ACCESS_READ)

# iterate over entries.
fout_master = open(out_file, "wb")
for entry in tab_gen(map_t, list(), token_general):
	
	# grab info.
	name = "test"
	tmp = entry[8].split(";")
	seq = tmp[2]
	mirna = tmp[1]

	# test complexity.
	test = kmer_complexity(mirna, 3, 8)
	
	if test == False:
		continue
		
	# write out entry.
	fout_master.write('\t'.join([str(x) for x in entry]) + "\n")

# close files.
fout_master.close()
map_t.close()
fin_t.close()


