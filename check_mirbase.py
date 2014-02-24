'''
compares mirbase hairpins against a fasta file of hairpins.
'''
### imports ###
from nw_align import *
import os
import sys
import logging
import tempfile
import shutil
import subprocess
from operator import itemgetter, attrgetter
import numpy as np
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )

### parameters ###
# defaults.
MIN_IDENTITY = .98
MIN_SIZE = 20

# options.
for i in range(1, len(sys.argv)):
	try:
		if sys.argv[i] == "-min_identity":
			MIN_IDENTITY = float(sys.argv[i+1])
		elif sys.argv[i] == "-min_size":
			MIN_SIZE = int(sys.argv[i+1])			
	except:
		logging.error("bad argument: %s %s" % (str(sys.argv[i]),str(sys.argv[i+1])))
		
# required.
REF_FILE = os.path.abspath(sys.argv[-3])
QUERY_FILE = os.path.abspath(sys.argv[-2])
OUT_FILE = os.path.abspath(sys.argv[-1])

### functions ###
def load_fasta(file_path):
	''' loads fasta file into dictionary'''
	
	# read file into memory.
	fin = open(file_path)
	lines = fin.readlines()
	fin.close()
	
	# build dictionary.
	data = dict()
	seq = ""
	for line in lines:
		
		# Skip blanks.
		if len(line) < 2: continue
		if line[0] == "#": continue

		# remove blanks.
		line = line.strip()

		# Check for ids.
		if line.count(">") > 0: 

			# Check if ending seq.
			if len(seq) > 0:
		
				# save.
				data[head] = seq
				
			# reset head.
			head = line.replace(">","")
			seq = ""
			
			# skip to next line.
			continue

		# Filter chars.
		seq += line
		
	# save the last one.
	data[head] = seq
	
	# return dictionary.
	return data

def coords_gen(file_path):
	''' generator for coords tokens '''

	# read in raw data.
	fin = open(file_path,"rb")
	lines = fin.readlines()
	fin.close()
	
	# yield tokens.
	for line in lines:
		if line.count("|") == 4 and line.count("S1") == 0:
			tmp = line.strip().replace("|","").split()
			tokens = dict()
			tokens["S1"] = int(tmp[0])
			tokens["E1"] = int(tmp[1])
			tokens["S2"] = int(tmp[2])
			tokens["E2"] = int(tmp[3])
			tokens["LEN1"] = int(tmp[4])
			tokens["LEN2"] = int(tmp[5])
			tokens["IDY"] = float(tmp[6])
			tokens["RNAME"] = tmp[7]
			tokens["QNAME"] = tmp[8]
			yield tokens

def short_bad(entries):
	''' removes short and bad'''

	# loop over each.
	tmp = list()
	for entry in entries:
		
		# check length.
		if min(entry["LEN1"], entry["LEN2"]) < MIN_SIZE:
			continue
			
		# check identity.
		if entry["IDY"] < MIN_IDENTITY:
			continue
			
		# keep it.
		tmp.append(entry)
	
	return tmp

### script ###

# load my odd data.
fin = open(QUERY_FILE)
query_lines = fin.readlines()
fin.close()

# create temporary file.
fout = tempfile.NamedTemporaryFile(delete=False)
tfile = fout.name
idx = 0
for line in query_lines:
	tmp = line.strip().split()
	seq = tmp[5].replace(";","")
	fout.write(">row_%d\n%s\n" % (idx, seq))
	idx += 1
fout.close()

# run nucmer and get coords file.
entries = list()
try:
	# make temporary directory for nucmer output.
	tmp_dir = tempfile.mkdtemp()

	# create nucmer command.
	cmd = ["nucmer", "--coords", "--maxmatch", REF_FILE, tfile]
	#print ' '.join(cmd)
	#sys.exit()
	# execute nucmer command.
	try:
		fout = tempfile.TemporaryFile()
		subprocess.call(cmd, cwd=tmp_dir, stdout=fout, stderr=fout)
		fout.close()
	except:
		logging.error("error running nucmer")
		sys.exit(1)
		
	# load the resulting coord file.
	try:

		# load the entries from coords file.
		entries = [x for x in coords_gen("%s/out.coords" % tmp_dir)]	
		
	except:
		logging.error("error in nucmer output")
		sys.exit(1)
		
finally:
	try:
		shutil.rmtree(tmp_dir) # delete directory
		os.unlink(tfile)
	except:
		logging.warning("could not clean up after nucmer")
	
	
# remove short and bad entries.
entries = short_bad(entries)

# load mirbase.
mirbase = load_fasta(REF_FILE)

# track best hits for each query.
bhits = list()
for i in range(len(query_lines)):
	bhits.append([0.0, 0.0, ""])

# match mature to each candidate.
for entry in entries:
	
	# tokenize QNAME.
	idx = int(entry['QNAME'].replace("row_",""))
	line = query_lines[idx].strip().split()
	mature = line[4]
	hairpin = line[5]
	
	# extract from mirbase.
	keys = sorted([entry['S1'], entry['E1']])
	start = keys[0]
	stop = keys[1]
	ref = mirbase[entry['RNAME']][start:stop]
	
	# score mature.
	aln1, aln2 = water(ref, mature)
	bad = aln2.count("-")
	score = len(mature) - bad
	mat_idy = float(score) / float(len(mature))
	
	# score hairpin.
	aln1, aln2 = water(ref, hairpin)
	bad = aln2.count("-")
	score = len(hairpin) - bad
	har_idy = float(score) / float(len(hairpin))
	
	# check against threshold.
	if score >= MIN_SIZE and mat_idy >= MIN_IDENTITY:
		
		# check if it beats best mature idy..
		if mat_idy > bhits[idx][0]:
		
			# add it as best.
			bhits[idx][0] = mat_idy
			bhits[idx][1] = har_idy
			bhits[idx][2] = entry['RNAME']
			
		elif mat_idy == bhits[idx][0]:
			
			# check if it beats hairpin.
			if har_idy > bhits[idx][1]:
		
				# add it as best.
				bhits[idx][0] = mat_idy
				bhits[idx][1] = har_idy
				bhits[idx][2] = entry['RNAME']		
	
	
# write out the best hits.
fout = open(OUT_FILE, "wb")
for i in range(len(query_lines)):
	tokens = query_lines[i].strip().split("\t")
	if len(tokens) < 9:
		tokens.append("")
	tokens.append(bhits[i][2])
	fout.write('\t'.join(tokens) + '\n')
fout.close()
	
