'''
takes reference and assembled contigs and creates scaffolding
test case. You can specify the following parameters:
min_identity: alignment must be equal to or greater than this value(0..1).
min_size: alignment must be equal to or greater than this size(0...).
overlap: number of non-overlapping bases.
norepeats: boolean disallow query to be found more than once.
'''
### imports ###
from data_structs.agps import save_agps
from data_structs.types import agp_dt

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
MIN_IDENTITY = .97
MIN_SIZE = 200
OVERLAP = sys.maxint
REPEATS = True

# options.
for i in range(1, len(sys.argv)):
	try:
		if sys.argv[i] == "-min_identity":
			MIN_IDENTITY = float(sys.argv[i+1])
		elif sys.argv[i] == "-min_size":
			MIN_SIZE = int(sys.argv[i+1])			
		elif sys.argv[i] == "-repeat_off":
			REPEAT_OFF = float(sys.argv[i+1])	
		elif sys.argv[i] == "-overlap":
			OVERLAP = int(sys.argv[i+1])	
		elif sys.argv[i] == "-repeats":
			REPEATS = False	
	except:
		logging.error("bad argument: %s %s" % (str(sys.argv[i]),str(sys.argv[i+1])))
		
# required.
REF_FILE = os.path.abspath(sys.argv[-6])
QUERY_FILE = os.path.abspath(sys.argv[-5])
TEST_FILE = os.path.abspath(sys.argv[-4])
AGP_FILE = os.path.abspath(sys.argv[-3])
INSERT_SIZE = int(sys.argv[-2])
STD_DEV = int(sys.argv[-1])

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

def overlapping(entries):
	''' removes overlapping '''

	# build list of (start,stop, entries index) for each ref.
	refs = dict()
	for i in range(len(entries)):
		
		# build tuple.
		tpl = (entries[i]["S1"], entries[i]["E1"], i)
		
		# add to list.
		if entries[i]["RNAME"] not in refs:
			refs[entries[i]["RNAME"]] = list()
		refs[entries[i]["RNAME"]].append(tpl)
			
	# remove overlapping.
	tmp = list()
	for rname in refs:
		
		# sort.
		active = sorted(refs[rname], key=itemgetter(0,1))
		
		# scan and look for overlapping.
		i = 0
		sz = len(active) - 1
		while i < sz:
			
			# compute overlap.
			keys1 = sorted(active[i][0:2])
			keys2 = sorted(active[i+1][0:2])
			olap = -1 * (keys2[0] - keys1[1])

			# check pair for overlap.
			if olap > OVERLAP:
			
				# delete i+1.
				active.pop(i+1)
				
				# reset counter.
				i = 0
				sz = len(active) - 1
				
			else:	
				# increment normally.
				i += 1
			
		# add back active.
		for x in active:
			tmp.append(entries[x[2]])

	return tmp

def repeats(entries):
	''' removes repeats '''

	# build list of (start,stop, entries index) for each query.
	refs = dict()
	for i in range(len(entries)):
		
		# build tuple.
		tpl = (entries[i]["S2"], entries[i]["E2"], i)
		
		# add to list.
		if entries[i]["QNAME"] not in refs:
			refs[entries[i]["QNAME"]] = list()
		refs[entries[i]["QNAME"]].append(tpl)
			
	# remove overlapping.
	tmp = list()
	for rname in refs:
		
		# sort.
		active = sorted(refs[rname], key=itemgetter(0,1))
		
		# scan and look for overlapping.
		i = 0
		sz = len(active) - 1
		while i < sz:
			
			# compute overlap.
			keys1 = sorted(active[i][0:2])
			keys2 = sorted(active[i+1][0:2])
			olap = -1 * (keys2[0] - keys1[1])

			# check pair for overlap.
			if olap > 0:
			
				# delete i+1.
				active.pop(i+1)
				
				# reset counter.
				i = 0
				sz = len(active) - 1
				
			else:	
				# increment normally.
				i += 1
			
		# add back active.
		for x in active:
			tmp.append(entries[x[2]])

	return tmp
	

def strip_gaps(agps):
	''' removes gaps '''
	
	# build scaffold subsets ignoring existing gaps.
	scafs = dict()
	for i in range(agps.shape[0]):
		if agps[i]['scaf_name'] not in scafs:
			scafs[agps[i]['scaf_name']] = list()
			
		if  agps[i]['comp_type'] == "W":
			scafs[agps[i]['scaf_name']].append(i)
		
		
	# count number without gaps.
	new_size = 0
	for skey in scafs:
		new_size += len(scafs[skey])
		
	# create new agps.
	new_agps = np.zeros(new_size, dtype=agp_dt)
	
	# copy contigs and add gaps.
	idx = 0
	for skey in scafs:
		
		scaf_idx = 1
		for i in range(len(scafs[skey])):
			
			j = scafs[skey][i]
			
			# copy existing.
			new_agps[idx] = agps[j]
			new_agps[idx]['scaf_idx'] = scaf_idx
			scaf_idx += 1
			idx += 1
			
	# return them modified agps.
	return new_agps

def add_gaps(agps):
	''' adds gaps to  agp file '''
		
	# build scaffold subsets ignoring existing gaps.
	scafs = dict()
	for i in range(agps.shape[0]):
		if agps[i]['scaf_name'] not in scafs:
			scafs[agps[i]['scaf_name']] = list()
			
		if  agps[i]['comp_type'] == "W":
			scafs[agps[i]['scaf_name']].append(i)
		
	# count number of missing gaps.
	new_size = 0
	for skey in scafs:
		new_size += len(scafs[skey])
		new_size += len(scafs[skey]) - 1
		
	# create new agps.
	new_agps = np.zeros(new_size, dtype=agp_dt)
	
	# fill out new agp.
	idx = 0
	for scaf in scafs.values():
		
		# zero index.
		scaf_idx = 1
		
		# copy data.
		for i in range(len(scaf)):
			
			# copy the contig.
			new_agps[idx] = agps[scaf[i]]
			idx += 1
			scaf_idx += 1
			
			# check if we add gap.
			if scaf[i] != scaf[-1]:
			
				# simplify.
				gstart = new_agps[idx - 1]['scaf_stop'] + 1
				gstop = agps[scaf[i + 1]]['scaf_start'] - 1
			
				# sanity check.
				if gstart >= gstop:
					
					if abs(gstart - gstop) > 10:
						print "BAD JUJU"
						sys.exit(1)
					
					# move up the start of next one.
					agps[scaf[i + 1]]['scaf_start'] = agps[scaf[i + 1]]['scaf_start'] + abs(gstart - gstop) + 2
				
				# add gap.
				new_agps[idx]['scaf_name'] = new_agps[idx - 1]['scaf_name']
				new_agps[idx]['scaf_idx'] = scaf_idx
				new_agps[idx]['scaf_start'] = new_agps[idx - 1]['scaf_stop'] + 1
				new_agps[idx]['scaf_stop'] = agps[scaf[i + 1]]['scaf_start'] - 1
				new_agps[idx]['comp_type'] = "N"
				new_agps[idx]['comp_name'] = "gap"
				new_agps[idx]['comp_start'] = 1 
				new_agps[idx]['comp_stop'] = new_agps[idx]['scaf_stop'] - new_agps[idx]['scaf_start']
				new_agps[idx]['comp_orien'] = 0
				
				if new_agps[idx]['scaf_start'] > new_agps[idx]['scaf_stop'] :
					print "BAD JUJU2"
					print new_agps[idx - 1]
					print new_agps[idx]
					print agps[scaf[i + 1]]
					sys.exit(1)
				
				idx += 1
				scaf_idx += 1
				
	# return them modified agps.
	return new_agps
			

def gap_split(agps):
	''' splits scaffolds at large gaps'''
	
	# build scaffold subsets ignoring existing gaps.
	scafs = dict()
	for i in range(agps.shape[0]):
		if agps[i]['scaf_name'] not in scafs:
			scafs[agps[i]['scaf_name']] = list()
			
		if  agps[i]['comp_type'] == "W":
			scafs[agps[i]['scaf_name']].append(i)
		
			
	# loop over scaffolds.
	new_scafs = list()
	for skey in scafs:
		
		# bootstrap list.
		cur = list()
		cur.append(scafs[skey][0])
		
		# loop till end.
		for i in range(1, len(scafs[skey])):
			
			# simplify.
			gap = agps[scafs[skey][i]]['scaf_start'] - agps[scafs[skey][i-1]]['scaf_stop']
			
			# check if we add to finalize.
			if gap > INSERT_SIZE + (5.6 * STD_DEV) and len(cur) > 0:
				
				# finalize.
				new_scafs.append(cur)
					
				# reset.
				cur = list()
				
			# just add.
			cur.append(scafs[skey][i])

		# finalize.
		if len(cur) > 0:
			new_scafs.append(cur)
			

	# count size of new agps.
	sz = 0
	for a in new_scafs:
		sz += len(a)
			
	# create new agp array.
	new_agps = np.zeros(sz, dtype=agp_dt)
			
	# build new scaffolds...
	scaf_count = 1
	idx = 0
	for scaf in new_scafs:
		
		# reset trackers.
		scaf_name = "scaffold_%i" % scaf_count
		scaf_idx = 1
		scaf_offset = agps[scaf[0]]['scaf_start']
		scaf_count += 1
		
		# start copying to new array.
		for i in scaf:
			
			# do whole copy.
			new_agps[idx]['scaf_name'] = scaf_name
			new_agps[idx]['scaf_idx'] = scaf_idx
			new_agps[idx]['scaf_start'] = agps[i]['scaf_start'] - scaf_offset + 1
			new_agps[idx]['scaf_stop'] = agps[i]['scaf_stop'] - scaf_offset + 1
			new_agps[idx]['comp_type'] = "W"
			new_agps[idx]['comp_name'] = agps[i]['comp_name']
			new_agps[idx]['comp_start'] = agps[i]['comp_start']
			new_agps[idx]['comp_stop'] = agps[i]['comp_stop']
			new_agps[idx]['comp_orien'] = agps[i]['comp_orien']
			
			# update index.
			scaf_idx += 1
			idx += 1
		
	# return new agps.
	return new_agps
	
### script ###
'''
# run nucmer and get coords file.
coord_lines = list()
try:
	# make temporary directory for nucmer output.
	tmp_dir = tempfile.mkdtemp()

	# create nucmer command.
	cmd = ["nucmer", "--coords", "--maxmatch", REF_FILE, QUERY_FILE]
	
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
		fin = open("%s/out.coords" % tmp_dir, "r")
		for line in fin:
			coord_lines.append(line)
		fin.close()
	except:
		logging.error("error in nucmer output")
		sys.exit(1)
		
finally:
	try:
		shutil.rmtree(tmp_dir) # delete directory
	except:
		logging.warning("could not clean up after nucmer")
'''

# load the entries from coords file.
entries = [x for x in coords_gen("out.coords")]
	
# remove short and bad entries.
entries = short_bad(entries)

# remove overlapping.
entries = overlapping(entries)

# remove repeats.
if REPEATS == False:
	entries = repeats(entries)

# load the reference and test into memory.
ref = load_fasta(REF_FILE)
test = load_fasta(QUERY_FILE)

# write out a fasta file.
i = 0
fout = open(TEST_FILE, "w")
for entry in entries:
	
	# sort key.
	key = sorted([entry["S2"], entry["E2"]])
	
	# get fasta info.
	name = "contig_%i" % i
	seq = test[entry["QNAME"]][key[0]:key[1]]
	i += 1
	
	# write it.
	fout.write(">%s\n%s\n" % (name, seq))
fout.close()

# group by reference.
refs = dict()
for i in range(len(entries)):
		
	# add to list.
	if entries[i]["RNAME"] not in refs:
		refs[entries[i]["RNAME"]] = list()
	refs[entries[i]["RNAME"]].append(entries[i])
	
# sort by start.
for rname in refs:
	refs[rname].sort(key=itemgetter('S1'))

# create the AGP structure.
agps = np.zeros(len(entries), dtype=agp_dt)

# copy data into it.
idx = 0
idx_cnt = dict()
ctg_cnt = dict()
for rname in refs:
	
	# get idx cnt.
	if entry['RNAME'] not in idx_cnt:
		idx_cnt[entry['RNAME']] = 1
	scaf_idx = idx_cnt[entry['RNAME']]
	idx_cnt[entry['RNAME']] += 1
	
	# loop over each entry.
	for entry in refs[rname]:
	
		# get orien.
		if entry['S2'] > entry['E2']:
			orien = 1
		else:
			orien = 0
		
		# copy to data.
		agps[idx]['scaf_name'] = entry['RNAME']
		agps[idx]['scaf_idx'] = scaf_idx
		agps[idx]['scaf_start'] = entry['S1']
		agps[idx]['scaf_stop'] = entry['E1']
		agps[idx]['comp_type'] = "W"
		agps[idx]['comp_name'] = "contig_%i" % idx
		agps[idx]['comp_start'] = 1 
		agps[idx]['comp_stop'] = agps[idx]['scaf_stop'] - agps[idx]['scaf_start']
		agps[idx]['comp_orien'] = orien
		idx += 1
	
# split at large gaps.
agps = gap_split(agps)

print agps
sys.exit()
# add the gaps back in.
agps = add_gaps(agps)


