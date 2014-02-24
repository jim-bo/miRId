#!/usr/bin/python
'''
Created on Apr 2, 2011

takes GFF files looking for hsit.

@author: eljimbo
'''
import os
import sys
import string
import subprocess
import mmap
import tempfile

if __name__ == '__main__':
    pass

# parameters.
gff_file = sys.argv[1]
out_file = sys.argv[2]



######## functions ###########

def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval

match_award      = 10
mismatch_penalty = -5
gap_penalty      = -15 # both for opening and extanding

def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty

def finalize(align1, align2):
    align1 = align1[::-1]    #reverse sequence 1
    align2 = align2[::-1]    #reverse sequence 2
    
    i,j = 0,0
    
    #calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:                
            symbol = symbol + align1[i]
            identity = identity + 1
            score += match_score(align1[i], align2[i])
    
        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-': 
            score += match_score(align1[i], align2[i])
            symbol += ' '
            found = 0
    
        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':          
            symbol += ' '
            score += gap_penalty
    
    identity = float(identity) / len(align1) * 100
    
    '''
    print 'Identity =', "%3.3f" % identity, 'percent'
    print 'Score =', score
    print align1
    print symbol
    print align2
    '''
    return align1, align2, score


def needle(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table
   
    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment 
    align1, align2 = '', ''
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: # end toching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1]):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    return finalize(align1, align2)
    

def water(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences
    
    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table
    pointer = zeros((m+1, n+1))    # to store the traceback path
    
    max_score = 0        # initial maximum score in DP table
    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i-1][j-1] + match_score(seq1[i-1], seq2[j-1])
            score_up = score[i][j-1] + gap_penalty
            score_left = score[i-1][j] + gap_penalty
            score[i][j] = max(0,score_left, score_up, score_diagonal)
            if score[i][j] == 0:
                pointer[i][j] = 0 # 0 means end of the path
            if score[i][j] == score_left:
                pointer[i][j] = 1 # 1 means trace up
            if score[i][j] == score_up:
                pointer[i][j] = 2 # 2 means trace left
            if score[i][j] == score_diagonal:
                pointer[i][j] = 3 # 3 means trace diagonal
            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j];
    
    align1, align2 = '', ''    # initial sequences
    
    i,j = max_i,max_j    # indices of path starting point
    
    #traceback, follow pointers
    while pointer[i][j] != 0:
        if pointer[i][j] == 3:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1
        elif pointer[i][j] == 1:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1

    return finalize(align1, align2)


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

	# create temporary files.
	temp_t = tempfile.NamedTemporaryFile(mode='w+b', delete=False)
	temp_o = tempfile.NamedTemporaryFile(mode='w+b', delete=False)
	temp_s = tempfile.NamedTemporaryFile(mode='w+b', delete=False)
	
	# write to temp file.
	temp_t.write(">%s\n%s\n" % (name, seq))
	
	# write script.
	temp_s.write("#!/bin/bash\nRNAfold -p -d2 --noLP -P vienna1.8.4.par < %s > %s" % (temp_t.name, temp_o.name))
	
	# close.
	temp_t.close()
	temp_o.close()
	temp_s.close()
	
	# run RNAfold command.
	fout = open("/dev/null", "wb")
	subprocess.call(['sh', temp_s.name], stderr=fout)
	fout.close()
	
	# read in results.
	fin = open(temp_o.name, "rb")
	idx = 0
	result = ""
	for line in fin:
		if idx == 2:
			result = line.strip().split()
			break
		idx += 1
	fin.close()
	
	# remove all files.
	os.unlink(temp_t.name)		
	os.unlink(temp_o.name)		
	os.unlink(temp_s.name)		
	
	# loop over window looking for pin.
	pinsz = 6
	matsz = 15
	matp = .5
	good_entry = False
	for i in range(pinsz, len(result[0]) - pinsz):
		
		# check pin.
		pin = result[0][i-pinsz:i+pinsz]
		bad = False
		for x in pin:
			if x != ".":
				bad = True
				break
		if bad == True:
			continue
			
		# check that % has match, left.
		mcnt = 0
		tcnt = 0
		for j in range(i-pinsz-matsz, i-pinsz):
			tcnt += 1
			if result[0][j] == "(":
				mcnt += 1
				
		# skip if not hit.
		if float(mcnt) / float(tcnt) < matp:
			continue

		# check that % has match, right.
		mcnt = 0
		tcnt = 0
		try:
			for j in range(i+pinsz, i+pinsz+matsz):
				tcnt += 1
				if result[0][j] == ")":
					mcnt += 1
		except:
			continue
			
		# skip if not hit.
		if float(mcnt) / float(tcnt) < matp:
			continue
	
		# note it.
		good_entry = True
		break
		
	# make sure mature is in right place.
	if good_entry == True:
		
		# align mature to hairpin.
		sz1 = len(seq)
		sz2 = len(mirna)
		hits = list()
		for i in range(0, sz1 - sz2 - 1):
			
			# align and score.
			#aln = nw.global_align(mirna, seq[i:i+sz2])
			aln1, aln2, sc = needle(mirna, seq[i:i+sz2])
			
			# score by dash cnt.
			#sc1 = aln[0].count("-")
			#sc2 = aln[1].count("-")
			#sc = max(sc1, sc2)
			
			# only add if 1< 
			if sc < 2:
				hits.append(i)
				
				
		# check match against location.
		for h in hits:
			
			# score by dot count.
			sc = result[0][h:h+sz2].count(".")
			
			# skip if more than 2 mismatch.
			if sc <= 2:
				
				# add alignment info to group field.
				entry[8] = "%s;%s" % (entry[8], result[0])
				
				# write it out.
				fout_master.write('\t'.join(entry[0:9]) + "\n")
				fout_master.flush()
				break


	# note it.
	#if good_entry == True:
	#	fout_master.write('\t'.join(entry[0:9]) + "\n")
	#	fout_master.flush()

# close files.
fout_master.close()
map_t.close()
fin_t.close()


