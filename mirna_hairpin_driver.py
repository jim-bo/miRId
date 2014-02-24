#!/usr/bin/python
'''
creates and executes the hairpin scripts.
'''
# imports.
import os
import sys
from string import Template
import subprocess


# parameters.
input_dir = os.path.abspath(sys.argv[1])
user_script = os.path.abspath(sys.argv[2])
output_dir = os.path.abspath(sys.argv[3])
working_dir = os.path.abspath(sys.argv[4])
log_dir = os.path.abspath(sys.argv[5])
script_dir = os.path.abspath(sys.argv[6])
bsub_script = os.path.abspath(sys.argv[7])


########## fucntions ###########

def write_script(tmp_dir, in_file, out_file, log_file, script_file):
	''' write the script '''
	
	'''
	txt = """#!/bin/bash
#BSUB -J struct
#BSUB -q normal
#BSUB -o $log_file

# make tmpdir.
if [ -d "${tmp_dir}" ]; then true
else
mkdir $tmp_dir
fi

# switch into tmp_dir.
cd ${tmp_dir}

# execute script.
python ${user_script} ${in_file} ${out_file}

"""	
	'''
	
	txt = """
	python ${user_script} ${in_file} ${out_file} ${tmp_dir}
"""
	
	# substitute.
	txt = Template(txt)
	txt = txt.substitute(log_file=log_file, tmp_dir=tmp_dir, user_script=user_script, in_file=in_file, out_file=out_file)

	# write to file.
	fout = open(script_file, "wb")
	fout.write(txt)
	fout.close()
	
	# make executable.
	subprocess.call(["chmod", "u+x", script_file])
	



########## scripts ############

# get files
files = os.listdir(input_dir)

# loop over files
idx = 0
for xfile in files:
	
	# create vars.
	tmp_dir = "%s/working_%i" % (working_dir, idx)
	input_file = "%s/%s" % (input_dir, xfile)
	output_file = "%s/output_%i.txt" % (output_dir, idx)
	script_file = "%s/script_%i.sh" % (script_dir, idx)
	log_file = "%s/log_%i.sh" % (log_dir, idx)
	
	# write script.
	write_script(tmp_dir, input_file, output_file, log_file, script_file)
	
	# submit script.
	subprocess.call([bsub_script, script_file])

	# increment idx.
	idx += 1
