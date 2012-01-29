#!/usr/bin/python

########################
# COPYRIGHT 2011 2012  #
# Alexander Kozik      #
# http://www.atgc.org/ #
# akozik@atgc.org      #
########################

def Lines_Extractor(in_name, out_name, tot_lines, out_lines, rand_opt):
	
	in_file  = open(in_name,  "rb")
	out_file = open(out_name, "wb")
	out_log  = open(out_name + '.log', "wb")
	
	item_array = {}
	random_entries = random.sample(range(1, tot_lines),out_lines)
	random_entries = sorted(random_entries)
	
	out_mod = out_lines/1000
	
	print " - - - - - - - - - - "
	
	k = 1
	for item in random_entries:
		item_array[item] = k
		### out_log.write(`k` + '\t' + `item` + '\n')
		if k <= 10:
			print `k` + '\t' + `item`
		if k == 11:
			print " . . . . . . . . . . "
		if k >= out_lines - 10:
			print `k` + '\t' + `item`
		k = k + 1
	
	print " - - - - - - - - - - "
	### DEBUGGING ###
	## print ""
	## print random_entries
	## print " - - - - - - - - - "
	## print item_array
	## print " - - - - - - - - - "
	## print ""
	
	l = 0
	m = 0
	f = 0
	while 1:
		t = in_file.readline()
		if t == '':
			break
		try:
			item_array[l]
			out_file.write(t)
			m = m + 1
			out_log.write(`m` + '\t' + `l` + '\n')
			print_update = math.fmod(m, out_mod)
			if print_update == 0:
				print " -= " + `m` + " lines are extracted out of " + `l` + " -= "
		except:
			f = f + 1
		l = l + 1
	
	print ""
	print " * * * * * * * * * * "
	print " * " + `m` + " lines are extracted out of " + `l` + " -= "
	print " * " + `f` + " lines are NOT extracted -= "
	print " * * * * * * * * * * "
	
	in_file.close()
	out_file.close()
	out_log.close()
	
import time
import math
import sys
import string
import random
if __name__ == "__main__":
	if len(sys.argv) <= 5 or len(sys.argv) > 6:
		print ""
		print "Program usage: "
		print "[input_file] [output_file1] [total_number_of_lines] [number_of_lines_to_extract] RANDOM"
		print ""
		print "The script will extract defined number of random entries (lines) from input file"
		print ""
		exit
	if len(sys.argv) == 6:
		in_name    = sys.argv[1]
		out_name   = sys.argv[2]
		tot_lines  = int(sys.argv[3])
		out_lines  = int(sys.argv[4])
		rand_opt   = sys.argv[5]
		if rand_opt != "RANDOM":
			print "last argument must be RANDOM"
			sys.exit()
		if in_name != out_name:
			Lines_Extractor(in_name, out_name, tot_lines, out_lines, rand_opt)
		else:
			print "Output should have different name than Input"
			sys.exit()
