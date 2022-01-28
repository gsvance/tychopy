# Greg Vance, 3/28/17
# Takes a list of TYCHO model files from the command line
# Spits out the time in seconds from each file's header

import sys

for arg in sys.argv[1:]:
	with open(arg, 'r') as f:
		for line in f:
			if line.find('time') > 0:
				break
	l = line.strip().split()
	print arg + ":", l[1]

