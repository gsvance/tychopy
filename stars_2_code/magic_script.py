#!/usr/bin/env python
# magic_script.py

# This is the requested "magic" python script that reads TYCHO model files
# MAKE SURE this script can find the file tycho_model.py that I also uploaded

# Greg Vance, 4/19/17

from tycho_model import TYCHO_Model as TM  # Uses tycho_model.py
import glob  # Bash-like wildcard (*) functionality in python

# This string should match all of the model files you want to read in
tycho_files = "Ea*"  # "Ea*" for example
# This is the name of the "standard text file" to be written out
output_file = "magic.txt"

# Number of seconds in one Myr, number of erg/s in one Lsun
Myr = 31557600 * 1e6
Lsun = 3.826e33 

# Open the output file
outfile = open(output_file, "w")

# Go over each file and write out the data we want
for model in sorted(glob.glob(tycho_files)):  # Sort to put them in order
	
	# Read in the TYCHO model file using the class
	modfile = TM(model)
	
	# Extract the set of values we want
	age_Myr = modfile.header["time"] / Myr
	L_Lsun = modfile.physics["luminosity"][-1] / Lsun
	core_H = modfile.composition["p"][0]
	core_He = modfile.composition["he4"][0]
	core_D = modfile.composition["d"][0]
	
	# Format the block of text for the output file
	lines = []
	lines.append(model)
	lines.append("  age (Myr)  " + str(age_Myr))
	lines.append("  L / Lsun   " + str(L_Lsun))
	lines.append("  core H     " + str(core_H))
	lines.append("  core He    " + str(core_He))
	lines.append("  core D     " + str(core_D))
	lines.append("")  # Extra newline at the end
	
	# Write it all out to the file
	outfile.write('\n'.join(lines))
	
	# Model files are big, so free up the memory explicitly
	del modfile

# Close the output file
outfile.close()