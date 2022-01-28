#!/usr/bin/env python
# Print out the mass of TYCHO model files
# Greg Vance, 4/23/17

from tycho_model import TYCHO_Model as TM
import sys

msun = 1.98855e33 # grams

for mod in sys.argv[1:]:
	model = TM(mod)
	mass = sum(model.physics["zone mass"])
	print mod, "\t", mass/msun, "Msun"
	del model

