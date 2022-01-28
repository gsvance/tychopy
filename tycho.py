# tycho.py

# Behold the *new* version of the code I originally wrote for Stars II in 2017!
# It has been updated for Python 3 and now uses the re package for parsing
# The main thing here is a class that reads TYCHO model files into Numpy arrays
# The code is now *significantly* less ugly and bootleg than the original was

# Last modified 20 Jul 2021 by Greg Vance


# TODO:
# Get dimensions for: rotation, angular momentum, convection speed, luminosity
# Is zone mass the just the mass coordinate value? Or something else?
# Figure out how to calculate internal energy


# IMPORT STATEMENTS
###################

from __future__ import division, print_function  # For Python 2 compatibility

import re
import numpy as np
from collections import OrderedDict
import astropy.units as u


# GLOBAL CONSTANTS
##################

# Regular expression patterns for matching ints, floats, isotopes, etc.
INT_PAT = r"[-+]?[0-9]+"  # Integer with possible sign
FLT_PAT = r"[-+]?(?:[0-9]+\.[0-9]*|\.[0-9]+)" \
	+ r"(?:[eE][-+]?[0-9]+|[-+][0-9]+)?"  # Allow for missing 'E'
ISO_PAT = r"[a-zY]{1,2}[0-9]{0,3}"  # Name of isotope, e.g., d, he3, ti44, Ye
STR_PAT = r"\S+(?: \S+)*"  # String with optional single spaces

# Regular expression patterns for matching blocks of data in a Tycho model file
FIRST_LINE = r"^TYCHO( +\S+)+ *$\n?"
HEADER_LINES = r"(^ *" + STR_PAT + r" {2,}" + STR_PAT + r" *$\n?){2,}"
BLANK_LINES = r"(^ *$\n?){2,}"
LABEL_LINE = r"^ *" + STR_PAT + r" *$\n?"
FLOAT_BLOCK = r"(^ *" + FLT_PAT + r"( +" + FLT_PAT + r")* *$\n?){2,}"
INT_BLOCK = r"(^ *" + INT_PAT + r"( +" + INT_PAT + r")* *$\n?){2,}"
ISOTOPE_BLOCK = r"(^ *" + ISO_PAT + r"( *" + ISO_PAT + r")* *$\n?){2,}"

# Combined pattern for matching and classifying anything in a Tycho model file
TYCHO_PATTERN = "|".join([
	r"(?P<FIRST>" + FIRST_LINE + r")",
	r"(?P<ISOTOPE>" + ISOTOPE_BLOCK + r")",
	r"(?P<HEADER>" + HEADER_LINES + r")",
	r"(?P<INT>" + INT_BLOCK + r")",
	r"(?P<FLOAT>" + FLOAT_BLOCK + r")",
	r"(?P<LABEL>" + LABEL_LINE + r")",
	r"(?P<BLANK>" + BLANK_LINES + r")",
])

# Units table for trying to keep track of how things are measured in Tycho.
# Patrick says everything is in CGS, but I'm not 100% sure about the dimensions
# of every quantity that might occur in here.
UNITS_TABLE = {
	"radius": u.centimeter,
	"velocity" : u.centimeter / u.second,
	"rotation" : None,  # Probably rad/s, not sure
	"angular momentum": None,  # Not sure, ask Patrick
	"temperature": u.Kelvin,
	"specific volume": u.centimeter**3 / u.gram,  # == 1/rho
	"pressure": u.gram / (u.centimeter * u.second**2),  # force per area
	"zone mass": None,  # Almost certainly grams, not 100% sure
	"convection speed": None,  # Probably cm/s, not sure
	"luminosity": None,  # Probably erg/s, not sure
	"nz": u.Dalton,  # Not really necessary...
 	"nn": u.Dalton,  # Not really necessary...
	"initial composition": u.dimensionless_unscaled,  # Mass fractions
}


# UTILITY FUNCTIONS
###################

def recover_float(string):
	"""Parse and return a floating point number from a string that Fortran may
	or may not have mangled.
	"""
	
	# Remove whitespace and try to parse the float the easy way
	string = string.strip()
	try:
		value = float(string)
	
	# If that strategy fails, then check that the 'E' is missing and fix it
	except ValueError:
		match = re.match(r"^([-+]?[0-9]+\.[0-9]+)([-+][0-9]+)$", string)
		if match:
			value = float(match.group(1) + "E" + match.group(2))
		else:
			msg = "recover_float failed to parse string: " + repr(string)
			raise ValueError(msg)
	
	return value

def convert(string):
	"""Convert the given string to an int or float if appropriate."""
	
	# Try converting to an int first, then to a float, then just give up
	try:
		value = int(string)
	except ValueError:
		try:
			value = recover_float(string)  # Never trust floats from Fortran
		except ValueError:
			value = str(string)  # Eliminate any weird subclasses of str
	
	return value


# MAIN TYCHOMODEL CLASS
#######################

class TychoModel(OrderedDict):
	"""An OrderedDict subclass for reading Tycho model files into an organized
	series of labeled Numpy arrays for easy data exploration and processing.
	The file's header data is stored in the secondary OrderedDict self.header,
	and all other data is recorded in self.items(). The object is mutable in
	the same way that an OrderedDict is mutable, but changing the data does not
	change the contents of the Tycho model file it was read from.
	"""
	
	def __init__(self, filename):
		"""Initialize a new TychoModel object by reading data from the
		specified model file.
		"""
		
		# Initialize self as an OrderedDict and set up the header OrderedDict
		super(TychoModel, self).__init__()
		self.header = OrderedDict()
		
		# Save the name of the Tycho model file and read its contents
		self.filename = filename
		with open(self.filename, "r") as modelfile:
			modelstring = modelfile.read()
		
		# Use the Tycho re pattern to classify all block of data in the file
		matches = re.finditer(TYCHO_PATTERN, modelstring, flags=re.MULTILINE)
		
		# Loop over the data blocks and process each of them depending on type
		label_without_data = None  # Record any label waiting for data
		for match in matches:
			
			# Extract the text of the match and how it was classified
			text = match.group()
			classified = match.lastgroup
			
			# If the last thing was a label, then this had better be floats
			if label_without_data is not None and classified != "FLOAT":
				msg = "was expecting FLOAT data for label, but instead got " \
					+ classified
				raise SyntaxError(msg)
			
			# Choose our parsing approach depending on the data classification
			if classified == "FIRST":
				self._parse_first(text)
			elif classified == "HEADER":
				self._parse_header(text)
			elif classified == "LABEL":
				label_without_data = text.strip()
			elif classified == "FLOAT":
				if label_without_data is not None:
					self._parse_float(text, label_without_data)
					label_without_data = None
				elif "initial composition" not in self:
					self._parse_float(text, "initial composition")
				else:
					msg = "found a secondary unlabeled block of floats"
					raise SyntaxError(msg)
			elif classified == "INT":
				if ("nz" not in self) and ("nn" not in self):
					self._parse_int(text, ["nz", "nn"])
				else:
					msg = "found an additional block of unlabeled ints"
					raise SyntaxError(msg)
			elif classified == "ISOTOPE":
				if "isotope" not in self:
					self._parse_isotope(text, "isotope")
				else:
					msg = "found a secondary block of isotope names"
					raise SyntaxError(msg)
			elif classified == "BLANK":
				pass  # Recognize blocks of blank lines, but ignore them
			else:
				msg = "unexpected data block classification: " \
					+ str(classified)
				raise ValueError(msg)
		
		# Make sure there isn't a remaining label without data
		if label_without_data is not None:
			msg = "model file ended while waiting for incoming data"
			raise SyntaxError(msg)
	
	def _parse_first(self, string):
		"""Parse the first line of the model file."""
		
		# Split the first line string into bit at any amount of whitespace
		# WEIRD EDGE CASE: the third thing on this line is the simulation file
		# prefix, which is usually just two characters. However, if the model
		# file is an imodel that was just generated, the prefix will be the
		# three-character string "new", which will abut the version number (the
		# second thing on this line). In this situation, we need to add some
		# whitespace in between before splitting the string. I originally tried
		# to solve this by using a zero-length splitter in re.split, but that
		# raises a FutureWarning in Python 3 apparently.
		string = re.sub(r"([0-9]+\.[0-9]+)(new )", r"\1 \2", string.strip())
		bits = re.split(r" +", string.strip())
		
		for i, bit in enumerate(bits):
			
			# Convert the bit to an int or a float if appropriate
			# Store the bits in the header dict by their ordered index
			self.header[i] = convert(bit)
	
	def _parse_header(self, string):
		"""Parse the header of the model file, not including the first line."""
		
		# Extract the variable name and value from each header line
		pattern = r"^ *(" + STR_PAT + r") {2,}(" + STR_PAT + r") *$"
		matches = re.finditer(pattern, string.strip(), flags=re.MULTILINE)
		
		for match in matches:
			
			name, value = match.group(1, 2)
			
			# Convert the value to an int or a float if appropriate
			# Save value to the header dict using the variable name as the key
			self.header[str(name)] = convert(value)
	
	def _parse_float(self, string, label):
		"""Parse a multiline block of floating point values and store the final
		Numpy array under self[label].
		"""
		
		# Use re to extract everything that looks like a float
		matches = re.finditer(FLT_PAT, string.strip(), flags=re.MULTILINE)
		floats = [match.group() for match in matches]
		
		# Index using the given label, eliminate any weird str subclasses
		self[str(label)] = np.array([recover_float(f) for f in floats])
	
	def _parse_int(self, string, labels):
		"""Parse a multiline block of integer values (possibly made up of 2 or
		more unlabeled blocks that have been conflated together), and store the
		final Numpy array(s) under self[labels[0]], self[labels[1]], etc.
		"""
		
		# Use re to extract everything that looks like an int
		matches = re.finditer(INT_PAT, string.strip(), flags=re.MULTILINE)
		ints = [match.group() for match in matches]
		
		# We need to split the block into proton numbers and neutron numbers
		# Sanity check: Make sure the number of ints is evenly divisible
		if len(ints) % len(labels) != 0:
			msg = "int block failed to split apart evenly"
			raise SyntaxError(msg)
		
		size = len(ints) // len(labels)
		for i, label in enumerate(labels):
			
			# Where to slice the list to split it properly
			start = i * size
			stop = (i + 1) * size
			
			# Convert sliced list of strings to a Numpy array of ints
			# Index using the given labels, eliminate any weird str subclasses
			self[str(label)] = np.array([int(j) for j in ints[start:stop]])
	
	def _parse_isotope(self, string, label):
		"""Parse a multiline block of isotope names (possibly abutting one
		another) and store the final Numpy array under self[label].
		"""
		
		# Use re to extract everything that looks like an isotope name
		matches = re.finditer(ISO_PAT, string.strip(), flags=re.MULTILINE)
		isotopes = [match.group() for match in matches]
		
		# Convert list of strings to a Numpy array of Unicode strings
		# Index using the given label, eliminate any weird str subclasses
		self[str(label)] = np.array([str(iso) for iso in isotopes], dtype="U")
	
	def __repr__(self):
		return "TychoModel(" + repr(self.filename) + ")"
	
	def __str__(self):
		return repr(self)
	
	def _repr_pretty_(self, p, cycle):
		"""IPython "pretty printer" override method to enforce __repr__."""
		
		# Without defining this extra method, the IPython terminal defers to
		# its default repr formatting for the OrderedDict class and flat-out
		# refuses to use the __repr__ method that I define above. When it does
		# that, it prints the entire contents of the model file including every
		# float in every Numpy array in a manner that is both excessive and
		# unhelpful. It also eats up a few seconds just to get it done.
		
		p.text(repr(self))
	
	def get_units(self, quantity):
		"""Given a key from the OrderedDict self, try to return appropriate
		units for that quantity using the astropy.units module.
		"""
		
		# Make sure the quantity exists in this file
		if quantity not in self.keys():
			msg = "get_units quantity " + repr(quantity) + " is not in file"
			raise KeyError(msg)
		
		# Try and look up the appropriate units in the table
		# If the quantity is an isotope or Ye, then it's dimensionless
		try:
			lookup = UNITS_TABLE[quantity]
		except KeyError:
			if quantity in self["isotope"]:
				return u.dimensionless_unscaled  # Mass fractions and Ye
			else:
				msg = "quantity " + repr(quantity) + \
					" not available for units lookup"
				raise KeyError(msg)
		
		# Return the lookup value only if it is not None
		if lookup is not None:
			return lookup
		else:
			msg = "sorry, the units of " + repr(quantity) + " are uncertain"
			raise KeyError(msg)
	
	def add_density(self):
		"""Add an entry for density to the OrderedDict self by calculating the
		multiplicative inverse of "specific volume."
		"""
		
		try:
			self["density"] = 1.0 / self["specific volume"]
		except KeyError:
			msg = "cannot calculate density, no specific volume in model file"
			raise KeyError(msg)
	
	def add_internal_energy(self):
		""""""
		
		raise NotImplemented


# TEST CODE
###########

def test():
	"""Test function to make sure the class works correctly."""
	
	#test_file = "../lanl/2021_summer/research/1d_models/from_patrick/1522334"
	test_file = "../lanl/2021_summer/research/1d_models/from_patrick/2018918"
	model = TychoModel(test_file)
	
	print("File:", model.filename)
	print([model.header[i] for i in model.header if type(i) is type(0)])
	print("Prefix:", model.header[2])
	print("Mass:", model.header[3])
	print("NNet:", model.header[4])
	print("Z:", model.header[5])
	print("kk:", model.header["kk"])
	print("xm(kk):", model.header["xm(kk)"])
	print("opatab:", model.header["opacity tables used"])
	print("core H:", model["p"][0])

# Only run the test code if this file is executed directly
if __name__ == "__main__":
	test()

