# tycho_model.py

# Small python package containing a class for reading in TYCHO model files easily
# This code is ugly, I am deeply sorry, please don't undertake editing it yourself
# I, and I alone, deserve to suffer the consequences of this crime against programming
# Admittedly, not all of it is my fault... Fortran makes silly formatting decisions

# Greg Vance, 4/20/17

class TYCHO_Model:
	"""Simple class for reading and exploring TYCHO model files."""
	
	def __init__(self, filename):
		"""Read the model file of the given name."""
		
		# Save the filename, open the file
		self.filename = filename
		model = open(filename, "r")
		
		# Establish some empty data structures to fill up
		self.header = {}
		self.physics = {}
		self.nuclei= {}
		self.composition = {}
		
		# Read the data in from the file in various phases
		# The file is split into several sections that need to be kept track of
		phase = "first line"
		while phase != "done":
			
			# Save position in file, then read one line from the file
			position = model.tell()
			line = model.readline()
			
			# The very first line in the file, basic bulk info
			if phase == "first line":
				# Save the first line, extract a few useful values, and move on
				self.first_line = line.strip()
				linelist = line.split()
				if linelist[1] == "8.00new":  # Catch "8.00" and "new" abutting one another
					linelist = [linelist[0]] + ["8.00", "new"] + linelist[2:]
				self.prefix = linelist[2]
				self.mass = float(linelist[3])
				self.network_size = int(linelist[4])
				self.metallicity = float(linelist[5])
				phase = "header"
			
			# The header of data before the bulk of the file starts
			elif phase == "header":
				# Check for end of the header phase
				if line.strip() == "":
					phase = "physics"
					continue
				# Split the variable name from the value, convert if possible
				magic_split = 32  # The 'column split' is about here
				variable = line[:magic_split].strip()
				value_str = line[magic_split:].strip()
				try:
					value = int(value_str)
				except ValueError:
					try:
						value = float(value_str)
					except ValueError:
						value = value_str
				self.header[variable] = value
			
			# First section of physical zone data before nuclear information
			elif phase == "physics":
				# Check for blank lines at end of header
				if line.strip() == "":
					continue
				# Check for ints at the end of this phase
				try:
					test = int(line.split()[0])
					phase = "nuclei"
					model.seek(position)  # Rewind the file one line
					continue
				except ValueError:
					pass
				# Read in the varaible name, set up the data list
				variable = line.strip()
				self.physics[variable] = []
				# Read all of the float values
				while len(self.physics[variable]) < self.header["kk"] + 2:
					linelist = list(reversed(model.readline().split()))
					while len(linelist) > 0:
						value = whyyy_float(linelist.pop())
						self.physics[variable].append(value)
			
			# Short section in the middle with nuclei and initial composition
			elif phase == "nuclei":
				# Rewind the file and read lines for the four quantities tabulated here
				model.seek(position)
				variables = ["nz", "nn", "name", "initial mass fraction"]
				for variable in variables:
					self.nuclei[variable] = []
					while len(self.nuclei[variable]) < self.header["netsize"] + 1:
						linelist = list(reversed(model.readline().split()))
						while len(linelist) > 0:
							if variable == "nz" or variable == "nn":
								value = int(linelist.pop())
							elif variable == "name":
								value = linelist.pop()
								# Catch the stupid long isotope names abutting one another
								while len(value) > 5:
									end = value[-5:]
									linelist.append(end)
									value = value[:-5]
							elif variable == "initial mass fraction":
								value = whyyy_float(linelist.pop())
							else:
								raise ValueError("bad variable in nuclei read")
							self.nuclei[variable].append(value)
				# Move on to the next phase of reading
				phase = "composition"
			
			# Most of the file is made up of this last section with composition data
			elif phase == "composition":
				# Check for blank lines signalling the end of the file
				if line.strip() == "":
					phase = "done"
					continue
				# Read in the isotope name, set up the data list
				isotope = line.strip()
				self.composition[isotope] = []
				# Read all of the float values
				while len(self.composition[isotope]) < self.header["kk"] + 2:
					linelist = list(reversed(model.readline().split()))
					while len(linelist) > 0:
						value = whyyy_float(linelist.pop())
						self.composition[isotope].append(value)
			
			# Something is wrong, raise an error
			else:
				raise ValueError("read phase not identified")
		
		# Close the file
		model.close()

def whyyy_float(string):
	"""Parse stupid floats with 'E' missing that Fortran has butchered."""
	try:
		value = float(string)
	except ValueError:
		i = -1
		while string[i].isdigit():
			i -= 1
		man = string[:i]
		exp = string[i:]
		value = float(man + 'E' + exp)
	return value

def main():
	"""Test function to make sure the class works properly."""
	tycho = TYCHO_Model("Ea00650")
	print "File:", tycho.filename
	print tycho.first_line
	print "Prefix:", tycho.prefix
	print "Mass:", tycho.mass
	print "NNet:", tycho.network_size
	print "Z:", tycho.metallicity
	print "kk:", tycho.header["kk"]
	print "xm(kk):", tycho.header["xm(kk)"]
	print "opatab:", tycho.header["opacity tables used"]
	print "core H:", tycho.composition["p"][0]

if __name__ == "__main__":
	main()
