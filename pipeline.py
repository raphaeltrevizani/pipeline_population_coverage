#!/usr/bin/python

import subprocess
from os.path import join, exists, getsize
from os import makedirs, rename
from shutil import rmtree

# -------------------------------------	
def parse_parameters(inputfile):

	'''
		Reads input file containing parameters (default: parameters.txt)
		and returns a dictionary {variables:values} as specified in
		the input file as variables = values
	'''

	parameters = dict()

	# Function that detects if a line contains alphanumeric characters
	with open(inputfile) as par:
		for line in par:
		
			# Removes non-alphanumeric chars
			line = line.strip('_\t\n\'\"').replace(' ', '').replace('-', '')

			# Skip comment lines
			if line.startswith('#'):
				continue
	
			# If line contains alphanumeric chars, is valid
			if line: 

				# File format: variable = value
				k,v = line.split('=') 

				# Removes non-alphanumeric extra characters from key
				k = k.strip('_\t\n\'\"').replace(' ', '').replace('-', '').lower()
				v = v.strip('_\t\n\'\"').replace(' ', '').replace('-', '')

				# Adds to dictionary of parameters
				parameters[k] = v


	return parameters

# -------------------------------------
def parse_hla_i_file(parameters):
	
	''' 
		Parses the HLA-I (default: IEDB reference file, data/hla_ref_set.class_i.txt)
		Returns a string of comma-separated HLAs and a string of comma separated sizes
		as required by the IEDB prediction tool
	'''
	
	hla_file = parameters['hlaifile'] 
	mhci_sizes = parameters['mhcisizes']
	
	# Store the results
	hla_list = list()
	len_list = list()

	valid_sizes = [int(size) for size in mhci_sizes.split(',')]

	with open(hla_file) as hf:

		for line in hf:
			line = line.rstrip()
			hla, size = line.split(',')

			if int(size) in valid_sizes:
				hla_list.append(hla)
				len_list.append(size)

	return ','.join(hla_list), ','.join(len_list)

# -------------------------------------
def parse_hla_ii_file(parameters):
	
	''' 
		Parses the HLA-I (default: IEDB reference file, data/hla_ref_set.class_i.txt)
		Returns a string of comma-separated HLAs and a string of comma separated sizes
		as required by the IEDB prediction tool
	'''
	
	hla_file = parameters['hlaiifile'] 
	mhcii_sizes = parameters['mhciisizes']
	
	with open(hla_file) as hf:
		hla_list = hf.readlines()

	hla_list = [hla.rstrip() for hla in hla_list]
	hla_list = ','.join(hla_list)
	
	return hla_list, mhcii_sizes

# -------------------------------------
def parse_hla_file(parameters, mhc_class):

	'''
		Parses IEDB reference HLA file for 
		MHC-I or MHC-II
	'''

	if mhc_class.lower() == 'i':
		return parse_hla_i_file(parameters)
	elif mhc_class.lower() == 'ii':
		return parse_hla_ii_file(parameters)
	else:
		print(f'Unknown MHC class: {mhc_class}.')
# -------------------------------------
def get_alleles_and_binders(prediction, parameters, mhc_class):

	'''
		Isolates from IEDB prediction output 
		the tuples of peptides and corresponding alleles
		that are above the binding threshold
	'''

	mhc_class = mhc_class.lower()

	# Decides whether it is MHC-I or MHC-II
	if mhc_class == 'i':
		peptide_column = 5
		
	if mhc_class == 'ii':
		peptide_column = 6
	
	cutoff = float(parameters['mhc' + mhc_class + 'threshold'])
	

	# Splits str(prediction) into a list of rows; ignores header; filters empty
	pred_lines = filter(None, prediction.split('\n')[1:])

	# Dictionary format {peptide:str(allele1,...,alleleN)}
	pred_results = dict()
	
	for line in pred_lines:

		# Parses tab-separated output
		fields = line.split('\t')
		allele = fields[0]
		peptide = fields[peptide_column]
		score = float(fields[-2])

		# Remove alpha chain in heterodimers of MHC-II
		if mhc_class == 'ii' and '/' in allele:
			allele = 'HLA-' + allele.split('/')[1]

		# Selects lines below the binding threshold
		if score <= cutoff:
			if peptide in pred_results:
				pred_results[peptide] += ',' + allele
			else:
				pred_results[peptide] = allele

	return pred_results

# -------------------------------------
def filter_locus(alleles, locus):
	''' From the alleles list (eg: HLA-A*02:01), keeps only those that pertain to arg:locus (eg: A)'''
	return ','.join([hla for hla in alleles.split(',') if locus in hla.split('*')[0].split('-')[1]])
# -------------------------------------
def create_pop_coverage_input(pred_results, parameters, filename, locus=False): 

	'''
		Creates the input file needed to call 
		the population coverage tool
		pred_results: dictionary {str(epitope):str('hla1,hla2')}
		Returns: created file name
	'''
	
	# File with population coverage input in tmpdir/:
	pop_input_name = join(parameters['temporarydirectory'], filename)

	# Writes to file
	with open(pop_input_name, 'a+') as f:
		for peptide, allele_list in pred_results.items():
			if locus:
				allele_list = filter_locus(allele_list, locus)
			if allele_list:
				line = peptide + '\t' + allele_list + '\n'
				f.write(line)

	# Returns file name
	return pop_input_name

# -------------------------------------
def run_mhc_prediction(inputdata, parameters, mhc_class):
	'''
		Uses the API to predict binders of arg:mhc_class for each 
		arg:inputdata
	'''

	# Ensures the mhc class is not capitalized
	mhc_class = mhc_class.lower()
	
	# Gets HLAs and respective sizes from HLA file
	hlas, sizes = parse_hla_file(parameters, mhc_class)

	if mhc_class == 'ii':
		nmers = list()
		for item in inputdata:
			nmers += split_item_nmers(item, int(parameters['mhciisizes']), int(parameters['nmerstep']))

		inputdata = nmers

	# Get peptides sequences from input data
	peptides = [item[2] for item in inputdata]

	# Add flanking characters needed by the API query
	peptides = ''.join(['%3Epeptide' + str(num) + '%0A' + pep.rstrip() + '%0A' for num, pep in enumerate(peptides, start = 1)])


	for attempt in range(1, 4): # 3 attempts to connect to the API

		success_connection_to_API = False

		command = "curl --data \"method=" + parameters['mhc'+mhc_class+'method'] + "&sequence_text="+peptides+"&allele=" + hlas + "&length="+ sizes +"\" http://tools-cluster-interface.iedb.org/tools_api/mhc"+mhc_class+"/"

		result = subprocess.run(command, shell=True, capture_output=True, text=True).stdout

		if result.strip('\t\n ') != '':
			success_connection_to_API = True
		
		if success_connection_to_API:			
			break
		else:
			print(f"Failed to connect to the IEDB API. Attempt {attempt}.")

		if attempt == 3:
			print('Failed to connect to the API after 3 attempts. Exiting.')
			exit(0)

	return result

# -------------------------------------
def parse_areas_file(parameters):
	''' Reads the input file (data/areas.txt) and 
		returns a list with the geographical regions 
	'''
	with open(parameters['areas']) as inputfile:
		return [line.rstrip() for line in inputfile]
	 
# -------------------------------------
def run_population_coverage(inputfile, parameters, mhc_class, areas=['World']):

	'''
		Calls the IEDB population coverage tool 
		from the dir specified in the parameters 
		file. 
		Args: mhci_prediction: return value of run_mhci_prediction()
	'''

	mhc_class = mhc_class.upper()

	if parameters['outputgraph'] != 'true':
		method_dir = join(parameters['populationcoveragedirectory'], 'no_graph')
	else:
		method_dir = join(parameters['populationcoveragedirectory'], 'graph')
		fig_dir = join(parameters['outputdirectory'], 'figures')
		makedirs(fig_dir, exist_ok=True)

	method_path = join(method_dir, 'calculate_population_coverage.py')
	
	py = parameters['pythonpath']

	# Runs the population tool for each sub-area
	coverage = dict()

	for area in areas:
		
		command = py + ' ' + method_path + ' -f ' + inputfile + ' -p ' + '"' + area + '"' +  ' -c ' + mhc_class + ' --plot ' + parameters['outputdirectory']

		result = subprocess.run(command, shell=True, capture_output=True, text=True)

		# Store the sequence and coverage in a dict
		try:
			coverage[area] = get_overall_results(result.stdout)

			if parameters['outputgraph'] == 'true':
				fmt_area = area.replace(' ', '_').lower()
				old_graphfile = join(parameters['outputdirectory'], 'popcov_' + fmt_area + '_' + mhc_class.lower() + '.png')
				fields = inputfile.split('.')   
				region = fields[1]
				locus  = fields[2]
				newfilename = fmt_area + '_MHC' + mhc_class + '_' + locus + '_' + region + '.png'
				new_graphfile = join(fig_dir, newfilename)
				rename(old_graphfile, new_graphfile)

		except:
			coverage[area] = '\t\t'

	return coverage

# -------------------------------------
def separate_hla_by_loci(hlas):

	dict_loci = dict()

	for hla in hlas.split(','):
		loci = hla.split('*')[0].split('-')[1]
		if loci in dict_loci:
			dict_loci[loci] += [hla]
		else:
			dict_loci[loci] = [hla]

	return {loci:','.join(hlas) for loci, hlas in dict_loci.items()}
# -------------------------------------
def map_peptides_to_regions(regions, peptides):
	dict_peptides_regions = dict()
	
	for peptide in peptides:
		for item in regions:
			if peptide in item[2]:
				name = item[0] + '-' + item[2]
				if name in dict_peptides_regions:
					dict_peptides_regions[name] += [peptide]
				else:
					dict_peptides_regions[name] = [peptide]

	return dict_peptides_regions
# -------------------------------------
def combine_cover_per_region(cover_per_region):
	d_hits = dict()
	d_pc90 = dict()
	d_coverage = dict()
	for epitope in cover_per_region:
		
		if epitope not in d_hits:
			d_hits[epitope] = dict()
			d_pc90[epitope] = dict()
			d_coverage[epitope] = dict()

		for peptide in cover_per_region[epitope]:
			for locus in cover_per_region[epitope][peptide]:
		
				if locus not in d_hits[epitope]:
					d_hits[epitope][locus] = dict()
					d_pc90[epitope][locus] = dict()
					d_coverage[epitope][locus] = dict()

				for area in cover_per_region[epitope][peptide][locus]:
					values = cover_per_region[epitope][peptide][locus][area].split('\t')

					coverage = values[0].strip('\t')
					hits = values[1].strip('\t')
					pc90 = values[2].strip('\t')

					if hits:
						hits = float(hits)
						pc90 = float(pc90)
						if area not in d_hits[epitope][locus]:
							d_hits[epitope][locus][area] = [hits]
							d_pc90[epitope][locus][area] = [pc90]
							d_coverage[epitope][locus][area] = coverage
						else:
							d_hits[epitope][locus][area] += [hits]
							d_pc90[epitope][locus][area] += [pc90]

	d_combined = dict()
	for epitope in d_hits:
		d_combined[epitope] = dict()
		for locus in d_hits[epitope]:
			d_combined[epitope][locus] = dict()
			for area in d_hits[epitope][locus]:
				hits = str(round(sum(d_hits[epitope][locus][area]),2))
				pc90 = str(round(sum(d_pc90[epitope][locus][area]),2))
				coverage = d_coverage[epitope][locus][area]
				d_combined[epitope][locus][area] = coverage + '\t' + hits +'\t' + pc90

	return d_combined
# -------------------------------------
def pop_coverage_single_region(epitope_regions, parameters, mhc_class, predictions):

	# Create the pop coverage input file
	pred_results = get_alleles_and_binders(predictions, parameters, mhc_class.lower())	
	if mhc_class.lower() == 'i':
		loci = ['A', 'B', 'any']
	if mhc_class.lower() == 'ii':
		loci = ['DP', 'DQ', 'DR', 'any']


	dict_regions_peptides = map_peptides_to_regions(epitope_regions, pred_results.keys())
	
	
	cover_per_region = dict()
	for region in  dict_regions_peptides:

		for locus in loci:
	
			for seq in dict_regions_peptides[region]:		
				hlas = pred_results[seq]
				num = region.split('-')[0]
				cover_per_locus = dict()
				
				if locus == 'any':
					# Creates the inputfile for the pop coverage tool for one loci
					inputfile = create_pop_coverage_input({seq:hlas}, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.' + num + '.'  + locus + '.input', locus=False)

				else:
					# Creates the inputfile for the pop coverage tool for all loci
					inputfile = create_pop_coverage_input({seq:hlas}, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.' + num + '.'  + locus + '.input', locus=locus)

			# Run the pop coverage input file 
			areas = parse_areas_file(parameters)
			coverage = run_population_coverage(inputfile, parameters, mhc_class.upper(), areas)

			# if locus not in cover_per_locus:
			cover_per_locus[locus] = coverage

			if region not in cover_per_region:
				cover_per_region[region] = cover_per_locus
			else:
				cover_per_region[region].update(cover_per_locus)

	return cover_per_region

# -------------------------------------
def get_overall_results(pop_coverage_output):

	'''
		Parse the output of the IEDB population coverage standalone tool
		str(pop_coverage_output) and returns a string containing the 
		percent value of the total coverage
	'''
	results = '\t'.join(pop_coverage_output.split('\n')[2].rstrip().split('\t')[1:])
	if not results:
		results = '\t\t'
	return results

# -------------------------------------
def output_to_table(mhci, mhcii, separator, filename):

	coverage_per_epitope = [k + separator + mhci[k] + separator + mhcii[k] for k in mhci]

	with open(filename, 'w') as out:

		# Header
		out.write('Region' + separator + 'MHC-I' + separator + 'MHC-II' + '\n')
		
		# Converts results to endline separated string
		out.write('\n'.join(coverage_per_epitope))

	return 

# -------------------------------------
def pop_coverage_all_regions(epitope_regions, parameters, mhc_class, prediction):


	# Isolate peptides and their respective binding alleles
	pred_results = get_alleles_and_binders(prediction, parameters, mhc_class)
	
	if mhc_class.lower() == 'i':
		loci = ['A', 'B']
	if mhc_class.lower() == 'ii':
		loci = ['DP', 'DQ', 'DR']

	# Get all major sub areas for the globe 
	areas = parse_areas_file(parameters)

	cover_per_locus = empty_dict(loci)

	for locus in loci:
		
		for seq in pred_results:
			hlas = pred_results[seq]

			# Creates the pop coverage input for each peptide individually
			inputfile = create_pop_coverage_input({seq:hlas}, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.allregions.'+locus+'.input', locus)
			
			# Run population coverage for the original sequences
			cover_per_locus[locus] = run_population_coverage(inputfile, parameters, mhc_class, areas)
	
	inputfile = False
	for seq in pred_results:
		hlas = pred_results[seq]
		inputfile = create_pop_coverage_input({seq:hlas}, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.allregions.input')

	if inputfile:
		cover_per_locus['any'] = run_population_coverage(inputfile, parameters, mhc_class, areas)
	else:
		# cover_per_locus['any'] = '\t\t' #EMPTY DICTIONARY
		cover_per_locus['any'] = empty_dict_regions()

	return cover_per_locus

# -------------------------------------
def save_prediction(prediction, parameters, mhc_class, name=''):
	with open(name, 'w') as pred_save:
		pred_save.write(prediction)

# -------------------------------------
def prediction_exists(file_path):

	# Check if file exists and is not empty
	if exists(file_path) and getsize(file_path) > 0:
		return True
	else:
		return False
# -------------------------------------
def parse_prediction_file(prediction_file_name):
	with open(prediction_file_name) as inputfile:
		return ''.join(inputfile.readlines())


# -------------------------------------
def run_API_prediction(parameters, mhc_class, epitope_items):

	# Check if the prediction files exist
	prediction_file_name = join(parameters['temporarydirectory'], 'MHC-' + mhc_class + '_' + epitope_items[0][2] + '.tsv')

	# Checks if user wants to reuse existing prediction
	answer = 'n'
	if prediction_exists(prediction_file_name):
		print('Prediction file', prediction_file_name, 'found. Do you want to use it (y/n)?', end=' ') 
		answer = input()

	# Reuses existing prediction
	if answer.lower() == 'y':
		prediction = parse_prediction_file(prediction_file_name)
	
	# Not reusing prediction: rerunning it.
	else:

		# Runs the MHC prediction tool specified in the parameter file
		prediction = run_mhc_prediction(epitope_items, parameters, mhc_class)
		# Save prediction to file
		save_prediction(prediction, parameters, mhc_class, name=prediction_file_name)
	
	return prediction
# -------------------------------------
def run(input_file, parameters, mhc_class):

	'''
		Coordinates the pipeline as parametrized
		by arg:parameters_file
	'''

	# Creates tmp dir and output dir
	makedirs(parameters['temporarydirectory'], exist_ok=True)
	
	# Gets the peptides from the csv input file
	epitope_items = parse_csv_input(input_file)

	# Use the API to run predictions
	prediction = ''
	while(prediction == ''):
		print(f'Using API to run predictions for MHC-{mhc_class}')
		prediction = run_API_prediction(parameters, mhc_class, epitope_items)

	# Run the pop coverage tool for each peptide separately
	coverage_dict = pop_coverage_single_region(epitope_items, parameters, mhc_class, prediction)

	# Run pop coverage tool for all regions combined and add to dictionary
	coverage_dict['all'] = pop_coverage_all_regions(epitope_items, parameters, mhc_class, prediction)
	
	return coverage_dict

# -------------------------------------
def merge_sequences(list_items):
	
	'''
		Auxiliary function to merge overlapping sequences
		Args: list_items = list of tuples that contain overlapping
		sequences. Format: (peptide number, start-end position, peptide sequence)
		Returns: a new list of tuples, with the overlapping sequences
		merged and start-end positions modified accordingly
	'''

	# List of tuples (position, aminoacid)
	positional_residues = list()

	for item in list_items:
		
		# Isolates start position, end position and AA sequence
		start = int(item[1].split('-')[0])
		end = int(item[1].split('-')[1])
		seq = item[2].rstrip()

		# Generates tuples (position, aminoacid) for all sequences
		positional_residues += list(enumerate(seq, start=start)) 

	# Removes redundancies; sorts list of tuples by position
	fullseq = sorted(list(set(positional_residues)))

	# Makes a new sequence (string) using aminoacids
	fullseq = ''.join([item[1] for item in fullseq])

	# # Gets start and end positions
	start = list_items[0][1].split('-')[0]
	end = list_items[-1][1].split('-')[1]

	return [list_items[0][0], start + '-' + end, fullseq]

# -------------------------------------
def group_items(list_items):

	'''
		Auxiliary function that receives a list of items
		and groups them into lists if the sequences overlap
	'''

	grouped_data = []
	current_group = []

	for sublist in list_items:

		# If the first element is not empty
		if sublist[0]:  

			# Appends to grouped data if current is not empty
			if current_group: 
				grouped_data.append(current_group) 

			current_group = [sublist]
		
		# If the first element is empty, append to current group
		else:
			current_group.append(sublist)
	
	# Append the last group if not empty
	if current_group:
		grouped_data.append(current_group)
	
	return grouped_data			

# -------------------------------------
def merge_items(data):

	'''
		Merges overlapping sequences
		Non-overlapping sequences are kept as is
	'''

	grouped_data = group_items(data)

	return [merge_sequences(item) for item in grouped_data]

# -------------------------------------
def split_item_nmers(item, nmer, step):

	split_list = list()

	num = item[0]
	pos = int(item[1].split('-')[0])
	seq = item[2]
	
	for i in range(0, len(seq)-nmer+1, step):
		split_list.append([num, pos, seq[i:i+nmer]])
		pos += step

	return split_list
# -------------------------------------
def parse_csv_input(csv_file):

	'''
		Parse the .csv input file. 
		Returns two lists of sequences:
		1. The merged overlapping sequences
		2. The combined regions split into nmers
	'''

	with open(csv_file) as csv:

		# Skips the first line (header)
		next(csv)

		separated = [line.rstrip().split(',') for line in csv.readlines()]
	
	# Merge overlapping peptides into one long sequence; outputs to fasta
	merged = merge_items(separated)

	return merged

# -------------------------------------
def empty_dict_regions():
	return {'Central Africa': '\t\t', 'Central America': '\t\t', 'East Africa': '\t\t', 'East Asia': '\t\t', 'Europe': '\t\t', 'North Africa': '\t\t', 'North America': '\t\t', 'Northeast Asia': '\t\t', 'Oceania': '\t\t', 'South Africa': '\t\t', 'South America': '\t\t', 'South Asia': '\t\t', 'Southeast Asia': '\t\t', 'Southwest Asia': '\t\t', 'West Africa': '\t\t', 'West Indies': '\t\t', 'World': '\t\t'}

# -------------------------------------
def empty_dict(loci):
	d = {}
	for locus in loci:
		d[locus] = empty_dict_regions()
	return d
# -------------------------------------
def output_to_files(coverage_mhci, coverage_mhcii, parameters):

	regions_mhci = set(coverage_mhci.keys())
	regions_mhcii = set(coverage_mhcii.keys())

	mhci_missing_regions = regions_mhcii - regions_mhci
	mhcii_missing_regions = regions_mhci - regions_mhcii

	for region in mhci_missing_regions:
		coverage_mhci[region] = empty_dict(['A', 'B', 'any'])
	
	for region in mhcii_missing_regions:
		coverage_mhcii[region] = empty_dict(['DP', 'DQ', 'DR', 'any'])


	for epitope in coverage_mhci:

		output_str = ''
		output_str += epitope + '\t' + 'Class-I' + 9*'\t' + 'Class-II' + '\n'

		loci_mhci_dict = coverage_mhci[epitope]
		loci_mhcii_dict = coverage_mhcii[epitope]
		
		loci_mhci_list = list(loci_mhci_dict.keys())
		loci_mhcii_list = [l for l in loci_mhcii_dict]
		
		formatted_mhci_list = ''.join([3*str(locus+'\t') for locus in loci_mhci_list])
		formatted_mhcii_list = ''.join([3*str(locus+'\t') for locus in loci_mhcii_list])	
		
		output_str += '\t' + formatted_mhci_list + formatted_mhcii_list  + '\n'
		output_str += 'Region\t' + 7*'Coverage\tAv.no hits\tPC90\t' + '\n'

		for region in loci_mhci_dict[loci_mhci_list[0]].keys():
			output_str += region  + '\t'
			for loci_mhci_list in loci_mhci_dict.keys():
				output_str += loci_mhci_dict[loci_mhci_list][region] + '\t'
			for loci_mhci_list in loci_mhcii_dict.keys():
				output_str += loci_mhcii_dict[loci_mhci_list][region] + '\t'
			output_str += '\n'

			with open(join(parameters['outputdirectory'], epitope+'.tsv'),'w') as outputfile:
				outputfile.write(output_str)
# -------------------------------------

if __name__ == '__main__':

	import argparse

	def parse_arguments():

		parser = argparse.ArgumentParser(description='Pipeline for population coverage of predicted MHC binders.')
		parser.add_argument('-i', required=True,  type=str, help='Input file')
		parser.add_argument('-p', required=False, type=str, help='Parameters file', default='parameters.md')
		args = parser.parse_args()

		return args

	args = parse_arguments()
	
	parameters = parse_parameters(args.p)
	outputdir = parameters['outputdirectory']
	makedirs(outputdir, exist_ok=True)
	
	# # Clears old tmp dir if it exists
	# if exists(parameters['temporarydirectory']):
	# 	rmtree(parameters['temporarydirectory'])

	coverage_mhci  = run(input_file = args.i, parameters = parameters, mhc_class = 'I')
	coverage_mhcii = run(input_file = args.i, parameters = parameters, mhc_class = 'II')
	output_to_files(coverage_mhci, coverage_mhcii, parameters)