#!/usr/bin/python

import subprocess
from os.path import join
from os import makedirs, rename

# -------------------------------------	
def parse_parameters(inputfile):

	'''
		Reads input file containing parameters (default: parameters.txt)
		and returns a dictionary {variables:values} as specified in
		the input file as variables = values
	'''

	parameters = dict()

	with open(inputfile) as par:
		for line in par:

			# Skip comment lines
			if line.startswith('#'):
				continue
   
			# Removes non-alphanumeric extra characters; converts to lowercase
			line = line.strip('\t\n\'\"').replace(' ', '').replace('-', '')

			# If line isn't empty after cleaning
			if line: 

				# File format: variable = value
				k,v = line.split('=') 

				# Removes non-alphanumeric extra characters
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
def create_pop_coverage_input(pred_results, parameters, filename): 

	'''
		Creates the input file needed to call 
		the population coverage tool
		pred_results: dictionary {str(epitope):str('hla1,hla2')}
		Returns: created file name
	'''
	
	# File with population coverage input in tmpdir/:
	pop_input_name = join(parameters['temporarydirectory'], filename)

	# Writes to file
	with open(pop_input_name, 'w') as f:
		for peptide, allele_list in pred_results.items():
			line = peptide + '\t' + allele_list + '\n'
			f.write(line)

	# Returns file name
	return pop_input_name

# -------------------------------------
def run_mhc_prediction(inputdata, parameters, mhc_class):
	'''
		Uses the API to predict binders
	'''

	# Ensures the mhc class is not capitalized
	mhc_class = mhc_class.lower()
	
	# Gets HLAs and respective sizes from HLA file
	hlas, sizes = parse_hla_file(parameters, mhc_class)
	
	# Get peptides sequences from input data
	peptides = [item[2] for item in inputdata]

	# Add flanking characters needed by the API query
	peptides = ''.join(['%3Epeptide' + str(num) + '%0A' + pep.rstrip() + '%0A' for num, pep in enumerate(peptides, start = 1)])

	command = "curl --data \"method=" + parameters['mhc'+mhc_class+'method'] + "&sequence_text="+peptides+"&allele=" + hlas + "&length="+ sizes +"\" http://tools-cluster-interface.iedb.org/tools_api/mhc"+mhc_class+"/"
	
	result = subprocess.run(command, shell=True, capture_output=True, text=True)

	return result.stdout

# -------------------------------------
def run_population_coverage(inputfile, parameters, mhc_class):

	'''
		Calls the IEDB population coverage tool 
		from the dir specified in the parameters 
		file. 
		Args: mhci_prediction: return value of run_mhci_prediction()
	'''

	mhc_class = mhc_class.upper()

	method_path = join(parameters['populationcoveragedirectory'], 'calculate_population_coverage.py')

	py = parameters['pythonpath']
	command = py + ' ' + method_path + ' -f ' + inputfile + ' -p ' + parameters['coveragearea'] + ' -c ' + mhc_class + ' --plot ' + parameters['outputdirectory']

	result = subprocess.run(command, shell=True, capture_output=True, text=True)
		
	return result.stdout

# -------------------------------------
def combine_peptides(joined_data, dict_pep_hla, parameters):

	# Empty dictionary where original sequences are keys that maps to a set of HLAs
	sequences = [item[2] for item in joined_data]
	combined_pep_hla = {seq:list() for seq in sequences}
	
	for key_ind in dict_pep_hla:
		for key_comb in combined_pep_hla:
			if key_ind in key_comb:
				combined_pep_hla[key_comb] += dict_pep_hla[key_ind].split(',')


	combined_pep_hla = {k:','.join(set(v)) for k,v in combined_pep_hla.items()}
	
	return combined_pep_hla

# -------------------------------------
def pop_coverage_mhc(input_file, parameters, mhc_class):

	# Runs the MHC-I prediction tool specified in the parameter file
	prediction = run_mhc_prediction(input_file, parameters, mhc_class)
	
	# Isolate peptides and their respective binding alleles
	pred_results = get_alleles_and_binders(prediction, parameters, mhc_class)
		
	# Combines the overlapping peptides and merges HLA sets
	pred_results_combined = combine_peptides(pred_results, parameters)

	# Creates the pop coverage input for each peptide individually
	inputfile = create_pop_coverage_input(pred_results_combined, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.original.input')
	
	# Run population coverage for the original sequences
	coverage_all_seqs = run_population_coverage(inputfile, parameters, mhc_class)

	return coverage_all_seqs

# -------------------------------------
def pop_coverage_single_region(joined_data, parameters, mhc_class, predictions):

	sequences = [item[2] for item in joined_data]

	# Create the pop coverage input file
	pred_results = get_alleles_and_binders(predictions, parameters, mhc_class.lower())	

	pred_results = combine_peptides(joined_data, pred_results, parameters)

	individual_cover = dict()
	
	for num, seq in enumerate(sequences, start=1):

		# Selects the prediction of one sequence
		pred = pred_results[seq]

		# Creates the inputfile for the pop coverage tool
		inputfile = create_pop_coverage_input({seq:pred}, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.' + str(num) + '.input')
		
		# Run the pop coverage input file 
		pop_seq = run_population_coverage(inputfile, parameters, mhc_class.upper())

		# Store the sequence and coverage in a dict
		coverage = get_overall_coverage(pop_seq)
		individual_cover[str(num) + '-' + seq] = coverage

	return individual_cover

# -------------------------------------
def get_overall_coverage(pop_coverage_output):

	'''
		Parse the output of the IEDB population coverage standalone tool
		str(pop_coverage_output) and returns a string containing the 
		percent value of the total coverage
	'''

	return pop_coverage_output.split('\n')[2].split('\t')[1]

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
def run(input_file, parameters, mhc_class):

	'''
		Coordinates the pipeline as parametrized
		by arg:parameters_file
	'''

	# Creates tmp dir and output dir
	makedirs(parameters['temporarydirectory'], exist_ok=True)
	
	separated_peptides, joined_peptides = parse_csv_input(input_file, parameters)

	# Runs the MHC-I prediction tool specified in the parameter file
	prediction = run_mhc_prediction(separated_peptides, parameters, mhc_class)
	
	# Isolate peptides and their respective binding alleles
	pred_results = get_alleles_and_binders(prediction, parameters, mhc_class)

	# Merges HLA sets according to overlapping peptides
	pred_results_combined = combine_peptides(joined_peptides, pred_results, parameters)

	# Creates the pop coverage input for each peptide individually
	inputfile = create_pop_coverage_input(pred_results_combined, parameters, 'pop_coverage_mhc' + mhc_class.lower() + '.original.input')

	# Run population coverage for the original sequences
	full_coverage_output = run_population_coverage(inputfile, parameters, mhc_class)
	full_coverage = get_overall_coverage(full_coverage_output)

	# Run the pop coverage tool for each peptide separately
	coverage_dict = pop_coverage_single_region(joined_peptides, parameters, mhc_class, prediction)
	coverage_dict['all'] = full_coverage
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
def parse_csv_input(csv_file, parameters):

	'''
		Parse the .csv input file. 
		Returns two lists of sequences:
		1. the original 15-mers and 
		2. the merged overlapping sequences
	'''

	with open(csv_file) as csv:

		# Skips the first 3 lines
		for _ in range(3):
				next(csv)

		separated = [line.rstrip().split(',') for line in csv.readlines()]
	
	# Merge overlapping peptides into one long sequence; outputs to fasta
	merged = merge_items(separated)
 
	return separated, merged


# -------------------------------------
if __name__ == '__main__':

	import argparse

	def parse_arguments():

		parser = argparse.ArgumentParser(description='Pipeline for population coverage of predicted MHC binders.')
		parser.add_argument('-i', required=True,  type=str,   help='Input file')
		parser.add_argument('-p', required=False,  type=str,   help='Parameters file', default='parameters.md')
		args = parser.parse_args()

		return args

	args = parse_arguments()
	
	parameters = parse_parameters(args.p)
	outputdir = parameters['outputdirectory']
	makedirs(outputdir, exist_ok=True)
	
	coverage_mhci  = run(input_file = args.i, parameters = parameters, mhc_class = 'I')
	coverage_mhcii = run(input_file = args.i, parameters = parameters, mhc_class = 'II')
	
	output_to_table(coverage_mhci, coverage_mhcii, separator='\t', filename=join(parameters['outputdirectory'],'final_table.tsv'))