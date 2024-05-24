#!/usr/bin/python

import subprocess
from os.path import join
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
def get_alleles_and_binders(prediction, parameters, mhci_or_mhcii):

	'''
		Isolates from IEDB prediction output 
		the tuples of peptides and corresponding alleles
		that are above the binding threshold
	'''

	# Decides whether it is MHC-I or MHC-II
	if mhci_or_mhcii == 'mhci':
		peptide_column = 5
		cutoff = float(parameters['mhcithreshold'])
	if mhci_or_mhcii == 'mhcii':
		peptide_column = 6
		cutoff = float(parameters['mhciithreshold'])
	

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
		if mhci_or_mhcii == 'mhcii' and '/' in allele:
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
def run_mhc_prediction(inputfile, parameters, mhc_class):
	'''
		Uses the API
	'''

	# Ensures the mhc class is not capitalized
	mhc_class = mhc_class.lower()
	
	# Gets HLAs and respective sizes from HLA file
	hlas, sizes = parse_hla_file(parameters, mhc_class)

	peptides = list()
	with open(inputfile) as inp:
		for line in inp:
			if not line.startswith('>'):
				peptides.append(line.rstrip())
	peptides = list(filter(None, peptides))
	peptides = ''.join(['%3Epeptide' + str(num) + '%0A' + pep.rstrip() + '%0A' for num, pep in enumerate(peptides, start = 1)])

	command = "curl --data \"method=" + parameters['mhc'+mhc_class+'method'] + "&sequence_text="+peptides+"&allele=" + hlas + "&length="+ sizes +"\" http://tools-cluster-interface.iedb.org/tools_api/mhc"+mhc_class+"/"
	
	result = subprocess.run(command, shell=True, capture_output=True, text=True)

	return result.stdout

# -------------------------------------
def run_population_coverage(inputfile, parameters, mhci_or_mhcii):

	'''
		Calls the IEDB population coverage tool 
		from the dir specified in the parameters 
		file. 
		Args: mhci_prediction: return value of run_mhci_prediction()
	'''

	method_path = join(parameters['populationcoveragedirectory'], 'calculate_population_coverage.py')

	py = parameters['pythonpath']
	command = py + ' ' + method_path + ' -f ' + inputfile + ' -p ' + parameters['coveragearea'] + ' -c ' + mhci_or_mhcii + ' --plot ' + parameters['outputdirectory']

	result = subprocess.run(command, shell=True, capture_output=True, text=True)
		
	return result.stdout

# -------------------------------------
def combine_peptides(dict_pep_hla, parameters):

	# Parse sequences from merged fasta file:
	merged_fasta_file = join(parameters['temporarydirectory'], 'merged.fasta.txt')
	with open(merged_fasta_file) as fasta:
		sequences = [line.rstrip() for line in fasta if not line.startswith('>')]
	
	# Empty dictionary where original sequences are keys that maps to a set of HLAs
	combined_pep_hla = {seq:list() for seq in sequences}
	
	for key_ind in dict_pep_hla:
		for key_comb in combined_pep_hla:
			if key_ind in key_comb:
				combined_pep_hla[key_comb] += dict_pep_hla[key_ind].split(',')


	combined_pep_hla = {k:','.join(set(v)) for k,v in combined_pep_hla.items()}
	
	return combined_pep_hla

# -------------------------------------
def save_prediction_to_file(prediction, parameters, filename):
	with open(join(parameters['temporarydirectory'], filename), 'w') as out:
		out.write(prediction)
	return
# -------------------------------------
def pop_coverage_mhci(input_file, parameters):

	# Runs the MHC-I prediction tool specified in the parameter file
	mhci_prediction = run_mhc_prediction(input_file, parameters, 'I')
	
	# Saves MHC-I prediction to a file for future use/debuggin
	save_prediction_to_file(mhci_prediction, parameters, 'mhc_i_prediction')

	# Isolate peptides and their respective binding alleles
	pred_results = get_alleles_and_binders(mhci_prediction, parameters, mhci_or_mhcii='mhci')
	
	# Creates the pop coverage input for each peptide individually
	inputfile = create_pop_coverage_input(pred_results, parameters, 'pop_coverage_mhci.split.input')
	
	# Run the population coverage tool for the peptides split in the sizes of parameters.md
	pop_mhci_separated = run_population_coverage(inputfile, parameters, 'I')

	# Changes the default figure name to reflect the 'split' nature of the peptides
	png_old = join(parameters['outputdirectory'], 'popcov_world_i.png')
	png_new = join(parameters['outputdirectory'], 'popcov_world_i_split.png')
	rename(png_old, png_new)
	
	# Run the population coverage tool for the peptides in their original sizes
	pred_results_combined = combine_peptides(pred_results, parameters)

	# Creates the pop coverage input for each peptide individually
	inputfile = create_pop_coverage_input(pred_results_combined, parameters, 'pop_coverage_mhci.original.input')
	
	# Run population coverage for the original sequences
	pop_mhci_full = run_population_coverage(inputfile, parameters, 'I')

	# Changes the default figure name to reflect the 'original' nature of the peptides
	png_old = join(parameters['outputdirectory'], 'popcov_world_i.png')
	png_new = join(parameters['outputdirectory'], 'popcov_world_i_original.png')
	rename(png_old, png_new)

	return pop_mhci_separated, pop_mhci_full

# -------------------------------------
def pop_coverage_mhcii(input_file, parameters):
	
	# Runs the MHC-II prediction tool specified in the parameter file
	mhcii_prediction = run_mhc_prediction(input_file, parameters, 'II')

	# Saves MHC-I prediction to a file for future use/debuggin
	save_prediction_to_file(mhcii_prediction, parameters, 'mhc_ii_prediction')

	# Isolate peptides and their respective binding alleles
	pred_results = get_alleles_and_binders(mhcii_prediction, parameters, mhci_or_mhcii='mhcii')

	# Outputs to file
	inputfile = create_pop_coverage_input(pred_results, parameters, 'pop_coverage_mhcii.split.input')

	# Run the population coverage tool in the inputfile
	pop_mhcii = run_population_coverage(inputfile, parameters, 'II')

	# Changes the default figure name to reflect the 'split' nature of the peptides
	png_old = join(parameters['outputdirectory'], 'popcov_world_ii.png')
	png_new = join(parameters['outputdirectory'], 'popcov_world_ii_split.png')
	rename(png_old, png_new)
	
	# Run the population coverage tool for the peptides in their original sizes
	pred_results_combined = combine_peptides(pred_results, parameters)

	# Creates the pop coverage input for each peptide individually
	inputfile = create_pop_coverage_input(pred_results_combined, parameters, 'pop_coverage_mhcii.original.input')
	
	# Run population coverage for the original sequences
	pop_mhcii_full = run_population_coverage(inputfile, parameters, 'II')

	# Changes the default figure name to reflect the 'original' nature of the peptides
	png_old = join(parameters['outputdirectory'], 'popcov_world_ii.png')
	png_new = join(parameters['outputdirectory'], 'popcov_world_ii_original.png')
	rename(png_old, png_new)

	return pop_mhcii, pop_mhcii_full

# -------------------------------------
def pop_coverage_mhcii_single_region(parameters):

	tmpdir = parameters['temporarydirectory']

	# Parse the merged fasta file TODO: isolate this into a sequence (used in more places)
	merged_fasta_file = join(tmpdir, 'merged.fasta.txt')
	with open(merged_fasta_file) as fasta:
		sequences = [line.rstrip() for line in fasta if not line.startswith('>')]

	# Parse the mhci input file
	with open(join(tmpdir, 'mhc_ii_prediction')) as predfile:
		predictions = '\n'.join(predfile.readlines())

	# Create the pop coverage input file
	pred_results = get_alleles_and_binders(predictions, parameters, mhci_or_mhcii='mhcii')	
	pred_results = combine_peptides(pred_results, parameters)

	individual_cover = dict()
	for num, seq in enumerate(sequences, start=1):

		# Selects the prediction of one sequence
		pred = pred_results[seq]

		# Creates the inputfile for the pop coverage tool
		inputfile = create_pop_coverage_input({seq:pred}, parameters, 'pop_coverage_mhcii.' + str(num) + '.input')
		
		# Run the pop coverage input file 
		pop_mhcii_seq = run_population_coverage(inputfile, parameters, 'II')

		# Store the sequence and coverage in a dict
		coverage = get_overall_coverage(pop_mhcii_seq)
		individual_cover[str(num) + '-' + seq] = coverage

	return individual_cover


# -------------------------------------
def pop_coverage_mhci_single_region(parameters):

	tmpdir = parameters['temporarydirectory']

	# Parse the merged fasta file TODO: isolate this into a sequence (used in more places)
	merged_fasta_file = join(tmpdir, 'merged.fasta.txt')
	with open(merged_fasta_file) as fasta:
		sequences = [line.rstrip() for line in fasta if not line.startswith('>')]

	# Parse the mhci input file
	with open(join(tmpdir, 'mhc_i_prediction')) as predfile:
		predictions = '\n'.join(predfile.readlines())

	# Create the pop coverage input file
	pred_results = get_alleles_and_binders(predictions, parameters, mhci_or_mhcii='mhci')	
	pred_results = combine_peptides(pred_results, parameters)

	individual_cover = dict()
	for num, seq in enumerate(sequences, start=1):

		# Selects the prediction of one sequence
		pred = pred_results[seq]

		# Creates the inputfile for the pop coverage tool
		inputfile = create_pop_coverage_input({seq:pred}, parameters, 'pop_coverage_mhci.' + str(num) + '.input')
		
		# Run the pop coverage input file 
		pop_mhci_seq = run_population_coverage(inputfile, parameters, 'I')

		# Store the sequence and coverage in a dict
		coverage = get_overall_coverage(pop_mhci_seq)
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
def run(input_file, parameters_file):

	'''
		Coordinates the pipeline as parametrized
		by arg:parameters_file
	'''

	# Parses parameter file, adds the input file as a parameter
	parameters = parse_parameters(parameters_file)
	parameters['inputfile'] = input_file

	# Creates tmp dir and output dir
	makedirs(parameters['temporarydirectory'], exist_ok=True)
	outputdir = parameters['outputdirectory']
	makedirs(outputdir, exist_ok=True)
	
	# TODO Attempt to detect file type
	separated_peptides, joined_peptides = convert_csv_to_fasta(input_file, parameters)
	
	# Runs pipeline for MHC-I; saves to file
	pop_mhci_separated, pop_mhci_full = pop_coverage_mhci(separated_peptides, parameters)
	with open(join(outputdir,'mhc_i_pop_coverage_split.txt'), 'w') as out:
		out.write(pop_mhci_separated)
	with open(join(outputdir,'mhc_i_pop_coverage_original.txt'), 'w') as out:
		out.write(pop_mhci_full)

	# Runs pipeline for MHC-II; saves to file
	pop_mhcii_separated, pop_mhcii_full = pop_coverage_mhcii(separated_peptides, parameters)
	with open(join(outputdir,'mhc_ii_pop_coverage_split.txt'), 'w') as out:
		out.write(pop_mhcii_separated)
	with open(join(outputdir,'mhc_ii_pop_coverage_original.txt'), 'w') as out:
		out.write(pop_mhcii_full)


	# Run the pop coverage tool for each MHC-I prediction separately
	pop_mhci_single_region = pop_coverage_mhci_single_region(parameters)

	# Run the pop coverage tool for each MHC-II prediction separately
	pop_mhcii_single_region = pop_coverage_mhcii_single_region(parameters)

	# Parse prediction for overall_coverage value and adds to dictionary
	pop_mhci_single_region['all'] = get_overall_coverage(pop_mhci_full)
	pop_mhcii_single_region['all'] = get_overall_coverage(pop_mhcii_full)

	
	# Create output table
	## DEPRECATED output_to_table(pop_mhci_single_region, parameters,  separator='\t', filename=join(outputdir,'table_mhci.tsv'))
	## DEPRECATED output_to_table(pop_mhcii_single_region, parameters, separator='\t', filename= join(outputdir,'table_mhcii.tsv'))
	output_to_table(pop_mhci_single_region, pop_mhcii_single_region, separator='\t', filename=join(outputdir,'final_table.tsv'))

	# Deletes temporary data dir 
	rm = parameters['removetemporarydirectory'].lower()
	if rm == 'true' or rm == 'y' or rm == 'yes':
		rmtree(parameters['temporarydirectory'])

	return
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
		and groups them if the sequences overlap
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
		Groups and merges overlapping sequences
		Non-overlapping sequences are kept as is
	'''

	grouped_data = group_items(data)
	return [merge_sequences(item) for item in grouped_data]

# -------------------------------------
def output_to_fasta(datalist, filename):

	'''
		Creates a fasta file from the list of tuples
		datalist format: (peptide number, start-end position, peptide sequence)
	'''

	with open(filename, 'w') as out:
		for item in datalist:
			number = item[0]
			positions = item[1]
			sequence = item[2]
			out.write('> ' + number + ' pos: ' + positions + '\n')
			out.write(sequence + '\n')
	return
# -------------------------------------
def convert_csv_to_fasta(csv_file, parameters):

	'''
		Automates converting csv files from the Sette lab
		to fasta files so the MHC prediction tools can be used
	'''

	with open(csv_file) as csv:

		# Skips the first 3 lines
		for _ in range(3):
				next(csv)

		lines = [line.split(',') for line in csv.readlines()]

	# Outputs to fasta each overlapping peptide
	fasta_individual = join(parameters['temporarydirectory'], 'individual.fasta.txt')
	output_to_fasta(lines, fasta_individual)

	# Merge overlapping peptides into one long sequence; outputs to fasta
	merged_lines = merge_items(lines)
	fasta_merged = join(parameters['temporarydirectory'], 'merged.fasta.txt')
	output_to_fasta(merged_lines, fasta_merged)
 
	return fasta_individual, fasta_merged

# -------------------------------------
def parse_arguments():

	import argparse

	parser = argparse.ArgumentParser(description='Pipeline for population coverage of predicted MHC binders.')
	parser.add_argument('-i', required=True,  type=str,   help='Input file')
	parser.add_argument('-p', required=False,  type=str,   help='Parameters file', default='parameters.md')
	args = parser.parse_args()

	return args

# -------------------------------------
if __name__ == '__main__':
	args = parse_arguments()
	run(input_file = args.i, parameters_file = args.p)