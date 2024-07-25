
# Adjust this accordingly
output directory = results_smalltest_final
 
# Set parameters 
MHC-I method  = netmhcpan_ba
MHC-II method  = netmhciipan_ba
MHC-I sizes = 9,10
MHC-II sizes = 15
MHC-I threshold = 500
MHC-II threshold = 500
N-mer Step = 5
output graph = false

# More advanced options
Population coverage directory = population_coverage_standalone
HLA-I file = data/hla_ref_set.class_i.txt
HLA-II file = data/hla_ref_set.class_ii.txt
areas = data/areas.txt
temporary directory = data/temp_pickle
python path=python
MHC-I predictor directory = /home/trevizani/prog/mhc_i
MHC-II predictor directory = /home/trevizani/prog/mhc_ii