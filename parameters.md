
# Adjust these accordingly
output directory = test/results

# Set parameters 
MHC-I method  = netmhcpan_ba
MHC-II method  = netmhciipan_ba
MHC-I sizes = 9,10
MHC-II sizes = 15
Coverage area = World
MHC-I threshold = 500
MHC-II threshold = 500

# More advanced options
Population coverage directory = population_coverage_standalone
HLA-I file = data/hla_ref_set.class_i.txt
HLA-II file = data/hla_ref_set.class_ii.txt
temporary directory = data/temp
python path=python
