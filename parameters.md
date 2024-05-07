
# Adjust these accordingly
MHC-I predictor directory	= /home/trevizani/prog/mhc_i
MHC-II predictor directory	= /home/trevizani/prog/mhc_ii
Population coverage directory = /home/trevizani/prog/population_coverage
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
HLA-I file = data/hla_ref_set.class_i.txt
HLA-II file = data/hla_ref_set.class_ii.txt
temporary directory = data/temp
python path=python
remove temporary directory = false