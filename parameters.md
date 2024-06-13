
# Adjust this accordingly
output directory = test/results

# Set parameters 
MHC-I method  = netmhcpan_ba
MHC-II method  = netmhciipan_ba
MHC-I sizes = 9,10
MHC-II sizes = 15
MHC-I threshold = 500
MHC-II threshold = 500
Coverage area = World, East Asia, Northeast Asia, South Asia, Southeast Asia, Southwest Asia, Europe, East Africa, West Africa, Central Africa, North Africa, South Africa, West Indies, North America, Central America, South America, Oceania

# More advanced options
Population coverage directory = population_coverage_standalone
HLA-I file = data/hla_ref_set.class_i.txt
HLA-II file = data/hla_ref_set.class_ii.txt
temporary directory = data/temp
python path=python
