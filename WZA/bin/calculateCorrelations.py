##  WZA prep for tetrasomes

## We need to calculate the correlation between allele frequency and the environment for each locus, storing the p-value from each

import numpy as np
import pandas as pd
import sys
import scipy.stats

# Read in the allele freq data
data  = pd.read_csv( sys.argv[1] )

# Read the env data
envs = pd.read_csv( sys.argv[2], sep = "\t" )

# Use Jared's 'Site' column as the DF index
envs = envs.set_index( "Site")
# Build dicts for each population
env_dict = envs.to_dict()

# Build a dict relating names of each pop in Jared's DF to the names I've used:
pop_name_dict = {"Atha-Buffalo_Prairie":"Buffalo Prairie",
                "Atha-Cabin_Creek":"Cabin Creek",
                "Atha-Deerlick_Creek":"Deerlick",
                "Atha-Groat_Creek":"Groat Creek",
                "Atha-Marsh_Head_Creek":"Marsh Head Creek",
                "Atha-Sakwatamau":"Sakwatamau River"}

# Get a list of the names of the environmental variables
list_of_envs = [ e for e in env_dict["Deerlick"].keys() if e!="Watershed" ]


metadat_indices = [o for o in range(len(list(data))) if list(data)[o] in ["CHROM","POS","p_bar","BASES"]]

population_indices = [o for o in range(len(list(data))) if list(data)[o] not in ["CHROM","POS","p_bar","BASES"]]

population_names = [list(data)[o] for o in range(len(list(data))) if list(data)[o] not in ["CHROM","POS","p_bar","BASES"]]

translated_population_names = [pop_name_dict[p] for p in population_names]

population_names = []

env_vector_dict = {}

for e in list_of_envs:
    env_vector_dict[e] =  np.array( [ float(env_dict[p][e]) for p in translated_population_names ] )
    print(env_vector_dict[e])

print("###############################################")
print()

count = 0

# Initialise the output file
output_file = open(sys.argv[3], "w")

header_part_1 =  [k for k in ["CHROM","POS","p_bar"]]
header_part_2 = []

for evt in list_of_envs:
    header_part_2.append(evt+"_cor")
    header_part_2.append(evt+"_pval")

header = header_part_1 + header_part_2

output_file.write( ",".join(header) + "\n")

for i, r in data.iterrows():

    allele_freq_vector = np.array([r[k] for k in population_indices])
    output = [str(r[m]) for  m in metadat_indices]
    output.append(allele_freq_vector.mean())
    for ev in list_of_envs:
        env_vector = env_vector_dict[ev]

        cor = scipy.stats.kendalltau( allele_freq_vector, env_vector )
        output.append( cor.correlation )
        output.append( cor.pvalue )
    output_file.write( ",".join( list(map(str, output) ) )+"\n" )

output_file.close()
