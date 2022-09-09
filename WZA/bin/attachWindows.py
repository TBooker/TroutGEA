# Attach windows to correlation files:

import numpy as np
import pandas as pd
import sys
import scipy.stats

# Read in the allele freq data
data  = pd.read_csv( sys.argv[1] )

print("naming the output",sys.argv[2])
max_pos = data.POS.max()

windows = "bin/window_positions_entire_genome.csv"

focal_chrom = "OmyA_19"

win_df = pd.read_csv( windows )



win_df = win_df[(win_df.Chromosome == focal_chrom)&(win_df["Window.Start"]<max_pos)].reset_index()

v = win_df.loc[:,"Window.Start":"Window.End"].apply(tuple,1).to_list()

idx = pd.IntervalIndex.from_tuples(v,closed = "neither")

print(idx.is_overlapping)

#print( idx.get_indexer(data.POS.values) )

#data["window"] = win_df.iloc[idx.get_indexer_non_unique(data.POS.values),"Window"].values


output = []
for i, r in data.iterrows():
#    print(i, r)
    interval_df = idx.get_indexer_non_unique([r.POS])
    if len( interval_df[0]) != 1:
        print("!")
        output.append(-99)
        continue
    output.append( win_df.Window.values[interval_df[0][0]] )


#    print( win_df.loc( idx.get_indexer_non_unique([r.POS])[0][0] ) )
#    print( win_df[(win_df.Chromosome == r.CHROM)&(win_df["Window.Start"] <= r.POS )&(win_df["Window.End"] >= r.POS )] )
data["Window"] = output

data.to_csv(sys.argv[2], index = False)
