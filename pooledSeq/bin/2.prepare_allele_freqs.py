import pandas as pd

list_of_files = ["TEMP_Buff", "TEMP_Groat", "TEMP_Cabin", "TEMP_Marsh", "TEMP_Deer", "TEMP_Sak"]

list_of_files =["Atha-Buffalo_Prairie_alleleFreqs.csv",  "Atha-Marsh_Head_Creek_alleleFreqs.csv","Atha-Cabin_Creek_alleleFreqs.csv", "Atha-Sakwatamau_alleleFreqs.csv","Atha-Deerlick_Creek_alleleFreqs.csv","Atha-Groat_Creek_alleleFreqs.csv"]


counter = 0
for l in list_of_files:
    counter += 1
    if counter == 1:
        big_df = pd.read_csv(l)
    else:
        tmp_df = pd.read_csv(l)
        big_df = big_df.merge(tmp_df, on = ["CHROM", "POS"], how = "outer")

    for i in list( big_df ):
        print(i)
    print()
tmp_df = 0
big_df.to_csv("Atha_alleleFreqs.csv", index = False)
