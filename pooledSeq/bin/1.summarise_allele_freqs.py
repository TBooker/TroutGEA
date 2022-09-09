# Combine allele freqs by population
import gzip
import sys


### MAF FILTER OF 5%

list_of_files_invar = ["/home/booker/work/Trout/pooledSeq/Atha/Atha-Buffalo_Prairie.alleleFreqs.mafs.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Cabin_Creek.alleleFreqs.mafs.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Deerlick_Creek.alleleFreqs.mafs.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Groat_Creek.alleleFreqs.mafs.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Marsh_Head_Creek.alleleFreqs.mafs.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Sakwatamau.alleleFreqs.mafs.gz"]

list_of_files = ["/home/booker/work/Trout/pooledSeq/Atha/Atha-Buffalo_Prairie.alleleFreqs.mafs.noInvar.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Cabin_Creek.alleleFreqs.mafs.noInvar.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Deerlick_Creek.alleleFreqs.mafs.noInvar.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Groat_Creek.alleleFreqs.mafs.noInvar.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Marsh_Head_Creek.alleleFreqs.mafs.noInvar.gz",
                "/home/booker/work/Trout/pooledSeq/Atha/Atha-Sakwatamau.alleleFreqs.mafs.noInvar.gz"]

pop_numbers = {"Atha-Buffalo_Prairie": 7,
                "Atha-Cabin_Creek": 33,
                "Atha-Deerlick_Creek":34,
                "Atha-Groat_Creek":35,
                "Atha-Marsh_Head_Creek":33,
                "Atha-Sakwatamau":28}

def overlap(x, y):
    return( max(x[0], y[0]) < min(x[1], y[1]))



def parseLine(raw, focal_chrom, switch = False):
    x = i.decode("UTF-8").strip().split()
    if x[0].startswith("chromo"): return None,None
    if switch:
        pass
    else:
        if not x[0].startswith(focal_chrom): return None,None
    chrom=x[0]
    pos = x[1]
    ref = x[2]
    alt = x[3]
    freq = float(x[4])
    inds = int(x[5])
    return chrom + ":" + pos, {"CHROM":chrom,
            "POS":int(pos),
            "REF":ref,
            "ALT":alt,
            "FREQ":freq,
            "INDS":inds}

chrom_ends = {"OmyA_13":74000000,
#            "OmyA_13":7400000,
            "OmyA_19":24000000,
            "OmyA_21":36000000,
            "OmyA_26":52000000}


tetrasomic_intervals = {"OmyA_13":[[2097076,15141548],[50600576,64186329]],
                        "OmyA_19":[[3345106,9316639]],
                        "OmyA_21":[[11851620,18630667]],
                        "OmyA_26":[[38659204,51104312]]}

print("Tetrasomic intervals:\n\t",tetrasomic_intervals)


all_positions = {}

for focal_chrom in ["OmyA_13", "OmyA_19", "OmyA_21","OmyA_26"]:
#for focal_chrom in ["OmyA_13"]:
    tetrasome_intervals =tetrasomic_intervals[focal_chrom]

    intervals_raw = [[j,j+2000000] for j in range(0,
                                                chrom_ends[focal_chrom],
                                                2000000)]

    intervals = []

    for ti in tetrasome_intervals:
        for inty in intervals_raw:
            if  overlap(ti, inty):
                intervals.append( inty )
    print("\nAnalysing Chromosome:",focal_chrom)

    print("\tGetting list of SNP positions\n")
    for current_file in list_of_files:
        pop = current_file.split("/")[-1].split(".")[0]
        print("\t\tGrabbing variant positions from:",pop)
        pop_total = pop_numbers[pop]
        count = 0
        pop_dict = {}
        for i in gzip.open(current_file, "rb"):
#            if len(all_positions.keys()) == 20000:break
            label, line_dict = parseLine(i, focal_chrom)
            if label == None:continue
            if line_dict["CHROM"]!=focal_chrom:continue
            if focal_chrom!="OmyA_13":
                if line_dict["POS"]>tetrasome_intervals[0][1]:break
                elif (line_dict["POS"] > tetrasome_intervals[0][0] and line_dict["POS"] < tetrasome_intervals[0][1] ):pass
                else:continue
            elif focal_chrom =="OmyA_13":
                if (line_dict["POS"] > tetrasome_intervals[0][0] and line_dict["POS"] < tetrasome_intervals[0][1] ) :
                    pass
                elif (line_dict["POS"] > tetrasome_intervals[1][0] and line_dict["POS"] < tetrasome_intervals[1][1] ) :
                    pass
                else:continue

#            elif line_dict["POS"] > inty[1]  :
#                break
#            else:continue
#            print(inty, line_dict)

            if line_dict["FREQ"]<0.01: continue
            if  line_dict["INDS"]/pop_total <0.5:continue
            all_positions[ label ] = 0
#            print(tetrasome_intervals,label)
#            if count == 30000:break
#            print(line_dict)
            count +=1


print("\n\nExtracting data for putative SNPs\n")
#    snp_DFs = []
#    first = True

#    for interval in intervals:
#        print(interval)
#        all = {}

generic_header = ["CHROM","POS","REF","ALT","PROP_INDS","FREQ"]


for current_file in list_of_files_invar:
    pop = current_file.split("/")[-1].split(".")[0]
    print("\t",pop)
    pop_total = pop_numbers[pop]
    header = generic_header[:2] + [pop+"."+g for g in generic_header[2:]]

    this_pop_this_chrom_output = open(pop+"_"+"alleleFreqs.csv", "w")
    this_pop_this_chrom_output.write(",".join(header) + "\n")
    for i in gzip.open(current_file, "rb"):
        label, line_dict = parseLine(i, focal_chrom, switch = True)
        if label == None:continue
#        if line_dict["CHROM"]!=focal_chrom:continue
        try:
            all_positions[ label ]
        except KeyError:
            continue
        print(line_dict)

        output_line = map(str,  [line_dict["CHROM"], line_dict["POS"], line_dict["REF"], line_dict["ALT"],line_dict["INDS"]/pop_total,  line_dict["FREQ"]])
        this_pop_this_chrom_output.write(",".join(output_line) + "\n")

    this_pop_this_chrom_output.close()
