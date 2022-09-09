import sys


def most_common(lst):
    return max(set(lst), key = lst.count)

output = open(sys.argv[1], "w")

counter = 0
for line in open("Atha_alleleFreqs.csv"):
    x = line.strip().split(",")

    if line.startswith("CHROM"):
        ref_indices = [i for i in range(len(x)) if x[i].endswith("REF")]
        alt_indices = [i for i in range(len(x)) if x[i].endswith("ALT")]
        freq_indices = [i for i in range(len(x)) if x[i].endswith("FREQ")]
        ind_indices = [i for i in range(len(x)) if x[i].endswith("INDS")]
        pop_names = [x[i].split(".")[0] for i in range(len(x)) if x[i].endswith("INDS")]
        output.write( "CHROM,POS,"+",".join(pop_names)+'\n')
        continue

    tmp = [x[i] for i in freq_indices if x[i] =='']
    # If any population has missing data, move to the next SNP
    if len(tmp)!=0:
        continue
    freqs = [float(x[i]) for i in freq_indices]
#    print(line)
#    print(x)
#    print(freqs)

    allele_freqs = []

    bases =  ["".join(sorted(x[i]+x[j])) for i,j in zip(ref_indices,alt_indices)]
    canonical_bases =  most_common( bases )
    ref_base = canonical_bases[0]

    snp_data = {}

# We have a lot of SNPs in this weird pooled-Seq analysis
# So let's not feel too bad about removing SNPs
    prop_inds = [float(x[k]) for k in ind_indices]
    if min(prop_inds) < 0.6: continue

    p_bar_count = 0
    freq_list = []

    for ref,alt,freq,ind in zip(ref_indices, alt_indices, freq_indices, ind_indices):
        pop_bases = "".join(sorted(x[ref]+x[alt]))
#        print(pop_bases == canonical_bases,canonical_bases,  pop_bases )
        if pop_bases == canonical_bases:
            if x[ref] == ref_base:
                freq_list.append( float( x[freq] ) )
            elif x[alt] == ref_base:
                freq_list.append( float( x[freq] ) )
        else :
            if x[ref] == ref_base:
                freq_list.append( 0 )
            else:
                pass
    if len(freq_list) < len(ref_indices):continue
    p_bar = sum(freq_list)/len(freq_list)

    if p_bar < 0.05: continue
    output.write( ",".join( [x[0], x[1]]+ list(map(str,freq_list)))+"\n" )
#    counter +=1
#    if counter == 100:break
output.close()
