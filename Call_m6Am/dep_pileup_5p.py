from Bio import SeqIO
import pandas as pd
import argparse
import pysam
import pandas as pd
from Bio.Seq import reverse_complement

def get_refer_base(key):
    global reference_genome
    chr, pos, strand = key
    if strand == "+":
        return reference_genome[chr][pos - 1]
    else:
        return reverse_complement(reference_genome[chr][pos -1])


if __name__ == "__main__":
    description = """
    Pileup 5' ends of reads:

    1. reads with soft clip in 5' were excluded
    2. reads with too many unconverted As were excluded
    3. use next A as control
    """

    parser = argparse.ArgumentParser(prog="pileup_reads_5p",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-r","--ref",dest="reference",required=True,help="reference fasta")
    group_required.add_argument("-l","--list",dest="fsite",required=True,help="site list file")
    group_required.add_argument("-b","--bam",dest="fbam",required=True,help="input bam, sorted")
    group_required.add_argument("-o","--output",dest="output",required=True,help="output")
    # Optional
    group_optional = parser.add_argument_group("Optional")
    group_optional.add_argument("--maxAs", dest="max_allowed_As", default=3,
                            help="maximum allowed number of As in reads to be count as signal ones, default=3")
    options = parser.parse_args()

    # Step1: Init
    ## parse reference
    reference_genome = {}
    for seq in SeqIO.parse(options.reference,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    ## init site dict
    output = {}
    with open(options.fsite, "r") as fsite:
        for line in fsite.readlines():
            chr, pos, strand = line.strip().split("\t")
            pos = int(pos)
            ID = (chr, pos, strand)

            output[ID] = {}
            output[ID]["Ref_base"] = get_refer_base(ID)
            output[ID]["A"] = 0
            output[ID]["T"] = 0
            output[ID]["C"] = 0
            output[ID]["G"] = 0
            output[ID]["n_raw"] = 0
    
    ## find the position of next A 
    with open(options.fsite, "r") as fsite:
        for line in fsite.readlines():
            chr, pos, strand = line.strip().split("\t")
            pos = int(pos)
            ID = (chr, pos, strand)

            # find next A
            next_pos = pos

            if strand == "+":
                while True:
                    next_pos += 1
                    key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "A":
                        output[ID]["Next_pos"] = next_pos
                        break
            else:
                while True:
                    next_pos -= 1
                    key = (chr, next_pos, strand)
                    # if key not in output:
                    next_base = reference_genome[chr][next_pos-1]
                    if next_base == "T":
                        output[ID]["Next_pos"] = next_pos
                        break

            output[ID]["Next_ref_base"] = get_refer_base((chr, next_pos, strand))
            output[ID]["Next_pos_A"] = 0
            output[ID]["Next_pos_T"] = 0
            output[ID]["Next_pos_C"] = 0
            output[ID]["Next_pos_G"] = 0

    # Step II: Pileup 
    with pysam.AlignmentFile(options.fbam, "rb") as input_BAM:
        # iterate by reads, pileup 5' end
        for read in input_BAM:
            aligned_pairs = read.get_aligned_pairs()
            aligned_query_pos = [i[0] for i in aligned_pairs]
            aligned_ref_pos = [i[1] for i in aligned_pairs] # 0-based reference pos

            # get aligned positon of the 5' end of reads
            if read.is_forward:
                pos = aligned_ref_pos[0]
                query_base = read.query_sequence[0]
                A_counts = read.query_sequence.count('A') # number of un-converted As in reads
            else:
                pos = aligned_ref_pos[-1]
                query_base = reverse_complement(read.query_sequence[-1])
                A_counts = read.query_sequence.count('T')
        
            chr = read.reference_name.replace("_AG_converted", "") # remove suffix
            strand = '+' if read.is_forward else '-'

            if pos is not None: # exclude 5' soft clipped reads
                pos = pos + 1 # convert to 1-based reference pos
                ID = (chr, pos, strand)

                if ID in output and query_base in ("A", "T", "C", "G"):
                    output[ID]["n_raw"] += 1
                    if A_counts <= options.max_allowed_As: # only count signal reads (un-converted As < max_allowed_As)
                        output[ID][query_base] += 1

                        # if next A pos is covered by this read, then count
                        next_pos = output[ID]["Next_pos"]
                        try: 
                            next_query_idx = aligned_ref_pos.index(next_pos-1)
                        except ValueError:
                            next_query_idx = None
                        
                        if next_query_idx is not None:
                            next_query_pos = aligned_query_pos[next_query_idx]
                            if next_query_pos is not None:
                                next_query_base = read.query_sequence[next_query_pos]
                                if read.is_reverse:
                                    next_query_base = reverse_complement(next_query_base)
                                if next_query_base in ("A", "T", "C", "G"):
                                    output[ID]["Next_pos_"+next_query_base] += 1
                    

df = pd.DataFrame.from_dict(output, orient='index')
df.to_csv(options.output)
        