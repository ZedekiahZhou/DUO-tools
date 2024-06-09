#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: June 9, 2024
Email: zzhou24@pku.edu.cn
Program: Merge and filter m6A sites from multiple samples
"""

import argparse, re
import pandas as pd
import polars as pl

parser = argparse.ArgumentParser(description="Merge and filter m6A sites from multiple samples")
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, help="Input files recording m6A sites (*totalm6A.FDR.csv)")
parser.add_argument("--prx", type=str, nargs='+', help="columns names (usually sample names) to used for each input")
parser.add_argument("-o", "--output", type=str, required=True, help="output file name")
parser.add_argument("-c", "--AGcov", type=int, default=15, help='minimum A+G coverage for TSS, default is 15')
parser.add_argument("-C", "--Acov", type=int, default=5, help='minimum A coverage for m6Am sites, default is 5')
parser.add_argument("-s", "--Signal_Ratio", type=float, default=0.8, 
                    help="minimum ratio of signal reads (eg. reads with unconverted As less than 3), default is 0.8")
parser.add_argument("-r", "--methyl_Ratio", type = float, default=0.1, help="minimum m6A level, default is 0.1")
parser.add_argument("-adp", "--FDR", type=float, default=0.001, help="FDR cutoff, default is 0.001")

args = parser.parse_args()

output = {}
if args.prx is None:
    prx=[re.match("(.*/)?([^/]+).totalm6A.FDR.csv$", s).group(2) for s in args.inputs]
else:
    prx=args.prx

for i in range(len(args.inputs)):
    print("Read sample " + prx[i], flush=True)
    with open(args.inputs[i]) as file:
        line = file.readline()  # skip first row which is column headers
        line = file.readline()
        while line:
            line_item = line.strip().split("\t")
            if len(line_item) == 12:
                chr, pos, strand, gene, geneCR, AGcov, Acov, Genecov, Signal_Ratio, ratio, Pvalue, p_adjust = line.strip().split("\t")
            elif len(line_item) == 11:
                chr, pos, strand, gene, geneCR, AGcov, Acov, Genecov, ratio, Pvalue, p_adjust = line.strip().split("\t")
                Signal_Ratio = 1
            else:
                raise ValueError("Input file format not correct!")
            ID = (chr, pos, strand)
            if ID not in output:
                    output[ID] = {"Chr": chr, "Sites": pos, "Strand": strand, "Gene": gene, "True_Any": 0}
            
            output[ID]["AGcov_" + prx[i]] = AGcov
            output[ID]["Acov_" + prx[i]] = Acov
            
            if int(AGcov) >= args.AGcov and int(Acov) >= args.Acov and float(Signal_Ratio) >= args.Signal_Ratio and \
                float(ratio) >= args.methyl_Ratio and float(p_adjust) < args.FDR:
                output[ID]["True_" + prx[i]] = int(1)
            else:
                output[ID]["True_" + prx[i]] = int(0)
                
            output[ID]["True_Any"] = int(output[ID]["True_Any"] or output[ID]["True_" + prx[i]])
            line = file.readline()

df = pl.from_pandas(pd.DataFrame.from_dict(output, orient='index'))
df = df.filter(pl.col("True_Any") == 1)
df.write_csv(args.output, separator='\t')
