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
parser.add_argument("-adp", "--FDR", type=float, default=0.05, help="FDR cutoff, default is 0.05")


args = parser.parse_args()
if args.prx is None:
    prx=[re.match("(.*/)?([^/]+).totalm6A.FDR.csv$", s).group(2) for s in args.inputs]
else:
    prx=args.prx


# read and merge file, keep sites passed in any file
for i in range(len(args.inputs)):
    print("Read sample " + prx[i], flush=True)

    df = pl.read_csv(args.inputs[i], separator="\t")
    if "Signal_Ratio" not in df.columns:
        df = df.with_columns(pl.lit(1).alias("Signal_Ratio"))
    df = df.select(
        pl.col("Chr", "Sites", "Strand", "Gene"),
        pl.col("AGcov").alias("AGcov_"+prx[i]),
        pl.col("Acov").alias("Acov_"+prx[i]),
        ((pl.col("AGcov") >= args.AGcov) & (pl.col("Acov") >= args.Acov) & (pl.col("Signal_Ratio") >= args.Signal_Ratio) & \
        (pl.col("Ratio") >= args.methyl_Ratio) & (pl.col("P_adjust") < args.FDR)).alias("Passed_"+prx[i])
    )

    if i == 0:
        df_merged = df
    else:
        df_merged = df_merged.join(
            df, on=["Chr", "Sites", "Strand", "Gene"], how="outer", coalesce=True
        )

# keep sites passed in any sample
df_merged = df_merged.fill_null(False).with_columns(
    pl.col("Sites").alias("Pos"),
    pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=pl.col("^Passed_.*$")).alias("Passed")
)

df_merged.write_csv(args.output + ".raw", separator='\t')
df_merged.filter(pl.col("Passed")).write_csv(args.output, separator='\t')
