#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: Dec 27, 2024
Email: zzhou24@pku.edu.cn
Program: Merge replicates into one sample
Recording the mean CR, mean Genecov, mean Signal_Ratio, best Pvalue and P_ajust
used polars version: 0.20.25 --> 1.5.0
"""

import argparse, re, time
import polars as pl

parser = argparse.ArgumentParser(description="Merge and filter m6A sites from multiple samples")
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, 
                    help="Input files recording m6A sites (*totalm6A.FDR.csv)")
parser.add_argument("-o", "--output", type=str, required=True, 
                    help="output file prefix, output file will be named as prefix.totalm6A.FDR.csv")
args = parser.parse_args()


# read and merge data --------
print("[%s] Pooling replicates ========" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
for i in range(len(args.inputs)):
    print("\t[%s] Read sample %s ========" % 
        (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), args.inputs[i]), 
        flush=True)

    df = pl.read_csv(args.inputs[i], separator="\t")
    if "Signal_Ratio" not in df.columns:
        df = df.with_columns(pl.lit(1).alias("Signal_Ratio"))

    # outer join samples, Ratio not updated
    if i == 0:
        df_merged = df
        used_cols = df.columns
        # Chr     Sites   Strand  Gene    CR      AGcov   Acov    Genecov Signal_Ratio    Ratio   Pvalue  P_adjust
    else:
        df_merged = df_merged.join(
            df, suffix="_right",
            on=["Chr", "Sites", "Strand", "Gene"], 
            how="full", coalesce=True 
        ).with_columns(
            CR=pl.mean_horizontal("CR", "CR_right"),
            AGcov=pl.sum_horizontal("AGcov", "AGcov_right"),
            Acov=pl.sum_horizontal("Acov", "Acov_right"),
            Genecov=pl.mean_horizontal("Genecov", "Genecov_right"),
            Signal_Ratio=pl.mean_horizontal("Signal_Ratio", "Signal_Ratio_right"),
            Pvalue=pl.min_horizontal("Pvalue", "Pvalue_right"),
            P_adjust=pl.min_horizontal("P_adjust", "P_adjust_right"),
        ).select(used_cols)

df_merged = df_merged.with_columns(
    Ratio=(pl.col("Acov") / pl.col("AGcov")).round(3)
).select(used_cols)
df_merged = df_merged.sort(["Chr", "Sites"])
df_merged.write_csv(args.output+".totalm6A.FDR.csv", separator="\t", null_value=".")
