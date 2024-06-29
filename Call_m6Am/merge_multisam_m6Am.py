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
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, help="Input files recording TSS sites (*_TSS.bed)")
parser.add_argument("--prx", type=str, nargs='+', help="columns names (usually sample names) to used for each input")
parser.add_argument("-o", "--output", type=str, required=True, help="output file prefix")

parser.add_argument("-c", "--AGcov", type=int, default=15, help='minimum A+G coverage for TSS, default is 15')
parser.add_argument("--tpm", type=float, default=1.0, help="minimum TPM value for TSS, default is 1.0")
parser.add_argument("--absDist", type=int, default=1000, help="maximum absolute distance to any annotated TSS from GTF file, default is 1000")
parser.add_argument("--prop", type=float, default=0.05, help="minimum proportion relative to the total TPM of a gene, default is 0.05")
parser.add_argument("--zscore", type=float, default=1.0, help="minimum Z-score (calculated within a gene) for TSS, default is 1.0")
parser.add_argument("-C", "--Acov", type=int, default=5, help='minimum A coverage for m6Am sites, default is 5')
parser.add_argument("-s", "--Signal_Ratio", type=float, default=0.8, 
                    help="minimum ratio of signal reads (eg. reads with unconverted As less than 3), default is 0.8")
parser.add_argument("-adp", "--FDR", type=float, default=0.05, help="FDR cutoff, default is 0.05")
args = parser.parse_args()

output = {}
if args.prx is None:
    prx=[re.match("(.*/)?([^/]+)_TSS.bed$", s).group(2) for s in args.inputs]
else:
    prx=args.prx

# I. merge TSS sites
# read and merge file, keep sites passed in any file
print("Merge TSS sites: ", flush=True)
for i in range(len(args.inputs)):
    print("Read sample " + prx[i], flush=True)

    tss = pl.read_csv(args.inputs[i], separator="\t")
    tss = tss.select(
        pl.col("Chr", "End", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype"),
        pl.col("Counts").alias("Counts_"+prx[i]),
        pl.col("TPM").alias("TPM_"+prx[i]),
        ((pl.col("Counts") >= args.AGcov) & (pl.col("TPM") >= args.tpm) & (pl.col("absDist") <= args.absDist) & \
         (pl.col("relSum") > args.prop) & (pl.col("zscore") > args.zscore)).alias("Passed_"+prx[i])
    )

    if i == 0:
        tss_merged = tss
    else:
        tss_merged = tss_merged.join(
            tss, on=["Chr", "End", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype"], how="outer", coalesce=True
        )

# Keep sites passed in any sample
tss_merged = tss_merged.fill_null(False).with_columns(
    pl.col("End").alias("Pos"),
    pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=pl.col("^Passed_.*$")).alias("Passed")
).filter(
    pl.col("Passed")
)


tss_merged.write_csv(args.output + "_TSS_merged.tsv" , separator='\t')

# m6Am_merged = tss_merged.select(
#     pl.col("Chr", "End", "Strand", "Base", "geneID", "txID", "Dist")
# )
## II. merged m6Am sites
print("Merge m6Am sites: ", flush=True)
for i in range(len(args.inputs)):
    print("Read sample " + prx[i], flush=True)

    m6Am = pl.read_csv(args.inputs[i].replace("_TSS.bed", "_m6Am_sites.tsv"), separator="\t")
    if "Signal_Ratio" not in m6Am.columns:
        m6Am = m6Am.with_columns(pl.lit(1).alias("Signal_Ratio"))
    m6Am = m6Am.select(
        pl.col("Chr", "Pos", "Strand"),
        pl.col("AG_cov").alias("AGcov_"+prx[i]),
        pl.col("A").alias("Acov_"+prx[i]),
        ((pl.col("AG_cov") >= args.AGcov) & (pl.col("TPM") >= args.tpm) & (pl.col("Dist") <= args.absDist) & \
        (pl.col("Signal_Ratio") >= args.Signal_Ratio) & (pl.col("A") >= args.Acov) & \
        (pl.col("FDR") < args.FDR)).alias("Passed_"+prx[i])
    )

    if i == 0:
        m6Am_merged = m6Am
    else:
        m6Am_merged = m6Am_merged.join(
            m6Am, on=["Chr", "Pos", "Strand"], how="outer", coalesce=True
        )

# Keep sites passed in any sample
m6Am_merged = m6Am_merged.fill_null(False).with_columns(
    pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=pl.col("^Passed_.*$")).alias("Passed")
).filter(
    pl.col("Passed")
)

filter_m6Am_merged = tss_merged.select(
        pl.col("Chr", "Pos", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype")
).join(
    m6Am_merged, on = ["Chr", "Pos", "Strand"], how="inner"
)

filter_m6Am_merged.write_csv(args.output + "_m6Am_merged.tsv" , separator='\t')