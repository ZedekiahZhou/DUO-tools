#!/usr/bin/env python
"""
Author: Zhe Zhou, Peking University, Yi lab
Date: June 9, 2024
Email: zzhou24@pku.edu.cn
Program: Merge and filter m6A sites from multiple samples
"""

# remove xxx 

import argparse, re
import polars as pl

parser = argparse.ArgumentParser(description="Merge and filter TSS and m6Am sites from multiple samples")
parser.add_argument("-i", "--inputs", type=str, nargs='+', required=True, 
                    help="Input files recording raw TSS sites (*_TSS_raw.bed.annotated.rmdup)")
parser.add_argument("--untreated", action="store_true", 
                          help="Untreated samaples (only merge and filter TSS)")
parser.add_argument("--prx", type=str, nargs='+', help="columns names (usually sample names) to used for each input")
parser.add_argument("-o", "--output", type=str, required=True, help="output file prefix")

parser.add_argument("-c", "--AGcov", type=int, default=15, help='minimum A+G coverage for TSS (default: 15)')
parser.add_argument("--tpm", type=float, default=1.0, help="minimum TPM value for TSS (default: 1.0)")
parser.add_argument("--absDist", type=int, default=1000, 
                    help="maximum absolute distance to any annotated TSS from GTF file (default: 1000)")
parser.add_argument("--prop", type=float, default=0.05, 
                    help="minimum proportion relative to the total TPM of a gene (default: 0.05)")
parser.add_argument("--zscore", type=float, default=1.0, 
                    help="minimum Z-score (calculated within a gene) for TSS (default: 1.0)")
parser.add_argument("-C", "--Acov", type=int, default=5, 
                    help='minimum A coverage for m6Am sites (default: 5)')
parser.add_argument("-s", "--Signal_Ratio", type=float, default=0.8, 
                    help="minimum ratio of signal reads, eg. reads with unconverted As less than 3 (default: 0.8)")
parser.add_argument("-adp", "--FDR", type=float, default=0.05, help="FDR cutoff (default: 0.05)")
parser.add_argument("-R", "--AG_Ratio", type=float, default=0.8, 
                    help="minimum ratio of (A+G reads)/total in this sites (default: 0.8)")
args = parser.parse_args()


if args.prx is None:
    prx=[re.match("(.*/)?([^/]+)_TSS_raw.bed.annotated.rmdup$", s).group(2) for s in args.inputs]
else:
    prx=args.prx
maxv=999999999

# I. merge TSS sites
# read and merge file, keep sites passed in any file
print("Merge TSS sites: ", flush=True)
for i in range(len(args.inputs)):
    print("Read sample " + prx[i], flush=True)

    tss = pl.read_csv(args.inputs[i], separator="\t", null_values=".")
    tss = tss.fill_null(".").with_columns(
        pl.col("absDist").fill_null(maxv),
        pl.col("Dist").fill_null(maxv),
        pl.col("relSum").fill_null(-maxv),
        pl.col("zscore").fill_null(-maxv),
    ).select(
        pl.col("Chr", "End", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype"),
        pl.col("Counts").alias("Counts_"+prx[i]),
        pl.col("TPM").alias("TPM_"+prx[i]),
        ((pl.col("Counts") >= args.AGcov) & (pl.col("TPM") >= args.tpm) & (pl.col("absDist") <= args.absDist) & \
        (pl.col("relSum") >= args.prop) & (pl.col("zscore") >= args.zscore)).alias("Passed_"+prx[i])
    )

    if i == 0:
        tss_merged = tss
    else:
        tss_merged = tss_merged.join(
            tss,
            on=["Chr", "End", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype"], how="outer", coalesce=True
            #tss.select(pl.exclude(["geneID", "txID", "Dist", "txBiotype"])), 
            #on=["Chr", "End", "Strand", "Base"], how="outer", coalesce=True
        )

# Keep sites passed in any sample
tss_merged = tss_merged.fill_null(False).with_columns(
    pl.col("End").alias("Pos"),
    pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=pl.col("^Passed_.*$")).alias("Passed")
).filter(
    pl.col("Passed")
).sort(["Chr", "Pos"])

tss_merged.write_csv(args.output + "_TSS_merged.tsv" , separator='\t', null_value=".")


## II. merged m6Am sites
if not args.untreated:
    print("Merge m6Am sites: ", flush=True)
    for i in range(len(args.inputs)):
        print("Read sample " + prx[i], flush=True)

        m6Am = pl.read_csv(
            args.inputs[i].replace("_TSS_raw.bed.annotated.rmdup", "_m6Am_sites_raw.tsv"), 
            separator="\t", null_values="."
        )
        if "Signal_Ratio" not in m6Am.columns:
            m6Am = m6Am.with_columns(pl.lit(1).alias("Signal_Ratio"))

        m6Am = m6Am.join(
            tss_merged.select(
                pl.col("Chr", "Pos", "Strand", "Passed_"+prx[i])
            ), 
            on = ["Chr", "Pos", "Strand"], how="left"
        )

        m6Am=m6Am.fill_null(".").with_columns(
            pl.col("Dist").fill_null(maxv)
        ).select(
            pl.col("Chr", "Pos", "Strand"),
            pl.col("Ref_base").alias("Base"),
            pl.col("geneID", "txID", "Dist", "txBiotype"),
            pl.col("AG_cov").alias("AGcov_"+prx[i]),
            pl.col("A").alias("Acov_"+prx[i]),
            ((pl.col("Passed_"+prx[i])) & (pl.col("AG_cov") >= args.AGcov) & (pl.col("A") >= args.Acov) & \
            (pl.col("Signal_Ratio") >= args.Signal_Ratio) & (pl.col("AG_Ratio") >= args.AG_Ratio) & \
            (pl.col("FDR") < args.FDR))
        )

        if i == 0:
            m6Am_merged = m6Am
        else:
            m6Am_merged = m6Am_merged.join(
                m6Am,
                on=["Chr", "Pos", "Strand", "Base", "geneID", "txID", "Dist", "txBiotype"], how="outer", coalesce=True
                #m6Am.select(pl.exclude(["geneID", "txID", "Dist", "txBiotype"])), 
                #on=["Chr", "Pos", "Strand", "Base"], how="outer", coalesce=True
            )

    # Keep sites passed in any sample
    m6Am_merged = m6Am_merged.fill_null(False).with_columns(
        pl.fold(acc=pl.lit(False), function=lambda acc, x: acc | x, exprs=pl.col("^Passed_.*$")).alias("Passed")
    ).filter(
        pl.col("Passed")
    ).sort(["Chr", "Pos"])

    m6Am_merged.write_csv(args.output + "_m6Am_merged.tsv" , separator='\t', null_value=".")