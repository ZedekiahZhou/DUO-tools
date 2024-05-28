import polars as pl
from scipy import stats
import argparse
from statsmodels.stats.multitest import multipletests

# args
parser = argparse.ArgumentParser(description="call m6Am sites from AGcounts file")
parser.add_argument("-i", "--input", type=str, required=True, help="AGcounts file from pileup_reads5p.py")
parser.add_argument("-o", "--output", type=str, required=True, help="output file name")
parser.add_argument("--Acov", type=int, default=5, help='minimum A coverage for m6Am sites, default is 5')
parser.add_argument("--FDR", type=float, default=0.001, help="FDR cutoff, default is 0.001")
parser.add_argument("--Signal_Ratio", type=float, default=0.8, 
                    help="minimum ratio of signal reads (eg. reads with unconverted As less than 3), default is 0.8")
parser.add_argument("--AG_Ratio", type=float, default=0.8, 
                    help="minimum ratio of (A+G reads)/total in this sites, default is 0.8")
args = parser.parse_args()

# input and output
df=pl.read_csv(args.input, separator="\t")
frawout=args.input.replace("_AGcount.tsv", "") + "_m6Am_sites_raw.tsv"


# binomtest function
def test_sites(x):
    if x["AG_cov"] == 0:
        return 1
    else:
        return stats.binomtest(x["A"], n=x["AG_cov"], p=x["Ctrl_Ratio"], alternative="greater").pvalue


# perform test
df = df.with_columns(
    (pl.col("A") + pl.col("T") + pl.col("C") + pl.col("G")).alias("Signal_cov"),  # signal reads (with unconverted As <= 3?)
    (pl.col("A") + pl.col("G")).alias("AG_cov"),  # A + G coverage (signal)
    (pl.col("Next_pos_A") + pl.col("Next_pos_G")).alias("Next_pos_AG")  # Ctrl A + G coverage
).with_columns(
    (pl.col("Signal_cov")/pl.col("Counts")).round(3).alias("Signal_Ratio"),  # signal ratio
    (pl.col("AG_cov")/pl.col("Counts")).round(3).alias("AG_Ratio"),  # AG ratio
    pl.when(pl.col("Next_pos_AG") == 0)
    .then(pl.lit(0))
    .otherwise(pl.col("Next_pos_A")/pl.col("Next_pos_AG"))
    .round(3).alias("Ctrl_Ratio"),  # Ctrol unconverted ratio
    pl.when(pl.col("AG_cov") == 0)
    .then(pl.lit(0))
    .otherwise(pl.col("A")/pl.col("AG_cov"))
    .round(3).alias("m6Am_Ratio")  # m6Am ratio
).with_columns(
    pl.struct(["A", "AG_cov", "Ctrl_Ratio"])
    .map_elements(test_sites, return_dtype=pl.Float64)
    .alias("Pvalue")  # binom_test p value
).with_columns(
    pl.col("Pvalue")
    .map_batches(lambda x: multipletests(x, method="fdr_bh")[1])
    .alias("FDR")  # BH adjusted
)

df.write_csv(frawout, separator="\t")


# filtering
used_col = ['Chr', 'Pos', 'Strand', 'A', 'AG_cov', 'm6Am_Ratio', 'Pvalue', 'FDR', 
            'Ref_base', 'Counts', 'TPM', 'geneID', 'txID', 'txBiotype', 'Dist']
df2 = df.filter(
    (pl.col("A") >= args.Acov) & (pl.col("FDR") < args.FDR) & (pl.col("Signal_Ratio") > args.Signal_Ratio) & (pl.col("AG_Ratio") > args.AG_Ratio)
).select(pl.col(used_col))
df2.write_csv(args.output, separator="\t")