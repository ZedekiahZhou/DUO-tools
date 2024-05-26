import polars as pl
import sys, argparse

parser = argparse.ArgumentParser(description="Annotate transcription start sites (TSSs)")
parser.add_argument("-i", "--input", type=str, required=True, help="Input TSS annotation files from bedtools intersect")
parser.add_argument("--cov", type=int, default=15, help='minimum A+G coverage for TSS, default is 15')
parser.add_argument("--tpm", type=float, default=1.0, help="minimum TPM value for TSS, default is 1.0")
parser.add_argument("--absDist", type=int, default=1000, help="maximum absolute distance to any annotated TSS from GTF file, default is 1000")
parser.add_argument("--zscore", type=float, default=1.0, help="minimum Z-score (calculated within a gene) for TSS, default is 1.0")
args = parser.parse_args()


priorities = {"snRNA":0, "snoRNA":1, "protein_coding":2, "miRNA":3, "lncRNA": 4}

def main():
    global args, fout, priorities
    frmdup=args.input + ".rmdup"
    fout=args.input.replace("_TSS_raw.bed.annotated", "_TSS.bed")

    df = pl.read_csv(args.input, separator="\t", has_header=False)
    colname = ["Chr", "Start", "End", "ID", "Counts", "Strand", "Base", "TPM", "DBname",
            "txChr", "txStart", "txEnd", "geneID", "txID", "txStrand", "txTSS", 
            "geneBiotype", "txBiotype"]
    df.columns = colname if df.shape[1] == 18 else colname[:8] + colname[9:]

    def get_pri(s: str) -> int:
        pri = priorities.get(s)
        if pri is None:
            pri = 99
        return pri

    # remove duplicate annotations
    df_rmdup = df.with_columns(
        pl.col("txBiotype").map_elements(get_pri, return_dtype=pl.Int64).alias("Priorities")
    ).with_columns(
        pl.when(pl.col("Strand") == "+")
        .then(pl.col("End")-pl.col("txTSS"))
        .otherwise(pl.col("txTSS")-pl.col("End"))
        .alias("Dist")
    ).with_columns(
        pl.col("Dist").abs().alias("absDist")
    ).group_by("ID").agg(
        pl.col("*").sort_by("absDist", "Priorities").first()
    ).sort("ID")

    df_rmdup = df_rmdup.with_columns(
        pl.col("TPM").sum().over("geneID").alias("totalTPM"),
        pl.col("TPM").count().over("geneID").alias("nTSS")
    ).with_columns(
        (pl.col("TPM")/pl.col("totalTPM")).alias("relSum"),
        pl.when(pl.col("TPM").count() <= 2)
        .then(pl.lit(99))
        .otherwise(((pl.col("TPM") - pl.col("TPM").mean()) / pl.col("TPM").std()))
        .over("geneID")
        .alias("zscore")
    )

    df_rmdup.write_csv(frmdup, separator="\t")

    # filter TSS
    df_clean = df_rmdup.filter(
        (pl.col("Counts") >= args.cov) & (pl.col("TPM") >= args.tpm) & (pl.col("absDist") <= args.absDist) & (pl.col("zscore") >= args.zscore)
    )
    df_clean.write_csv(fout, separator="\t")

if __name__ == '__main__':
    main()