import polars as pl
import sys, argparse

parser = argparse.ArgumentParser(description="Annotate transcription start sites (TSSs)")
parser.add_argument("-i", "--input", type=str, required=True, help="Input TSS annotation files from bedtools intersect")
parser.add_argument("-o", "--output", type=str, required=False, help="output file name, default is input.rmdup")
args = parser.parse_args()


priorities = {"snRNA":0, "snoRNA":1, "protein_coding":2, "miRNA":3, "lncRNA": 4}

def main():
    global args, priorities
    if args.output is None:
        args.output=args.input + ".rmdup"

    df = pl.read_csv(args.input, separator="\t", has_header=False, null_values=".")
    colname = ["Chr", "Start", "End", "ID", "Counts", "Strand", "Base", "TPM", "DBname",
            "txChr", "txStart", "txEnd", "geneID", "txID", "txStrand", "txTSS", 
            "geneBiotype", "txBiotype"]
    df.columns = colname if df.shape[1] == 18 else colname[:8] + colname[9:]
    if "DBname" in df.columns:
        df = df.select(pl.col("*").exclude("DBname"))

    def get_pri(s: str) -> int:
        pri = priorities.get(s)
        if pri is None:
            pri = 99
        return pri

    # remove duplicate annotations
    df_rmdup = df.with_columns(
        pl.col("txBiotype").map_elements(get_pri, return_dtype=pl.Int64).alias("Priorities"),
        pl.when(pl.col("Strand") == "+")
        .then(pl.col("End")-pl.col("txTSS"))
        .otherwise(pl.col("txTSS")-pl.col("End"))
        .alias("Dist")
    ).with_columns(
        pl.col("Dist").abs().alias("absDist")
    ).group_by("ID").agg(
        pl.col("*").sort_by("absDist", "Priorities").first()
    ).sort("ID")

    # group sum, mean, and std
    df_agg=df_rmdup.group_by("geneID").agg(
        pl.col("TPM").sum().alias("raw_totalTPM"),
        pl.col("TPM").count().alias("raw_nTSS"),
        pl.col("TPM").mean().alias("meanTPM"),
        pl.col("TPM").std().alias("stdTPM")
    )

    df_rmdup = df_rmdup.join(df_agg, on="geneID", how="left").with_columns(
        (pl.col("TPM")/pl.col("raw_totalTPM")).round(3).alias("relSum"),
        pl.when((pl.col("raw_nTSS") <= 3))
        .then(pl.lit(99))
        .when(pl.col("stdTPM") == 0)
        .then(pl.lit(0))
        .otherwise(((pl.col("TPM") - pl.col("meanTPM")) / pl.col("stdTPM")))
        .round(3).alias("zscore"), 
        pl.col("raw_totalTPM").round(3),
        pl.col("meanTPM").round(3),
        pl.col("stdTPM").round(3)
    )

    df_rmdup.write_csv(args.output, separator="\t",null_value=".")


if __name__ == '__main__':
    main()