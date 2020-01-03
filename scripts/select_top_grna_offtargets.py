import pandas as pd


def select(infile, outfile, outfile_fa, num):
    df = pd.read_csv(
        infile, sep="\t", index_col=False, header=0, quotechar='"', doublequote=True
    )
    # get rid of rows wothout a efficacy value
    df["top10OfftargetTotalScore"] = pd.to_numeric(df["top10OfftargetTotalScore"], errors="coerce")
    #df = df.dropna()
    # sort
    df = df.sort_values(by=["top10OfftargetTotalScore", "gRNAefficacy"], ascending=[True, False])
    # select top num
    df = df.head(num)
    df.to_csv(outfile, sep="\t", index=False)
    with open(outfile_fa, "w") as out:
        for index, row in df.iterrows():
            out.write(">{}\n{}\n".format(row["names"], row["gRNAsPlusPAM"]))


select(
    snakemake.params["infile"],
    snakemake.output[0],
    snakemake.output[1],
    snakemake.params["num"]
)

