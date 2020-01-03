import pandas as pd


def select(infile, outfile, outfile_fa, num, minimum):
    df = pd.read_csv(
        infile, sep="\t", index_col=False, header=0, quotechar='"', doublequote=True
    )
    # get rid of rows wothout a efficacy value
    df["gRNAefficacy"] = pd.to_numeric(df["gRNAefficacy"], errors="coerce")
    df = df.dropna()
    # sort
    df = df.sort_values(by="gRNAefficacy", ascending=False)
    # subselect based on value of efficacy
    df = df.loc[df["gRNAefficacy"] >= minimum]
    # select top num
    df = df.head(num)
    df.to_csv(outfile, sep="\t", index=False)

    with open(outfile_fa, "w") as out:
        for index, row in df.iterrows():
            out.write(">{}\n{}\n".format(row["name"], row["gRNAplusPAM"]))


select(
    snakemake.params["in_eff"],
    snakemake.output[0],
    snakemake.output[1],
    snakemake.params["num"],
    snakemake.params["min"],
)

