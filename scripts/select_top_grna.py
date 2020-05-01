import sys, csv
import pandas as pd


def get_loc_strand(infile_loc, location):
    strand = "+"
    with open(infile_loc, "r") as fin:
        d = {}
        rdr = csv.reader(fin, delimiter="\t")
        for a in rdr:
            if a[0] == location:
                strand = a[2]
    return strand


def main(infile, infile_loc, outfile, outfile_fa, num, minimum, same_strand, location):
    df = pd.read_csv(
        infile, sep="\t", index_col=False, header=0, quotechar='"', doublequote=True
    )
    # get rid of rows wothout a efficacy value
    df["gRNAefficacy"] = pd.to_numeric(df["gRNAefficacy"], errors="coerce")
    df = df.dropna()

    # if same_strand = True, retain only desired strand
    if same_strand:
        strand = get_loc_strand(infile_loc, location)
        df = df.loc[df["strand"] == strand]

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


main(
    snakemake.params["in_eff"],
    snakemake.input[1],
    snakemake.output[0],
    snakemake.output[1],
    snakemake.params["num"],
    snakemake.params["min"],
    snakemake.params["same_strand"],
    snakemake.params["location"],
)
