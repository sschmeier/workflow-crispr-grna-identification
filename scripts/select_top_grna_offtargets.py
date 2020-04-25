# 1. Select all gRNA withput offtargets, sort by best efficacy
# 2. Select from those the top number
# 3. If number cannot be filled as here are not enough gRNAs
#    withoput offtargets, do
#
#    - Calculate number of offtargets for the ones with offtargets
#    - Select the ones with minimum offtargets
#    - stop if we have the requred number selected or running out of gRNAs
import sys, csv
import pandas as pd


def select(infile, infile_grna, outfile, num):
    df = pd.read_csv(
        infile, sep="\t", index_col=False, header=0, quotechar='"', doublequote=True
    )

    outfile = open(outfile, "w")
    # need to extract gRNAs from gRNA file in order of res
    reader = csv.reader(open(infile_grna, "r"), delimiter="\t")
    header = next(reader)
    outfile.write(
        "{}\tofftarget_num\tofftarget_max_score\tofftarget_max_efficacy\n".format(
            "\t".join(header)
        )
    )

    # if not gRNAefficacy column found, analysis did not work correct, add info re this
    if "gRNAefficacy" not in df.columns or "top10OfftargetTotalScore" not in df.columns:
        for a in reader:
            outfile.write(
                "{}\t{}\t{}\t{}\n".format(
                    "\t".join(a),
                    "offtarget_analysis_unsuccessful",
                    "offtarget_analysis_unsuccessful",
                    "offtarget_analysis_unsuccessful",
                )
            )

    else:
        # column to numeric replace NaN with 0 = > NaN are good!!
        df["top10OfftargetTotalScore"] = pd.to_numeric(
            df["top10OfftargetTotalScore"], errors="coerce"
        ).fillna(0)

        # subselext gRNA names + offtarget + efficiancy value
        df = df[["names", "top10OfftargetTotalScore", "gRNAefficacy"]]

        # no: no_offtargetsm => sort by based efficacy
        no = df.loc[df["top10OfftargetTotalScore"] == 0]
        no = no.sort_values(by=["gRNAefficacy"], ascending=[False])
        no_rows = no.shape[0]
        num_still_to_get = num - no_rows
        res = list(no["names"].head(num))

        res2 = []
        # if we cannot fill the num with gRNAs from the pool without off targets
        if num_still_to_get > 0:
            # sort the gRNAs with of targets and take top
            # criteria?

            # wo: with_offtargets
            wo = df.loc[df["top10OfftargetTotalScore"] != 0]
            # number of off target sites
            wo_ot_counts = (
                wo.groupby("names")["top10OfftargetTotalScore"].count().to_frame()
            )
            wo_ot_counts.columns = ["ot_count"]
            wo_ot_counts = wo_ot_counts.sort_values(by="ot_count", ascending=True)
            res2 = list(wo_ot_counts.index)[0:num_still_to_get]

            # also get maxima for top10OfftargetTotalScore and efficacy
            wo_ot_max = wo.groupby("names").max()

        d = {}
        for a in reader:
            name = a[1]
            assert name not in d
            d[name] = a

        # iterate over results
        for name in res:
            outfile.write("{}\t0\tNA\tNA\n".format("\t".join(d[name])))

        # iterate over results2
        for name in res2:
            num_ot = wo_ot_counts.loc[name]["ot_count"]
            max_ot_score = wo_ot_max.loc[name]["top10OfftargetTotalScore"]
            max_ot_eff = wo_ot_max.loc[name]["gRNAefficacy"]
            outfile.write(
                "{}\t{}\t{}\t{}\n".format(
                    "\t".join(d[name]), num_ot, max_ot_score, max_ot_eff
                )
            )


select(
    snakemake.params["infile_summary"],
    snakemake.input["grna_file"],
    snakemake.output[0],
    snakemake.params["num"],
)
