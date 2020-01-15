import sys
import os, os.path
import csv
import pandas as pd
from collections import OrderedDict



def stats(infiles, info_file, bed_file, outfile):

    outfile = open(outfile, "w")
    reader = csv.reader(open(bed_file, "r"), delimiter="\t")
    dtx2loc = {}
    for a in reader:
        assert a[3] not in dtx2loc
        dtx2loc[a[3]] = "chr{}:{}-{},{}".format(a[0], a[1], a[2], a[5])
        
    reader = csv.reader(open(info_file, "r"), delimiter="\t")
    d = OrderedDict()
    for a in reader:
        d[a[0]] = a

    a_num_grna = []
    a_num_avg_ot = []
    for f in infiles:
        loc = os.path.basename(f).split(".")[0]

        tx = d[loc][1]
        try:
            location = dtx2loc[tx]
        except KeyError:
            sys.stderr.write("Tx not found. EXIT.")

        d[loc] += [location]

        df = pd.read_csv(f, sep="\t", index_col=False, header=0)
        df["offtarget_num"] = pd.to_numeric(df["offtarget_num"], errors="coerce")

        # number without ot
        num_zero_ot = df[df["offtarget_num"] == 0].shape[0]

        # with ot
        y = df[df["offtarget_num"] > 0]
        num_with_ot = y.shape[0]
        avg_num_ot = y["offtarget_num"].mean()

        d[loc] += [str(num_zero_ot), str(num_with_ot), str(avg_num_ot)]

    out_header = "loc\ttx\tregion_searched\tnum_grna\tnum_grna_without_offtargets\tnum_grna_with_offtargets\tavg_num_offtargets\n"

    outfile.write(out_header)
    for k in d:
        outfile.write("{}\n".format("\t".join(d[k])))
    outfile.close()
        

stats(snakemake.input["infiles"],
      snakemake.input["info"],
      snakemake.input["bed"],
      snakemake.output[0])

