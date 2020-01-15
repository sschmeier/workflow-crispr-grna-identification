import sys
import os, os.path
import csv
import gzip
import bz2
import zipfile


def load_file(filename):
    """ LOADING FILES """
    if filename in ["-", "stdin"]:
        filehandle = sys.stdin
    elif filename.split(".")[-1] == "gz":
        filehandle = gzip.open(filename, "rt")
    elif filename.split(".")[-1] == "bz2":
        filehandle = bz2.open(filename, "rt")
    elif filename.split(".")[-1] == "zip":
        filehandle = zipfile.ZipFile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def combine(infiles, stats_file, outfile):
    outfile = open(outfile, "w")
    
    reader = csv.reader(load_file(stats_file), delimiter="\t")
    d = {}
    for a in reader:
        d[a[0]] = a[0:3]

    out_header = "loc\ttx\tregion_searched\trank\tgRNAplusPAM\tname\tstart\tstrand\textendedSequence\tgRNAefficacy\tofftarget_num\tofftarget_max_score\tofftarget_max_efficacy\n"
    outfile.write(out_header)
    for f in infiles:
        rank = 1
        loc = os.path.basename(f).split(".")[0]
        assert loc in d
        tmp = d[loc]
        reader = csv.reader(load_file(f), delimiter="\t")
        head = next(reader)
        for a in reader:
            res = tmp + [str(rank)] + a
            outfile.write("{}\n".format("\t".join(res)))
            rank += 1
    # close res file        
    outfile.close()
            

combine(snakemake.input["infiles"],
        snakemake.input["stats"],
        snakemake.output[0])
