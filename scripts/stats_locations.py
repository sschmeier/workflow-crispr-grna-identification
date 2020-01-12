#!/usr/bin/env python
"""
NAME: split_seqs.py
===================

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.1    2019    Initial version.

LICENCE
=======
2019, copyright Sebastian Schmeier
s.schmeier@pm.me // https://www.sschmeier.com

template version: 2.0 (2018/12/19)
"""
__version__ = "0.0.1"
__date__ = "2019"
__email__ = "s.schmeier@pm.me"
__author__ = "Sebastian Schmeier"

import sys
import os
import argparse
import csv
import gzip
import bz2
import zipfile
import time
import re
import pandas as pd
from collections import OrderedDict

# non-standard lib: For color handling on the shell
try:
    from colorama import init, Fore

    # INIT color
    # Initialise colours for multi-platform support.
    init()
    reset = Fore.RESET
    colors = {
        "success": Fore.GREEN,
        "error": Fore.RED,
        "warning": Fore.YELLOW,
        "info": "",
    }
except ImportError:
    sys.stderr.write(
        "colorama lib desirable. " + 'Install with "conda install colorama".\n\n'
    )
    reset = ""
    colors = {"success": "", "error": "", "warning": "", "info": ""}


def alert(atype, text, log, repeat=False, flush=False):
    if repeat:
        textout = "{} [{}] {}\r".format(
            time.strftime("%Y%m%d-%H:%M:%S"), atype.rjust(7), text
        )
    else:
        textout = "{} [{}] {}\n".format(
            time.strftime("%Y%m%d-%H:%M:%S"), atype.rjust(7), text
        )

    log.write("{}{}{}".format(colors[atype], textout, reset))
    if flush:
        log.flush()
    if atype == "error":
        sys.exit(1)


def success(text, log=sys.stderr, flush=True):
    alert("success", text, log, flush=flush)


def error(text, log=sys.stderr, flush=True):
    alert("error", text, log, flush=flush)


def warning(text, log=sys.stderr, flush=True):
    alert("warning", text, log, flush=flush)


def info(text, log=sys.stderr, repeat=False, flush=True):
    alert("info", text, log, repeat=repeat, flush=flush)


def parse_cmdline():
    """ Parse command-line args. """
    # parse cmd-line ----------------------------------------------------------
    description = "Read files, make stats"
    version = "version {}, date {}".format(__version__, __date__)
    epilog = "Copyright {} ({})".format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("--version", action="version", version="{}".format(version))

    parser.add_argument("str_file_bed", metavar="BED-FILE", help="Bed-file")

    parser.add_argument("str_file_info", metavar="INFO-FILE", help="INFO-file")

    parser.add_argument("infiles", nargs="+", metavar="INFILE", help="Input files.")

    parser.add_argument(
        "-o", "--outfile", dest="outfile_name", metavar="OUTFILE", help="Outfile."
    )

    # if no arguments supplied print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args, parser


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


def stats(infiles, info_file, bed_file, outfile):

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
            error("Tx not found. EXIT.")

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


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ["-", "stdout"]:
        outfileobj = sys.stdout
    elif args.outfile_name.split(".")[-1] == "gz":
        outfileobj = gzip.open(args.outfile_name, "wt")
    elif args.outfile_name.split(".")[-1] == "bz2":
        outfileobj = bz2.BZ2File(args.outfile_name, "wt")
    else:
        outfileobj = open(args.outfile_name, "w")

    stats(args.infiles, args.str_file_info, args.str_file_bed, outfileobj)
    # ------------------------------------------------------

    return


if __name__ == "__main__":
    sys.exit(main())
