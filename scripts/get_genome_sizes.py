#!/usr/bin/env python
"""
NAME: get_genome_sizes.py
=========================

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
    description = "Read gtf file and extract a set of unique 5 prime coordinates bed-style. Multiple files can be pipe'd into this script with (z)cat."
    version = "version {}, date {}".format(__version__, __date__)
    epilog = "Copyright {} ({})".format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("--version", action="version", version="{}".format(version))

    parser.add_argument("str_file", metavar="GENOME-FASTA-FILE", help="Delimited file.")

    parser.add_argument(
        "-o",
        "--out",
        metavar="STRING",
        dest="outfile_name",
        default=None,
        help='Out-file. [default: "stdout"]',
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


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    try:
        fileobj = load_file(args.str_file)
    except IOError:
        error('Could not load file "{}". EXIT.'.format(args.str_file))

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

    outfileobj.write("chrom\tsize\n")
    from Bio import SeqIO

    for record in SeqIO.parse(fileobj, "fasta"):
        outfileobj.write("{}\t{}\n".format(record.id, len(record.seq)))

    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == "__main__":
    sys.exit(main())
