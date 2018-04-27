#!usr/bin/ python

import os
import datetime
import argparse
import bio_sample_mine_accs
import sra_filter

''' 
TODO This file will act as a wrapper for the data-analysis pipeline.
All data passed from the ruby interface will be stored and managed here, then
passed into the python and R pipelines
'''

# parse commandline arguments
arg_parser = argparse.ArgumentParser("Master Wrapper for  \
                                     input to analysis script")
# species
arg_parser.add_argument("--species", required=True, 
                        dest='spec', metavar='SPECIES',
                        help="species to query biosample for")
# retain raw 
arg_parser.add_argument("--raw", required=True, dest="raw",metavar="RAW",
                        help="Specify whether or not to keep an uncompressed \
                        list of biosample ids, on top of an archived version",
                        choices=["yes", "no"])
# aligned, unaligned, or both
arg_parser.add_argument("--aligned", required=True, dest="alnd",
                        metavar="ALIGNED", help="Specify whether the sra \
                        filtering should filter for just aligned or aligned \
                        genomes, or both", choices=["aligned", "unaligned",
                        "both"])
# biomolecule (RNA or DNA)
arg_parser.add_argument("--biomol", required=True,
                        dest="biomol", metavar="BIOMOL",
                        help="Specify whether to filter for just RNA or DNA",
                        choices=["DNA", "RNA"])
arg_parser.add_argument("--email", required=True,
                        dest="email", metavar="EMAIL", help="Email is \
                        required for purposes of job-logging and querying \
                        the ncbi via entrez utils")

# process commandline arguments
args = arg_parser.parse_args()
spec = args.spec
raw = args.raw
alnd = args.alnd
biomol = args.biomol
email = args.email

if raw == 'yes':
    raw = True
else:
    raw = False

# search biosample database    
bio_sample_mine_accs.filter(spec, raw, email)

# biosample IDs -> filtered SRA ids
sra_filter.filter(spec, alnd, biomol)

