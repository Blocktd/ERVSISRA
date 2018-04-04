#!/usr/bin python
import sys
import os.path
import subprocess
import random
import datetime
from optparse import OptionParser

# parse commandline arguments
parser = OptionParser()
parser.add_option("-sp", "--species", dest="species", help="species to query", metavar="SPECIES")
parser.add_option("-ss", "--sample-size", dest="sample_size", help="number of desired acc to sample", metavar="SAMPLESIZE")
(options, args) = parser.parse_args()
species = options.species
sample_size = options.sample_size
acc_filename = "acc_list_" + species + ".csv"
# remove underscore from species name
species = species.replace("_", " ")

# if a file of accession numbers already exists, don't waste resources and time redownloading
if not os.path.isfile(acc_filename):
    cmd = "wget -O" + acc_filename + " \'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=acclist&" \
          "term=" + '\"' + species + '\"' + "[Organism] AND " + '\"' + "cluster_public" + '\"' + "[prop] AND " + '\"' + \
          "strategy wgs" + '\"' + "[Properties]\'"

    sub = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output, error = sub.communicate()
else:
    print "acc_list file already exists, will not download another"

# read and sample 100 random accession numbers
with open(acc_filename, "rb") as acc_file:
    rand_accs = random.sample(acc_file.readlines(), sample_size)

# strip newlines
rand_accs = map(lambda s: s.strip(), rand_accs)
print rand_accs

# write a time-stamped file of random accession numbers
time_stamp = str(datetime.datetime.now().isoformat())
rand_acc_file = "rand_accs_" + time_stamp + ".txt"
with open(rand_acc_file, "wb") as write_file:
    for acc in rand_accs:
        write_file.write(acc)
        write_file.write("\n")

