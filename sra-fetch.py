#!/usr/bin python
import sys
import subprocess
from optparse import OptionParser

# parse commandline arguments
parser = OptionParser()
parser.add_option("-s", "--species", dest="species", help="species to query", metavar="SPECIES")
(options, args) = parser.parse_args()
species = options.species

cmd = "wget -O ./acc_list.csv \'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=acclist&" \
      "term=" + '\"' + species + '\"' + "[Organism] AND " + '\"' + "cluster_public" + '\"' + "[prop] AND " + '\"' + \
      "strategy wgs" + '\"' + "[Properties]\'"

sub = subprocess.Popen(cmd, stdout=subprocess.PIPE)
output, error = sub.communicate()

