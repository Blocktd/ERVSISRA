#!/usr/bin python
import os.path
import datetime
from Bio import Entrez
import argparse

parser = argparse.ArgumentParser(description="primary filtering of biosample data")
parser.add_argument('-sp', dest='sp', help="query species", required=True)
parser.add_argument('-bm', dest='bm', help="biomolecule type", choices=['dna', 'rna'], required=True)
parser.add_argument('-al', dest='al', help="return only aligned", choices=['yes','no'], required=True)
args = parser.parse_args()

# Second level of filterng that can only be done via the SRA
# Only want SRA IDs that fit criteria

# retrieve unfilitered SRA IDs from sra-id file
unfiltered = []
os.chdir("./id_archives")
for file in os.listdir('.'):
    if "sra-ids" in file:
        with open(file, "rb") as idfile:
            for line in idfile:
                line = line.replace("\n", "")
                unfiltered.append(line)

# print "unfiltered is of length %d" % len(unfiltered)
unfiltered_ids = set(unfiltered)
# print "the set is of length %d" % len(unfiltered_ids)

aligned = ""

# process commandline args
biomol = args.bm
if args.al == "yes":
    aligned = " AND sra_nuccore_alignment[Filter]"
    aln_param = "aligned"
else:
    aligned = ""
    aln_param = "unaligned"

# define species
species = args.sp
species = species.replace("_", " ")

# Perform esearch query
Entrez.email = "scr6045@g.rit.edu"
esearch_term = '\"' + species + '\"' + "[ORGANISM] AND " + '\"' + species + '\"' + "[orgn] AND cluster_public[prop]" \
               " AND biomol " + biomol + "[Properties]" + aligned

print "Quering with term:\n"
print esearch_term

esearch_query = Entrez.esearch(restart=0, db="sra", term=esearch_term, retmax=100000, idtype="acc")
esearch_result = Entrez.read(esearch_query)

# return id #s
sra_ids = esearch_result["IdList"]
# exhaustively search entire database (retmax is capped at 100000)
if len(sra_ids) == 100000:
    ret_len = 100000
    iter_ctr = 1
    while ret_len == 100000:
        esearch_query = Entrez.esearch(retstart=(100000*iter_ctr), db="sra",
                                       term=esearch_term, retmax=100000, idtype="acc")
        esearch_result = Entrez.read(esearch_query)
        temp_ids = esearch_result["IdList"]
        for i in temp_ids:
            sra_ids.append(i)
        ret_len = len(temp_ids)
        iter_ctr += 1   

filtered_ids = unfiltered_ids.intersection(sra_ids)

curr_time = datetime.datetime.now().isoformat()
filtered_fn = args.sp + "-" + aln_param + "-" + biomol + "-accessions-" + curr_time + ".txt"

# print information about results
print "\nFound %d hits in the SRA database for esearch query \n" % len(sra_ids)
print "Found %d valid accession numbers from secondary filter of SRA accession " % len(filtered_ids)
print "numbers using filter: %s" % esearch_term

# write filtered accession numbers to a file
print "\nWriting file of filtered acc..."
os.chdir('..')
with open(filtered_fn, 'wb') as ofn:
    for acc in filtered_ids:
        ofn.write(acc)
        ofn.write("\n")
ofn.close()
print "...done"

exit(0)
