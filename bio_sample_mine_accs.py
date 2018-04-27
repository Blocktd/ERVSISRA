#!/usr/bin python

import os
import datetime
from Bio import Entrez
from time import sleep
import tarfile

# simplify calling of Entrez.post, for use with batch calling of 1000 ids 
# later in script
# returns a tuple of webenv and query_key
def post_to_ncbi(ids, database, history):
    epost_query = Entrez.epost(db=database, id=",".join(ids), 
                               method="post", history=history)
    epost_results = Entrez.read(epost_query)
    webenv = epost_results["WebEnv"]
    query_key = epost_results["QueryKey"]
    return webenv, query_key


# user-friendly wrapper for Entrez.elink
def elink_from_history(webenv, query_key):
    elink_query = Entrez.elink(dbfrom="BioSample", db="sra", webenv=webenv, 
                               query_key=query_key, idtype="acc", retmax=1000)
    elink_result = Entrez.read(elink_query)
    sra_ids = [link["Id"] for link in elink_result[0]["LinkSetDb"][0]["Link"]]
    return sra_ids


# takes a list of BioSample ids and returns corresponding SRA ids
def get_sra_from_bio_samp(id_list):
    post_info = post_to_ncbi(id_list, "BioSample", "y")
    webenv = post_info[0]
    query_key = post_info[1]
    sra_ids = elink_from_history(webenv, query_key)
    return sra_ids


# Split a list into multiple lists of size n
def split_list(list, size):
    for i in xrange(0, len(list), size):
        yield list[i:i + size]


# performs elink querying on a very large set of fromdb uids
def perf_large_elink(fromdb_id_list):
    fail_ctr = 0
    success_ctr = 0
    sra_ids_final = []
    # split list into chunks of 1k ids
    chunked_list = split_list(fromdb_id_list, 1000)
    # attempt to elink each chunk, if fails once, retries once and only once
    for chunk in chunked_list:
        print "Retrieving chunk..."
        try:
            sra_ids_final.append(get_sra_from_bio_samp(chunk))
            success_ctr += 1
        except:
            print "chunk download failed, retrying..."
            try:
                sra_ids_final.append(get_sra_from_bio_samp(chunk))
                print "download recovered"
                success_ctr += 1
            except:
                print "download unrecoverable. Will continue, data will be \
                       incomplete"
                fail_ctr += 1
        sleep(.5)
        print "...done"
    print "successfully retrieved %d chunks, failed to retrieve %d chunks" % \
    (success_ctr, fail_ctr)
    return sra_ids_final


# archives ids for future use and version ctl
def archive_data(biosamp_data, sra_data, dir, raw, species):
    # get current datetime in ISO format
    curr_time = datetime.datetime.now().isoformat()
    os.chdir("./id_archives")
    bio_samp_arch_fn = "biosamp-ids-" + curr_time + ".txt"
    sra_id_arch_fn = "sra-ids-" + curr_time + ".txt"
    print bio_samp_arch_fn
    print sra_id_arch_fn

    # write text file of biosample ids
    with open(bio_samp_arch_fn, "wb") as ba:
        for id in biosamp_data:
            ba.write(id)
            ba.write("\n")
    ba.close()

    # write text file of sra ids
    with open(sra_id_arch_fn, "wb") as sa:
        for id in sra_data:
            sa.write(id)
            sa.write("\n")
    sa.close()

    # create a tar archive of the data and if del_og == true, 
    # delete non-archived files
    tar_fn = species +"-data-archive-" + curr_time + ".tar.gz"
    # tar interpets : in a strange and undesirable manner
    tar_fn = tar_fn.replace(":", "-")
    tarchive = tarfile.open(tar_fn, 'w:gz')
    for file in os.listdir(os.getcwd()):
        if file.endswith(".txt"):
            tarchive.add(file)
            if not raw:
                os.remove(file)
    tarchive.close()


# top-level filtering command, ONLY THIS SHOULD BE CALLED IN OTHER FILES
def filter(species, raw, email):
    # esearch entrez util for:
    # find specified ids with disease attribute, public, sra linked
    Entrez.email = email
    esearch_term = '\"' + species + '\"' + "[ORGANISM] AND public[filter] \
    AND attribute disease[filter] AND biosample sra[filter]"
    print esearch_term
    esearch_query = Entrez.esearch(db="BioSample", term=esearch_term, 
                                   retmax=100000, idtype="BioSampleId")
    esearch_result = Entrez.read(esearch_query)
    
    biosamp_ids = esearch_result["IdList"]
    
    # Retrieve all linked SRA archives from biosample ids
    sra_ids = perf_large_elink(biosamp_ids)
    # un nest sra ids
    sra_ids = [val for sublist in sra_ids for val in sublist]
    
    # Generate an archive of the data
    if not os.path.isdir("id_archives"):
        os.mkdir("id_archives")
        archive_data(biosamp_ids, sra_ids, "id_archives", raw=raw, 
                     species=species)
    else:
        archive_data(biosamp_ids, sra_ids, "id_archives", raw=raw, 
                     species=species)
   
    print "...done filtering."
    return(0)
