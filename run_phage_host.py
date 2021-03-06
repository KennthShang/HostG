import os
import sys
import Bio
import logging
import argparse
import subprocess
import scipy as sp
import numpy as np
import pandas as pd
import pickle as pkl
import networkx as nx
import scipy.stats as stats
import scipy.sparse as sparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Defined folder
bacteria_in = "bacteria/"
blast_database_out = "blast_db/"
blast_tab_out = "blast_tab/"
Knowledge_graph = "Cyber_data/"

################################################################################
############################  Check the folder #################################
################################################################################
def check_folder(file_name):
    if not os.path.exists(file_name):
        _ = os.makedirs(file_name)
    else:
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name):
            try:
                _ = subprocess.check_call("rm -rf {0}".format(file_name), shell=True)
                _ = os.makedirs(file_name)
                print("Dictionary cleaned")
            except:
                print("Cannot clean your folder... permission denied")
                exit(1)

check_folder(blast_database_out)
check_folder(blast_tab_out)



#  combine phage file 
_ = subprocess.check_call("cat dataset/nucl.fasta single_contig/* > out/query.fa", shell=True)



################################################################################
###############################  Run BLASTN   ##################################
################################################################################


genome_list = os.listdir(bacteria_in)
for genome in genome_list:
    make_blast_cmd = 'makeblastdb -in '+ bacteria_in + genome +' -dbtype nucl -parse_seqids -out '+ blast_database_out + genome.split(".")[0]
    print("Creating blast database...")
    _ = subprocess.check_call(make_blast_cmd, shell=True)
    blast_cmd = 'blastn -query out/query.fa -db blast_db/'+genome.split(".")[0]+' -outfmt 6 -out '+ blast_tab_out + genome.split(".")[0]+'.tab -num_threads 8'
    print("Running blastn...")
    _ = subprocess.check_call(blast_cmd, shell=True)


################################################################################
###############################  phage-host   ##################################
################################################################################


tab_file_list = os.listdir(blast_tab_out)
phage_set = set()
bacteria2phage = {}
for file in tab_file_list:
    bacteria_id = file.split('.')[0]
    phage_id_list = []
    with open(blast_tab_out+file) as file_in:
        for line in file_in.readlines():
            tmp = line.split('\t')
            phage_id = tmp[0]
            try:
                bacteria2phage[bacteria_id].append(phage_id)
            except:
                bacteria2phage[bacteria_id] = [phage_id]

# De-duplication
for key in bacteria2phage:
    bacteria2phage[key] = list(set(bacteria2phage[key]))


# Save the phage-host graph
with open("out/phage_host.ntw", 'w') as file_out:
    for bacteria in bacteria2phage:
        for phage in bacteria2phage[bacteria]:
            _ = file_out.write(bacteria + "," + phage + "\n")
