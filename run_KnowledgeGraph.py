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


phage_phage_ntw = "phage_phage.ntw"
phage_host_ntw = "phage_host.ntw"




################################################################################
############################  Edge construction  ###############################
################################################################################

# Add phage-phage
G = nx.Graph()
with open(phage_phage_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0]
        node2 = tmp[1]
        G.add_edge(node1, node2, weight = 1)


# Add phage-host
with open(phage_host_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0]
        node2 = tmp[1]
        G.add_edge(node1, node2, weight = 1)


################################################################################
############################  Nodes construction ###############################
################################################################################

