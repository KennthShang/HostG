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

# required contig2id: Dict type
# required phage2id: Dict type
# required id2contigF: List type
# required id2phageF: List type
node_feature = []
for node in G.nodes():
    # if query node
    if node in contig2id.keys():
        node_feature.append(id2contigF[contig2id[node]])
    # if phage node
    elif node in phage2id.keys():
        node_feature.append(id2pahgeF[phage2id[node]])
    # if bacteria node
    else:
        for neighbor in G.edges(nodes):
            cnt = 0
            tmpF = 0
            # if neighbor is query node
            if neighbor in contig2id.keys():
                tmpF += id2contigF[contig2id[node]]
            # if neighbor is phage node
            elif neighbor in phage2id.keys():
                tmpF += id2pahgeF[phage2id[node]]
            cnt+=1
        node_feature.append(tmpF/cnt)


################################################################################
############################  Label construction ###############################
################################################################################


# required Label.csv (label_df): DataFrame type
# required taxa: String type
node_label = []
for node in G.nodes():
    # if query node
    if node in contig2id.keys():
        neighbor_label = []
        for _, neighbor in G.edges(node):
            if neighbor not in contig2id.keys():
                tmpL = label_df[label_df['accession'] == node][neighbor].values[0]
                neighbor_label.append(tmpL)

        if len(set(neighbor_label)) == 1:
            node_label.append(neighbor_label[0])
        else:
            node_label.append(-1)

    # if phage or host node
    else:
        tmpL = label_df[label_df['accession'] == node][taxa].values[0]
        node_label.append(tmpL)


# Convert string label to int label
label_set = list(set(node_label))
label2int = {tmpL:tmpid for tempL, tmpid in enumerate(label_set)}


test_id = []
node_int_label = []
for label in node_label:
    # labeled nodes
    if label != -1:
        node_int_label.append(label2int[label])
    else:
        node_int_label.append(-1)
        test_id.append(idx)
    idx+=1