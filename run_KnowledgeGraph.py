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




################################################################################
###############################  Input Params  #################################
################################################################################

phage_phage_ntw = "out/phage_phage.ntw"
phage_host_ntw = "out/phage_host.ntw"
parser = argparse.ArgumentParser(description='manual to this script')
inputs = parser.parse_args()

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

phage2id = pkl.load(open("dataset/phage2id.dict",'rb'))
contig2id = pkl.load(open("Cyber_data/contig2id.dict",'rb'))
id2phageF = pkl.load(open("dataset/phage.F",'rb'))
id2contigF = pkl.load(open("Cyber_data/contig.F",'rb'))

node_feature = []
for node in G.nodes():
    # if query node
    if node in contig2id.keys():
        node_feature.append(id2contigF[contig2id[node]])
    # if phage node
    elif node in phage2id.keys():
        node_feature.append(id2phageF[phage2id[node]])
    # if bacteria node
    else:
        for _, neighbor in G.edges(node):
            cnt = 0
            tmpF = 0
            # if neighbor is query node
            if neighbor in contig2id.keys():
                tmpF += id2contigF[contig2id[neighbor]]
            # if neighbor is phage node
            elif neighbor in phage2id.keys():
                tmpF += id2phageF[phage2id[neighbor]]
            cnt+=1
        node_feature.append(tmpF/cnt)

node_feature = np.array(node_feature)

################################################################################
############################  Label construction ###############################
################################################################################

def taxa_label(taxa):
    cnt = 0
    label_df = pd.read_csv("dataset/label.csv")
    node_label = []
    for node in G.nodes():
        # if query node
        if node in contig2id.keys():
            neighbor_label = []
            for _, neighbor in G.edges(node):
                if neighbor not in contig2id.keys():
                    tmpL = label_df[label_df['accession'] == neighbor.split('.')[0]][taxa].values[0]
                    neighbor_label.append(tmpL)
            if len(set(neighbor_label)) == 1:
                node_label.append(neighbor_label[0])
                cnt += 1
            else:
                node_label.append(-1)
        # if phage or host node
        else:
            tmpL = label_df[label_df['accession'] == node.split('.')[0]][taxa].values[0]
            node_label.append(tmpL)
    # Convert string label to int label
    label_set = list(set(node_label))
    label2int = {tempL:tmpid for tmpid, tempL in enumerate(label_set)}
    # Record the test_id
    idx = 0
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
    # store features for GCN
    id2node = {idx: node for idx, node in enumerate(G.nodes())}
    adj = nx.adjacency_matrix(G)
    pkl.dump(adj, open("GCN_data/"+taxa+"_contig.graph", "wb" ))
    pkl.dump(test_id, open("GCN_data/"+taxa+"_contig.test_id", "wb" ))
    pkl.dump(node_feature, open("GCN_data/"+taxa+"_contig.feature", "wb" ))
    pkl.dump(node_int_label, open("GCN_data/"+taxa+"_contig.label", "wb" ))
    pkl.dump(id2node, open("GCN_data/"+taxa+"_id2node.dict", "wb" ))
    pkl.dump(label2int, open("GCN_data/"+taxa+"_label2int.dict", "wb" ))


taxa_list = ["phylum", 'class', 'order', 'family', 'genus']
for taxa in taxa_list:
    taxa_label(taxa)

