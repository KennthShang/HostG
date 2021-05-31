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
from cluster import make_diamond_db, run_diamond
from cluster import make_protein_clusters_mcl, load_mcl_clusters, build_clusters
from cluster import build_pc_matrices, to_clusterer, create_network


# Defined folder
out_fn = "out/"
contig_in = "input/"
contig_out = "single_contig/"
file_in_fn = "single_contig/"
file_out_fnn = "all_proteins/"
Knowledge_graph = "Cyber_data/"

################################################################################
############################  Check the folder #################################
################################################################################
if not os.path.exists(out_fn):
    _ = os.makedirs(out_fn)
else:
    print("folder {0} exist... cleaning dictionary".format(out_fn))
    if os.listdir(out_fn):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(out_fn), shell=True)
            _ = os.makedirs(out_fn)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists(file_in_fn):
    _ = os.makedirs(file_in_fn)
else:
    print("folder {0} exist... cleaning dictionary".format(file_in_fn))
    if os.listdir(file_in_fn):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(file_in_fn), shell=True)
            _ = os.makedirs(file_in_fn)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists(file_out_fnn):
    _ = os.makedirs(file_out_fnn)
else:
    print("folder {0} exist... cleaning dictionary".format(file_out_fnn))
    if os.listdir(file_out_fnn):
        try:
            _ = subprocess.check_call("rm -rf {0}".format(file_out_fnn), shell=True)
            _ = os.makedirs(file_out_fnn)
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)



################################################################################
############################  Rename the files  ################################
################################################################################

# Split contigs files into single contig per file.
file_list = sorted(os.listdir(contig_in))
seq = []
old_file_id = 0
contig_id = 0
with open("name_list.csv",'w') as list_out:
    list_out.write("contig_name,idx\n")
    for file_n in file_list:
        for record in SeqIO.parse(contig_in+file_n, "fasta"):
            name = str(old_file_id) + "_" + str(contig_id)
            contig_id += 1
            list_out.write(record.id + "," + name + "\n")
            _ = SeqIO.write(record, contig_out+name+".fasta", "fasta")
        old_file_id += 1


################################################################################
###################### Translate contigs into 6 ORFs ###########################
################################################################################

def return_protien(contig):
    frame = Bio.Seq.translate(contig)
    protein_list = frame.split("*")
    sorted_list = sorted(protein_list, key=len, reverse=True)
    sorted_list = [item for item in sorted_list if len(item) > 10 ]
    return sorted_list


file_list = os.listdir(file_in_fn)
for file in file_list:
    old_file_name = file.rsplit(".", 1)[0]
    contig_id = int(old_file_name.split("_")[-1])
    label_id = int(old_file_name.split("_")[0])
    for record in SeqIO.parse(file_in_fn+file, "fasta"):
        protein_record = []
        contig = str(record.seq)
        # conver into protein
        frame1 = return_protien(contig)
        frame2 = return_protien(contig[1:])
        frame3 = return_protien(contig[2:])
        rev_contig = Bio.Seq.reverse_complement(contig)
        frame4 = return_protien(rev_contig)
        frame5 = return_protien(rev_contig[1:])
        frame6 = return_protien(rev_contig[2:])
        proteins = np.concatenate([frame1, frame2, frame3, frame4, frame5, frame6])
        for i in range(len(proteins)):
            rec = SeqRecord(Seq(proteins[i]), id=str(label_id)+ "_" + str(contig_id) + "_" + str(i), description="")
            protein_record.append(rec)
        _ = SeqIO.write(protein_record, file_out_fnn+file, "fasta")

all_protein_f = out_fn+"all_translate_proteins.fa"
_ = subprocess.check_call("cat {0} > {1}".format(file_out_fnn+"*", all_protein_f), shell=True)


################################################################################
############################## Run diamond BLASTp  #############################
################################################################################

print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
print("Creating Diamond database and running Diamond...")




proteins_aa_fp = all_protein_f
db_fp = "database/database.dmnd"
diamond_out_fnn = '{}.diamond.tab'.format(os.path.basename(proteins_aa_fp).rsplit('.', 1)[0])
diamond_out_fnp = os.path.join(out_fn, diamond_out_fnn)
# Create database
_ = make_diamond_db(8)
# Run BLASTP
similarity_fp = run_diamond(proteins_aa_fp, db_fp, 8, diamond_out_fnp)

# capture the query, referencde, e-value from diamond output
contig_abc_fp = out_fn + diamond_out_fnn + ".abc"
abc_fp = out_fn+"merged.abc"
_ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fnp, contig_abc_fp), shell=True)
_ = subprocess.check_call("cat database/database.self-diamond.tab.abc {0} > {1}".format(contig_abc_fp, abc_fp), shell=True)

# Generating gene-to-genome.csv: protein_id, contig_id, keywords
blastp = pd.read_csv(contig_abc_fp, sep=' ', names=["contig", "ref", "e-value"])
protein_id = sorted(list(set(blastp["contig"].values)))
contig_id = [item.rsplit("_", 1)[0] for item in protein_id]
description = ["hypothetical protein" for item in protein_id]
gene2genome = pd.DataFrame({"protein_id": protein_id, "contig_id": contig_id ,"keywords": description})
gene2genome.to_csv(out_fn+"contig_gene_to_genome.csv", index=None)

# Counting for the mapped proteins
mapped_record = []
for record in SeqIO.parse(all_protein_f, "fasta"):
    if record.id in protein_id:
        mapped_record.append(record)
        
# Save the mapped proteins
mapped_fp = out_fn+"mapped_protein.fa"
with open(mapped_fp, 'w') as file_out:
    _ = SeqIO.write(mapped_record, file_out, "fasta")


# Combining the gene-to-genomes files
_ = subprocess.check_call("cat database/database_gene_to_genomes.csv {0}contig_gene_to_genome.csv > {1}gene_to_genome.csv".format(out_fn, out_fn), shell=True)



################################################################################
################################## Run MCL #####################################
################################################################################
print("\n\n" + "{:-^80}".format("Protein clustering"))
print("Loading proteins...")
gene2genome_fp = out_fn+"gene_to_genome.csv"
gene2genome_df = pd.read_csv(gene2genome_fp, sep=',', header=0)


# Parameters for MCL
pc_overlap, pc_penalty, pc_haircut, pc_inflation = 0.8, 2.0, 0.1, 2.0
pcs_fp = make_protein_clusters_mcl(abc_fp, out_fn, pc_inflation)
print("Building the cluster and profiles (this may take some time...)")


# Dump MCL results
protein_df, clusters_df, profiles_df, contigs_df = build_clusters(pcs_fp, gene2genome_df)
print("Saving files")
dfs = [gene2genome_df, contigs_df, clusters_df]
names = ['proteins', 'contigs', 'pcs']
output_dir = out_fn

for name, df in zip(names, dfs):
    fn = "Cyber_{}.csv".format(name)
    fp = os.path.join(output_dir, fn)
    index_id = name.strip('s') + '_id'
    if not os.path.exists(fp):
        df.set_index(index_id).to_csv(fp)
    else:
        print("File {} exists and will be used. Use -f to overwrite.".format(fn))

profiles_fn = "Cyber_profiles.csv"
profiles_fp = os.path.join(out_fn, profiles_fn)
if not os.path.exists(profiles_fp):
    profiles_df.to_csv(profiles_fp, index=False)
else:
    print("File {} exists and will be used. Use -f to overwrite.".format(profiles_fn))




################################################################################
################################## Run MCL #####################################
################################################################################

# Loding dataset
contigs_df = pd.read_csv("out/Cyber_contigs.csv")
clusters_df = pd.read_csv("out/Cyber_pcs.csv")
profiles_df = pd.read_csv("out/Cyber_profiles.csv")

# Replace names
contigs_csv_df = contigs_df.copy()
print("Read {} entries from {}".format(len(contigs_csv_df), os.path.join(output_dir, '{}_contigs.csv'.format(name))))
contigs_csv_df.index.name = "pos"
contigs_csv_df.reset_index(inplace=True)

pcs_csv_df = clusters_df.copy()
profiles = profiles_df.copy()

# Filtering the PC profiles that appears only once
before_filter = len(profiles)
cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

# get the number of contigs for each pcs and add it to the dataframe
cont_by_pc.columns = ["pc_id", "nb_proteins"]
pcs_csv_df = pd.merge(pcs_csv_df, cont_by_pc, left_on="pc_id", right_on="pc_id", how="left")
pcs_csv_df.fillna({"nb_proteins": 0}, inplace=True)

# Drop the pcs that <= 1 contig from the profiles.
pcs_csv_df = pcs_csv_df[pcs_csv_df['nb_proteins'] > 1]  # .query("nb_contigs>1")
at_least_a_cont = cont_by_pc[cont_by_pc['nb_proteins'] > 1]  # cont_by_pc.query("nb_contigs>1")
profiles = profiles[profiles['pc_id'].isin(at_least_a_cont.pc_id)]
print("Read {} entries (dropped {} singletons) from {}".format(len(profiles), (before_filter - len(profiles)), profiles_fp))
pcs_csv_df = pcs_csv_df.reset_index(drop=True)
pcs_csv_df.index.name = "pos"
pcs_csv_df = pcs_csv_df.reset_index()

matrix, singletons = build_pc_matrices(profiles, contigs_csv_df, pcs_csv_df)
profiles_csv = {"matrix": matrix, "singletons": singletons}
merged_df = contigs_csv_df
merged_fp = os.path.join(output_dir, 'merged_df.csv')
merged_df.to_csv(merged_fp)

ntw = create_network(matrix, singletons, thres=1, max_sig=300)
fi = to_clusterer(ntw, out_fn+"intermediate.ntw", merged_df.copy())




################################################################################
############################### Dump the graph #################################
################################################################################

G = nx.Graph()
# Create graph
with open(out_fn+"intermediate.ntw") as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(" ")
        node1 = tmp[0]
        node2 = tmp[1]
        G.add_edge(node1, node2, weight = 1)


graph = "phage_phage.ntw"
with open(graph, 'w') as file_out:
    for node1 in G.nodes():
        for node2 in G.edges(node):
            file_out.write(node1+","+node2)