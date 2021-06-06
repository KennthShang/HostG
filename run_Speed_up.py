import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
import pandas as pd
import subprocess
import argparse
import re


#####################################################################
##########################  Input Params  ###########################
#####################################################################

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'test_contigs.fa')
parser.add_argument('--len', type=int, default=8000)
parser.add_argument('--gpus', type=int, default = 0)
parser.add_argument('--t', type=float, default=0.0)
inputs = parser.parse_args()


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

check_folder("input")
check_folder("pred")
check_folder("Split_files")



#####################################################################
#########################  Start Program  ###########################
#####################################################################

def special_match(strg, search=re.compile(r'[^ACGT]').search):
    return not bool(search(strg))


cnt = 0
file_id = 0
records = []
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    if cnt !=0 and cnt%1000 == 0:
        SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    seq = seq.upper()
    if special_match(seq):
        if len(record.seq) > inputs.len:
            records.append(record)
            cnt+=1

SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta")
file_id+=1

for i in range(file_id):
    cmd = "mv Split_files/contig_"+str(i)+".fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Moving file Error for file {0}".format("contig_"+str(i)))
        continue

    cmd = "python run_CNN.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Pre-trained CNN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue
        

    cmd = "python run_phage_phage.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("phage_phage Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = "python run_phage_host.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("phage_host Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = "python run_KnowledgeGraph.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Knowledge Graph Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = "python run_GCN_ECE.py --t "  + inputs.t + " --gpus " + inputs.gpus
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("GCN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    # Clean files
    cmd = "rm input/*"
    out = subprocess.check_call(cmd, shell=True)



    # combine taxa label
    phylum_df = pd.read_csv("tmp_pred/prediction_phylum.csv")
    class_df = pd.read_csv("tmp_pred/prediction_class.csv")
    order_df = pd.read_csv("tmp_pred/prediction_order.csv")
    family_df = pd.read_csv("tmp_pred/prediction_family.csv")
    genus_df = pd.read_csv("tmp_pred/prediction_genus.csv")

    tmp_pred = pd.merge(phylum_df, class_df, on='contig_names')
    tmp_pred = pd.merge(tmp_pred, order_df, on='contig_names')
    tmp_pred = pd.merge(tmp_pred, family_df, on='contig_names')
    tmp_pred = pd.merge(tmp_pred, genus_df, on='contig_names')


    name_list = pd.read_csv("name_list.csv")
    prediction = tmp_pred.rename(columns={'contig_names':'idx'})
    contig_to_pred = pd.merge(name_list, prediction, on='idx')
    contig_to_pred.to_csv("pred/contig_"+str(i)+".csv", index = None)

    cmd = "rm name_list.csv prediction.csv"
    out = subprocess.check_call(cmd, shell=True)

    cmd = "rm tmp_pred/*"
    out = subprocess.check_call(cmd, shell=True)


cmd = "cat pred/* > final_prediction.csv"
out = subprocess.check_call(cmd, shell=True)

