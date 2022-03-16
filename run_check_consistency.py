import pandas as pd
import numpy as np
import math
from ete3 import NCBITaxa
ncbi = NCBITaxa()


def return_taxa_label(taxa, item):
    label = 'unknown'
    try:
        name2taxid = ncbi.get_name_translator([item])
        taxid = list(name2taxid.values())[0][0]
        lineage = ncbi.get_lineage(taxid)
        for item in lineage:
            rank_id = ncbi.get_rank([item])
            if list(rank_id.values())[0] == taxa:
                label_dict = ncbi.get_taxid_translator([list(rank_id.keys())[0]])
                label = list(label_dict.values())[0]
    except:
        pass
    return label




HostG = pd.read_csv("final_prediction.csv")
Phylum = HostG['phylum'].values
Class  = HostG['class'].values
Order  = HostG['order'].values
Family = HostG['family'].values
Genus  = HostG['genus'].values
accession = HostG['contig_name'].values
tmp_name = HostG['idx'].values
new_Phylum = []
new_Class  = []
new_Order  = []
new_Family = []
new_Genus  = []

for phy, cla, orde, fam, gen in zip(Phylum, Class, Order, Family, Genus):
    if phy != return_taxa_label('phylum', str(cla)):
        new_Phylum.append(phy)
        new_Class.append('-')
        new_Order.append('-')
        new_Family.append('-')
        new_Genus.append('-')
    elif cla != return_taxa_label('class', str(orde)):
        new_Phylum.append(phy)
        new_Class.append(cla)
        new_Order.append('-')
        new_Family.append('-')
        new_Genus.append('-')
    elif orde != return_taxa_label('order', str(fam)):
        new_Phylum.append(phy)
        new_Class.append(cla)
        new_Order.append(orde)
        new_Family.append('-')
        new_Genus.append('-')
    elif fam != return_taxa_label('family', str(gen)):
        new_Phylum.append(phy)
        new_Class.append(cla)
        new_Order.append(orde)
        new_Family.append(fam)
        new_Genus.append('-')
    else:
        new_Phylum.append(phy)
        new_Class.append(cla)
        new_Order.append(orde)
        new_Family.append(fam)
        new_Genus.append(gen)




new_pred = pd.DataFrame({'contig_name': accession, 'idx': tmp_name, 'Phylum':new_Phylum, 'Class': new_Class, 'Order': new_Order, "Family": new_Family, "Genus": new_Genus})
new_pred.to_csv("final_prediction_consistency.csv", index=None)
