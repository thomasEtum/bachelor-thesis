#!/usr/bin python3

import pickle
import networkx as nx
import csv
import time

starttime = time.time()
tf_target_dict = open("path/to/dict_tf_targets.pickle", "rb")
tf_target_dict = pickle.load(tf_target_dict)

G = nx.read_graphml("Cpath/to/G_ppi_reg_lrr.graphml")

# make cmap_drug -> uniprot dictionary or load from pickle
#
# drug_to_uni = csv.reader(open("path/to/STITCH_TTD_GTP_BDB_uniprot_drugbank_chembl.txt"), delimiter="\t")
# drug_uni_dict = dict()
# for line in drug_to_uni:
#     uni = line[0]
#     drug = line[1]
#     if drug in drug_uni_dict:
#         drug_uni_dict[drug].append(uni)
#     else:
#         drug_uni_dict[drug] = [uni]
#
# path_to_name_dict = "path/to/dict_cmap_name_drugbank_atc_chembl.pickle" # dictionary for cmap names of drugs
#
# drug_name_dict = pickle.load(open(path_to_name_dict, "rb"))
#
# big_dict = dict()
# for drug in drug_name_dict:
#     names = drug_name_dict[drug]
#     all_uni = []
#     if "drugbank_id" in names:
#         db = names["drugbank_id"]
#         db = db.strip("drugbank.")
#         if db in drug_uni_dict:
#             all_uni.extend(drug_uni_dict[db])
#     if "chembl_id" in names:
#         chembl = names["chembl_id"]
#         if chembl in drug_uni_dict:
#             all_uni.extend(drug_uni_dict[chembl])
#     # atc is a list
#     if "atc" in names:
#         atc = names["atc"]
#         for i in atc:
#             if i in drug_uni_dict:
#                 all_uni.extend(drug_uni_dict[i])
#     if all_uni:
#         big_dict[drug] = all_uni
#
# pickle.dump(big_dict, open("path/to/new/drugname_uniprot_dict.pickle", "wb"))

big_dict = pickle.load(open("C:\\Users\\thoma\\Documents\\Uni\\bachelor\\gsea\\network_proximity\\drugname_uniprot_dict.pickle", "rb"))

# s = candidate tf, T = [associated proteins], t = T[i]
# d(s, T) = 1/len(T)*sum(d(s, t))
# d(s, t) = number of edges between s and t

pairs_to_test = open("path/to/pairs/of/drug/and/tf.tsv")
pair_reader = csv.reader(pairs_to_test, delimiter="\t")

def checksum(G, big_dict, drug, uniIDls, uniID):
    targets = big_dict[drug]
    sumPath = 0
    checker = False
    for t in targets:
        t = "uniprot." + t
        if t in G:
            lenls = []
            for i in uniIDls:
                if i in G:
                    if i == uniID:
                        lenls.append(nx.shortest_path_length(G, i, t))
                    else:
                        lenls.append(nx.shortest_path_length(G, i, t)+1)
            sumPath += min(lenls)
            checker = True
    if checker:
        d = 1/len(targets)*sumPath
    else:
        d = False
    return d

dists = []

for line in pair_reader:
    drug = line[0]
    origuniID = line[2]
    uniID = "uniprot." + origuniID
    if drug in big_dict:
        if uniID in tf_target_dict:
            uniIDls = tf_target_dict[uniID]
            uniIDls.append(uniID)
            uniIDls = list(set(uniIDls))
        else:
            uniIDls = [uniID]
        checker = False
        for i in uniIDls:
            if i in G:
                checker = True
                break
        if checker:
            distance = checksum(G, big_dict, drug, uniIDls, uniID)
        if distance:
            dists.append([drug, origuniID, distance])


result_writer = csv.writer(open("path/to/network/proximity/file.tsv", "w", newline=''), delimiter="\t")
for result in dists:
    result_writer.writerow(result)
