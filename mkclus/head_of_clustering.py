# -------------------------------------------------------------------------------------------
"""
Input Files Required:
1.  repin_with_1k_flanks (from prepare_datasets.py)
    Path = bank_path/repin_with_1k_flanks.p
2.  blastgenomedb (from prepare_datasets.py)
    Path = bank_path/genomes_blastdb/allgenomes

Output Files:
1.  File with REPIN clusters
    Path = ./output/clusters_{date}

setup()
    Sets up repin_with_1k flanks and blastgenomedb

setup_flank_matches()
    For each flanking region for each REPIN, it performs a blast search and clusters
    REPINs of similar flanking region
    closest_repin_criteria indicates how close/far a REPIN must/can be from a flanking gene
    Ex: On blasting left flank of A, we get sequences for B, C, D
        We identify which repins are nearby and if they are to the left or right of this sequence
        Then cluster them as:
            mixed_clusters [A_left, B_left, C_left, D_right ....]
    This is equivalent to the mixed clustering I was previously doing in CD-HIT

prepare_repin_dict(mixed_clusters)
    Takes the mixed_clusters and establishes
        repin_dict = {rep: {'L': -1, 'R': -1} for rep in clusters}
    and left_clusters and right_clusters, which help group REPINs together in the next step
    Clusters are made in cases where left and right flank both match

flankclusterer(clusters)
    Middle Gene Deletion Case is identified and prevented
    Cluster B is merged with Cluster A IF:
        Left flank of B = Left flank of A
            Right flank of B DOESNT match with any other cluster
    OR:
        Right flank of B = Right flank of A
            Left flank of B DOESNT match with any other cluster

"""
# -------------------------------------------------------------------------------------------
import os
import time
import pickle
import networkx as nx
from copy import deepcopy
from itertools import combinations as itercomb
# import numpy as np
# filepath = "/".join(__file__.split("/")[:-1])
# os.chdir(filepath)
from mkclus import prepare_datasets

todaysdate = time.strftime("%b%d")
genome_blastdb = "genomes_blastdb/allgenomes"
temp_files = "dumpyard/"

genomes_list = []
repins_with_1k_flanks = {}
clusters = []
repin_dict, leftclus, rightclus = {}, {}, {}
repin_names = []
repins_per_genome = {}
closest_repin_criteria = 500
fast_mode_testin_only = False
all_parameters = {}
flank_gene_param = {}
allreplength = []
logging_mergers = []
logging_repin_dict = []
logging_trans = []


def progress_bar(current, total, bar_length=70):
    percent = float(current) * 100 / total
    arrow = '#' * int(percent / 100 * bar_length - 1)
    spaces = ' ' * (bar_length - len(arrow))
    if current == total:
        print('Progress: [%s%s] %d %%' % (arrow, spaces, percent))
    else:
        print('Progress: [%s%s] %d %%' % (arrow, spaces, percent), end='\r')


def prepare_repin_dict(mixed_clusters):
    global repin_dict, leftclus, rightclus, clusters, logging_mergers, logging_repin_dict
    repin_dict = {rep: {'L': -1, 'R': -1} for rep in repin_names}
    inverse_repdict = {}
    leftclus = {}
    rightclus = {}

    for i in range(len(mixed_clusters)):
        for rep in mixed_clusters[i]:
            repin_dict[rep[0]][rep[-1]] = i

    for rep in repin_names:
        newkey = (repin_dict[rep]['L'], repin_dict[rep]['R'])
        if newkey not in inverse_repdict.keys():
            inverse_repdict[newkey] = []
        inverse_repdict[newkey].append(rep)
        if repin_dict[rep]['L'] not in leftclus.keys():
            leftclus[repin_dict[rep]['L']] = []
        leftclus[repin_dict[rep]['L']].append(repin_dict[rep]['R'])
        if repin_dict[rep]['R'] not in rightclus.keys():
            rightclus[repin_dict[rep]['R']] = []
        rightclus[repin_dict[rep]['R']].append(repin_dict[rep]['L'])
    for key in leftclus.keys():
        leftclus[key] = list(set(leftclus[key]))
    for key in rightclus.keys():
        rightclus[key] = list(set(rightclus[key]))

    clusters = deepcopy(inverse_repdict)

    for key, val in inverse_repdict.items():
        if len(set(val)) > 1:
            combval = '\n'.join(val)
            logstring = f">{key[0]}_{key[1]} self\n{combval}"
            logging_mergers.append(logstring)

    logging_repin_dict = {}
    for key, val in inverse_repdict.items():
        for entry in val:
            logging_repin_dict[entry] = key


def nearby_repins(gen, posa, posb):
    near_reps = []
    for rep in repins_per_genome[gen]:
        x = int(rep.split(" ")[1])
        y = int(rep.split(" ")[2])
        diff_left = min(abs(posa - x), abs(posa - y))
        diff_right = min(abs(posb - x), abs(posb - y))
        # ------------------------------
        # DISCLAIMER
        # [rep, 'left'] => ---repin-sequence--------posa--flanking-gene--posb-----
        # Here [rep, 'left'] means that the repin is to the left of the position of interest
        # [rep, 'right'] => -----posa--flanking-gene--posb--------repin-sequence---
        # ------------------------------
        if diff_left <= closest_repin_criteria:
            near_reps.append([rep, 'left'])
        if diff_right <= closest_repin_criteria:
            near_reps.append([rep, 'right'])
    return near_reps


def setup_flank_matches():
    global left_flank, right_flank, unique_id_cards
    mixed_clusters = []
    switch_dir = {'left': 'R', 'right': 'L'}
    for key in repins_with_1k_flanks.keys():
        progress_bar(int(len(mixed_clusters) / 2) + 1, allreplength * 1.1)
        repin = repins_with_1k_flanks[key]
        lhs = repin[2]
        lhs_hits = prepare_datasets.search_blastdb(
            all_parameters['bank'], lhs, flank_gene_param)
        mixed_clusters.append([])
        for hit in lhs_hits:
            # hit[1] = hit[1][:3].lower() + hit[1][3:]
            near_reps = nearby_repins(hit[1], hit[4], hit[5])
            for rep in near_reps:
                # ------------------------------
                # DISCLAIMER
                # Here, rep_left => the flanking region is the left flank of this repin
                # Opposite of what is said in near_repins
                # ------------------------------
                mixed_clusters[-1].append([rep[0], switch_dir[rep[1]]])

        rhs = repin[4]
        rhs_hits = prepare_datasets.search_blastdb(
            all_parameters['bank'], rhs, flank_gene_param)
        mixed_clusters.append([])
        for hit in rhs_hits:
            # hit[1] = hit[1][:3].lower() + hit[1][3:]
            near_reps = nearby_repins(hit[1], hit[4], hit[5])
            for rep in near_reps:
                # ------------------------------
                # DISCLAIMER
                # Here, rep_left => the flanking region is the left flank of this repin
                # Opposite of what is said in near_repins
                # ------------------------------
                mixed_clusters[-1].append([rep[0], switch_dir[rep[1]]])

    pickle.dump(mixed_clusters, open(
        temp_files + f"mixed_clusters_{todaysdate}.p", 'wb'))
    prepare_repin_dict(mixed_clusters)


def setup(bank_path):
    global repins_with_1k_flanks, clusters, repin_names, repins_per_genome
    global all_parameters, genomes_list, flank_gene_param
    global allreplength, temp_files, genome_blastdb

    all_parameters = pickle.load(open(f"{bank_path}/all_parameters.p", "rb"))

    temp_files = all_parameters['bank'] + temp_files
    genome_blastdb = all_parameters['bank'] + genome_blastdb

    genomes_list = all_parameters['genomes']
    flank_gene_param = {
        'pident': all_parameters['pident'], 'lengthmatch': all_parameters['coverage']}

    prepare_datasets.setup_blastdb(all_parameters['bank'])
    rw1k_path = f"{all_parameters['bank']}/repin_with_1k_flanks.p"
    repins_with_1k_flanks = pickle.load(open(rw1k_path, "rb"))
    clusters = [key for key in repins_with_1k_flanks.keys()]
    repin_names = [key for key in repins_with_1k_flanks.keys()]
    repins_per_genome = {
        gen.split("/")[-1].split(".")[0]: [] for gen in genomes_list}
    for key in clusters:
        gen = key.split(" ")[0]
        repins_per_genome[gen].append(key)

    allreplength = len(list(repins_with_1k_flanks.keys()))
    if not fast_mode_testin_only:
        setup_flank_matches()
    else:
        mixclus_file = f"mixed_clusters_{todaysdate}.p"
        mixed_clusters = pickle.load(open(mixclus_file, "rb"))
        prepare_repin_dict(mixed_clusters)


def flankclusterer():
    global clusters, logging_mergers
    progress_bar(allreplength * 1.05, allreplength * 1.1)
    # ---------------------Block Starts---------------------
    # Finding Middle Gene Deletion Cases
    flank_gene_pair = []
    for rep in repin_names:
        a = repin_dict[rep]['L']
        b = repin_dict[rep]['R']
        if a != b:
            flank_gene_pair.append(tuple(set([a, b])))
        else:
            flank_gene_pair.append((a, b))
    flank_gene_graph = nx.Graph()
    flank_gene_graph.add_edges_from(flank_gene_pair)
    all_cliques = nx.enumerate_all_cliques(flank_gene_graph)
    triad_cliques = [x for x in all_cliques if len(x) == 3 and -1 not in x]
    triad_cliques = sum([list(itercomb(x, 2)) for x in triad_cliques], [])
    # ---------------------Block Ends---------------------

    to_merge = []
    for key1 in clusters.keys():
        for key2 in clusters.keys():
            if key1 == key2:
                continue
            if key1 in triad_cliques or key2 in triad_cliques:
                logging_trans.append(f"{key1}_{key2}:Triad")
                continue
            if (key1[1], key1[0]) in triad_cliques or (key2[1], key2[0]) in triad_cliques:
                logging_trans.append(f"{key1}_{key2}:Triad")
                continue
            if key1[0] == key2[0] and key1[0] != -1:
                if len(rightclus[key1[1]]) > 1:
                    logging_trans.append(
                        f"{key1}_{key2}:TS:{rightclus[key1[1]]}")
                    continue
                if len(rightclus[key2[1]]) > 1:
                    logging_trans.append(
                        f"{key1}_{key2}:TS:{rightclus[key2[1]]}")
                    continue
                combval = "\n".join(clusters[key1] + clusters[key2])
                logstring = f">{key1[0]}_{key1[1]} with {key2[0]}_{key2[1]}\n{combval}"
                logging_mergers.append(logstring)
                to_merge.append((key1, key2))
            if key1[1] == key2[1] and key1[1] != -1:
                if len(leftclus[key1[0]]) > 1:
                    logging_trans.append(
                        f"{key1}_{key2}:TS:{leftclus[key1[0]]}")

                    continue
                if len(leftclus[key2[0]]) > 1:
                    logging_trans.append(
                        f"{key1}_{key2}:TS:{leftclus[key2[0]]}")
                    continue
                combval = "\n".join(clusters[key1] + clusters[key2])
                logstring = f">{key1[0]}_{key1[1]} with {key2[0]}_{key2[1]}\n{combval}"
                logging_mergers.append(logstring)
                to_merge.append((key1, key2))

    merge_graph = nx.Graph()
    merge_graph.add_edges_from(to_merge)
    merge_graph = [list(x) for x in list(nx.connected_components(merge_graph))]

    to_merge = merge_graph
    in_merge_queue = [j for i in to_merge for j in i]
    new_clusters = []
    for key in clusters.keys():
        if key not in in_merge_queue:
            new_clusters.append(clusters[key])
    for group in to_merge:
        mergedclus = []
        for key in group:
            mergedclus += clusters[key]
        new_clusters.append(mergedclus)

    clusters = deepcopy(new_clusters)


def print_out_clusters():
    progress_bar(allreplength * 1.1, allreplength * 1.1)

    outfile = open(f"{all_parameters['out']}/clusters_{todaysdate}.txt", "w")
    meta_output = open(
        f"{all_parameters['out']}/meta_cluster_{todaysdate}.txt", "w")
    for i in range(len(clusters)):
        for rep in clusters[i]:
            a = repins_with_1k_flanks[rep][0]
            b = repins_with_1k_flanks[rep][1]
            c = repins_with_1k_flanks[rep][3]
            outfile.write(f"{i} {a} {b} {c}\n")
            kiss = logging_repin_dict[a]
            meta_output.write(
                f"{i} {a} leftflank_{kiss[0]} rightflank_{kiss[1]}\n")
    outfile.close()
    meta_output.close()

    with open(f"{all_parameters['out']}/path_making_{todaysdate}.txt", "w") as f:
        f.write("\n".join(logging_mergers))

    with open(f"{all_parameters['out']}/missed_hits_{todaysdate}.txt", "w") as f:
        f.write("\n".join(logging_trans))


def main(bank_path):
    setup(bank_path)
    flankclusterer()
    print_out_clusters()


if __name__ == "__main__":
    st = time.time()
    input_params = {
        "repin": "./input/repins.txt",
        "genomes": "./input/genomes",
        "out": "./output/",
        "win": 250,
        "fsize": 1000,
        "pident": 90,
        "coverage": 90
    }
    pickle.dump(input_params, open("./bank/all_parameters.p", "wb"))
    main(".")
    print("Runtime: {:.2}s".format(time.time() - st))
