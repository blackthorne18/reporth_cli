# -------------------------------------------------------------------------------------------
"""
Input Files Required:
1.  repin_with_1k_flanks (from prepare_datasets.py)
    Path = ./bank/repin_with_1k_flanks.p
2.  blastgenomedb (from prepare_datasets.py)
    Path = ./bank/genomes_blastdb/allgenomes

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
import prepare_datasets
import networkx as nx
from copy import deepcopy
from itertools import combinations as itercomb
# import numpy as np

todaysdate = time.strftime("%b%d")
genome_blastdb = "./bank/genomes_blastdb/allgenomes"
temp_files = "./bank/dumpyard/"

genomes_list = []
repins_with_1k_flanks = {}
clusters = []
repin_dict, leftclus, rightclus = {}, {}, {}
repin_names = []
move_tracker = []
repins_per_genome = {}
closest_repin_criteria = 500
fast_mode_testin_only = False
all_parameters = {}
flank_gene_param = {}


def prepare_repin_dict(mixed_clusters):
    global repin_dict, leftclus, rightclus, clusters
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
        print("Setting up match", int(len(mixed_clusters) / 2))
        repin = repins_with_1k_flanks[key]
        lhs = repin[2]
        lhs_hits = prepare_datasets.search_blastdb(lhs, flank_gene_param)
        mixed_clusters.append([])
        for hit in lhs_hits:
            hit[1] = hit[1][:3].lower() + hit[1][3:]
            near_reps = nearby_repins(hit[1], hit[4], hit[5])
            for rep in near_reps:
                # ------------------------------
                # DISCLAIMER
                # Here, rep_left => the flanking region is the left flank of this repin
                # Opposite of what is said in near_repins
                # ------------------------------
                mixed_clusters[-1].append([rep[0], switch_dir[rep[1]]])

        rhs = repin[4]
        rhs_hits = prepare_datasets.search_blastdb(rhs, flank_gene_param)
        mixed_clusters.append([])
        for hit in rhs_hits:
            hit[1] = hit[1][:3].lower() + hit[1][3:]
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


def setup():
    global repins_with_1k_flanks, clusters, repin_names, repins_per_genome
    global all_parameters, genomes_list, flank_gene_param

    all_parameters = pickle.load(open("./bank/all_parameters.p", "rb"))
    genomes_list = [gen[:-4]
                    for gen in next(os.walk(all_parameters['genomes']), (None, None, []))[2]]
    flank_gene_param = {
        'pident': all_parameters['pident'], 'lengthmatch': all_parameters['coverage']}

    prepare_datasets.setup_blastdb()
    rw1k_path = "./bank/repin_with_1k_flanks.p"
    repins_with_1k_flanks = pickle.load(open(rw1k_path, "rb"))
    clusters = [key for key in repins_with_1k_flanks.keys()]
    repin_names = [key for key in repins_with_1k_flanks.keys()]
    repins_per_genome = {gen: [] for gen in genomes_list}
    for key in clusters:
        gen = key.split(" ")[0]
        repins_per_genome[gen].append(key)

    if not fast_mode_testin_only:
        setup_flank_matches()
    else:
        mixclus_file = "./bank/dumpyard/mixed_clusters_Jan08.p"
        mixed_clusters = pickle.load(open(mixclus_file, "rb"))
        prepare_repin_dict(mixed_clusters)


def flankclusterer():
    global clusters
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
                continue
            if (key1[1], key1[0]) in triad_cliques or (key2[1], key2[0]) in triad_cliques:
                continue
            if key1[0] == key2[0] and key1[0] != -1:
                if len(rightclus[key2[1]]) > 1:
                    continue
                move_tracker.append([key1, key2])
                to_merge.append((key1, key2))
            if key1[1] == key2[1] and key1[1] != -1:
                if len(leftclus[key2[0]]) > 1:
                    continue
                move_tracker.append([key1, key2])
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
    outfile = open(all_parameters['out'], "w")
    for i in range(len(clusters)):
        for rep in clusters[i]:
            a = repins_with_1k_flanks[rep][0]
            b = repins_with_1k_flanks[rep][1]
            c = repins_with_1k_flanks[rep][3]
            outfile.write(f"{i} {a} {b} {c}\n")
    outfile.close()


def rmcleanup():
    os.system("rm -rf {}".format("./bank/"))
    os.system("mkdir {}".format("./bank/"))
    os.system("mkdir {}".format("./bank/dumpyard"))
    os.system("mkdir {}".format("./bank/genomes_blastdb"))


def main():
    setup()
    flankclusterer()
    print_out_clusters()
    rmcleanup()


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
    main()
    print("Runtime: {:.2}s".format(time.time() - st))
