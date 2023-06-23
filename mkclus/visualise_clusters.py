import os
import click
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.widgets import Slider
from math import comb


ISOLATION_LIMIT = 0.9
CACHE = True


def make_network(cnum):
    node_color1, node_color2 = [], []

    G1 = nx.DiGraph(just_pairs['left'][cnum]).to_undirected()
    for pair in G1.edges:
        G1.edges[pair]["weight"] = matrix_tags['left'][pair[0]][pair[1]]
    if cnum in isolates['left']:
        for node in isolates['left'][cnum]:
            G1.add_node(node)
    for node in G1:
        if cnum in isolates['left']:
            if node in isolates['left'][cnum]:
                node_color1.append('red')
            else:
                node_color1.append('blue')
        else:
            node_color1.append('blue')
    nx.relabel_nodes(G1, unique_ids, copy=False)

    G2 = nx.DiGraph(just_pairs['right'][cnum]).to_undirected()
    for pair in G2.edges:
        G2.edges[pair]["weight"] = matrix_tags['right'][pair[0]][pair[1]]
    if cnum in isolates['right']:
        for node in isolates['right'][cnum]:
            G2.add_node(node)
    for node in G2:
        if cnum in isolates['right']:
            if node in isolates['right'][cnum]:
                node_color2.append('red')
            else:
                node_color2.append('blue')
        else:
            node_color2.append('blue')
    nx.relabel_nodes(G2, unique_ids, copy=False)

    try:
        edges1, weights1 = zip(*nx.get_edge_attributes(G1, 'weight').items())
        edges2, weights2 = zip(*nx.get_edge_attributes(G2, 'weight').items())
    except ValueError:
        edges1, weights1 = None, None
        edges2, weights2 = None, None

    xo1 = (G1, edges1, weights1, node_color1)
    xo2 = (G2, edges2, weights2, node_color2)

    return (xo1, xo2)


def process():
    global infile_name, outdir, infile, clusters, matrix_tags, just_pairs, unique_ids, isolates, cliques, connects

    # if os.path.isdir(outdir) and CACHE:
    #     clusters = pickle.load(open("./bankclusters", "rb"))
    #     matrix_tags = pickle.load(open("./bankmatrix_tags", "rb"))
    #     just_pairs = pickle.load(open("./bankjust_pairs", "rb"))
    #     unique_ids = pickle.load(open("./bankunique_ids", "rb"))
    #     isolates = pickle.load(open("./bankisolates", "rb"))
    #     cliques = pickle.load(open("./bankcliques", "rb"))
    #     connects = pickle.load(open("./bankconnects", "rb"))
    #     return 0

    with open(infile_name, "r") as f:
        infile = f.read().split("\n")
        infile = [x.split() for x in infile]

    unique_ids = {}
    id_by_gens = {}
    isolates = {'left': {}, 'right': {}}
    clusters = {}
    matrix_tags = {'left': {}, 'right': {}}
    just_pairs = {'left': {}, 'right': {}}
    left_cliques, right_cliques = [], []
    left_connects, right_connects = [], []
    for num, line in enumerate(infile):
        if len(line) < 1:
            continue
        try:
            a = int(line[0])
            memory = " ".join(line[1:])
            if a not in clusters:
                clusters[a] = []
            clusters[a].append(memory)
            gen = memory.split(" ")[0]
            p1 = int(memory.split(" ")[1])
            if gen not in id_by_gens.keys():
                id_by_gens[gen] = {}
            id_by_gens[gen][memory] = p1
        except Exception:
            if memory not in matrix_tags[line[0]]:
                matrix_tags[line[0]][memory] = {}
            matrix_tags[line[0]][memory] = [float(x) for x in line[1:]]

    for num, reps in clusters.items():
        for key in just_pairs:
            just_pairs[key][num] = []
        for rep in reps:
            for key, val in matrix_tags.items():
                matrix_tags[key][rep] = dict(zip(reps, val[rep]))
            for key, val in matrix_tags.items():
                for k2, v2 in matrix_tags[key][rep].items():
                    if (rep, k2) not in just_pairs[key][num] and rep != k2:
                        if matrix_tags[key][rep][k2] < ISOLATION_LIMIT:
                            if num not in isolates[key]:
                                isolates[key][num] = []
                            isolates[key][num].extend([rep, k2])
                            continue
                        just_pairs[key][num].append((rep, k2))

    for key, val in id_by_gens.items():
        v2 = dict(sorted(val.items(), key=lambda item: item[1]))
        ct = 1
        for k2, v2_1 in v2.items():
            unique_ids[k2] = k2.split(" ")[0] + f"_{ct}"
            ct += 1

    for x in range(len(clusters)):
        xo1, xo2 = make_network(x)

        if xo1[1] is None:
            left_cliques.append(1)
            left_connects.append(1)
            right_cliques.append(1)
            right_connects.append(1)
            continue

        left_cliques.append(len(list(nx.find_cliques(xo1[0]))))
        right_cliques.append(len(list(nx.find_cliques(xo2[0]))))

        left_connects.append(xo1[0].number_of_edges() /
                             comb(xo1[0].number_of_nodes(), 2))
        right_connects.append(
            xo2[0].number_of_edges() / comb(xo2[0].number_of_nodes(), 2))

    cliques = left_cliques, right_cliques
    connects = left_connects, right_connects

    # os.system(f"mkdir {outdir}")
    # pickle.dump(cliques, open(outdir + "cliques", "wb"))
    # pickle.dump(connects, open(outdir + "connects", "wb"))
    # pickle.dump(clusters, open(outdir + "clusters", "wb"))
    # pickle.dump(matrix_tags, open(outdir + "matrix_tags", "wb"))
    # pickle.dump(just_pairs, open(outdir + "just_pairs", "wb"))
    # pickle.dump(unique_ids, open(outdir + "unique_ids", "wb"))
    # pickle.dump(isolates, open(outdir + "isolates", "wb"))


def draw_connecting_graphs(cnum):
    xo1, xo2 = make_network(cnum)

    if xo1[1] is None:
        print("Cluster has only one REPIN. Nothing to draw")
        return 400

    G1, edges1, weights1, node_color1 = xo1
    G2, edges2, weights2, node_color2 = xo2

    gs = gridspec.GridSpec(2, 2, height_ratios=[30, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])

    plt.axes(ax0)
    plt.title("Left Flanking Region")
    # nx.draw_circular(G1, edgelist=edges1, edge_color=weights1, with_labels=True, edge_cmap=plt.cm.cool)
    pos1 = nx.spring_layout(G1)
    nx.draw_networkx(G1, pos1, edgelist=edges1, edge_color=weights1,
                     with_labels=True, edge_cmap=plt.cm.cool)
    # labels = nx.get_edge_attributes(G1,'weight')
    # nx.draw_networkx_edge_labels(G1,pos1,edge_labels=labels)
    plt.axes(ax2)
    w1 = list(weights1)
    w1.sort()
    plt.title("{} to {}".format(min(w1), max(w1)))
    plt.imshow([w1], cmap="cool", aspect="auto")
    plt.axis('off')

    plt.axes(ax1)
    plt.title("Right Flanking Region")
    # nx.draw_circular(G2, edgelist=edges2, edge_color=weights2, with_labels=True, edge_cmap=plt.cm.cool)
    pos2 = nx.spring_layout(G2)
    nx.draw_networkx(G2, pos2, edgelist=edges2, edge_color=weights2,
                     with_labels=True, edge_cmap=plt.cm.cool)
    # labels = nx.get_edge_attributes(G2,'weight')
    # nx.draw_networkx_edge_labels(G2,pos2,edge_labels=labels)
    plt.axes(ax3)
    w2 = list(weights2)
    w2.sort()
    plt.title("{} to {}".format(min(w2), max(w2)))
    plt.imshow([w2], cmap="cool", aspect="auto")
    plt.axis('off')

    plt.suptitle(f"Sequence Similarity >= {ISOLATION_LIMIT}")
    plt.tight_layout()
    plt.draw()
    plt.show()


def plot_clique_plots(summary_types):
    # left = [random.randint(2,5) for x in range(1100)]
    # right = [random.randint(1,2) for x in range(1100)]

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)

    if summary_types == 1:
        plt.ylabel("Number of Cliques")
        plt.title("Number of Cliques in the Flanking Sequences for Each Cluster")
        left, right = cliques
    elif summary_types == 2:
        plt.ylabel("Percentage Connections")
        plt.title(
            "Percentage of Total Possible Connections found amongst the Flanking Sequences for Each Cluster")
        left, right = connects

    length_list = [len(clusters[x]) for x in clusters]

    xaxis = np.arange(len(left))
    plt.plot(xaxis, left, label="Left")
    plt.plot(xaxis, right, label="Right")
    plt.plot(xaxis, length_list, label="Size")
    plt.xticks(xaxis, rotation=45)
    plt.legend()
    BUFFER = 20

    # ymin = (min(left + right) * 0.85)
    ymin = -10
    ymax = max(left + right + length_list) * 1.5
    plt.axis([0, BUFFER, ymin, ymax])

    axpos = plt.axes([0.2, 0.1, 0.65, 0.03])

    spos = Slider(axpos, 'Cluster', 0, max(xaxis) - (BUFFER * 0.8))

    def update(val):
        pos = spos.val
        ax.axis([pos, pos + BUFFER, ymin, ymax])
        fig.canvas.draw_idle()

    spos.on_changed(update)

    plt.show()


def summary_plots(summary_types, outpath, toprint):
    plt.figure(figsize=(8, 6))
    if summary_types == 2:
        left, right = connects
        plt.suptitle(
            "Distribution of the Connectivity amongst the Flanking Sequences of a Cluster")
        x_lab = "Percentage Connectivity within a cluster"
    elif summary_types == 1:
        left, right = cliques
        plt.suptitle(
            "Distribution of the Clique Sizes amongst the Flanking Sequences of a Cluster")
        x_lab = "Number of Cliques withn a cluster"
    else:
        exit("Incorrect value for argument --summary_types. Exiting...")

    if summary_types == 2:
        left = np.array(left) * 100
        right = np.array(right) * 100
    unique_left = len(list(set(left)))
    unique_right = len(list(set(right)))
    leftbins = np.arange(unique_left + 1) - 0.5
    rightbins = np.arange(unique_right + 1) - 0.5
    percent_bins = np.arange(0, 110, 10)

    plt.subplot(1, 2, 1)
    if summary_types == 2:
        plt.hist(left, percent_bins)
        plt.xticks(percent_bins)
    else:
        plt.hist(left, leftbins, density=True)
        plt.xticks(range(unique_left))
    plt.title("Left Flanking Sequnce")
    plt.xlabel(x_lab)
    plt.ylabel("Number of occurrences")

    plt.subplot(1, 2, 2)
    if summary_types == 2:
        plt.hist(right, percent_bins)
        plt.xticks(percent_bins)
    else:
        plt.hist(right, rightbins, density=True)
        plt.xticks(range(unique_right))
    plt.title("Right Flanking Sequnce")
    plt.xlabel(x_lab)
    plt.ylabel("Number of occurrences")

    plt.tight_layout()
    if toprint:
        plt.savefig(
            f"{outpath}/summary_histogram_{summary_types}.png", dpi=200)
    else:
        plt.show()


def mkclus_main(metafilename, outpath, cluster=-1, summary=0):
    global infile_name, outdir, CACHE
    infile_name = metafilename
    outdir = "./bank/"
    # To remove any cache
    CACHE = 1
    process()
    if cluster >= 0:
        draw_connecting_graphs(cluster)
    elif summary == 0:
        # For cliques based output
        summary_plots(summary_types=1, outpath=outpath, toprint=True)
        # For percentage connectivity based output
        summary_plots(summary_types=2, outpath=outpath, toprint=True)
    else:
        # For cliques based output
        plot_clique_plots(summary_types=1)


@click.command()
@click.option('--file', prompt='Please Input REPIN Cluster Meta File')
@click.option('--cluster', default=-1, help="Choose cluster number to display connected graphs")
@click.option('--summary', default=0, help="Choose '0' for Histogram summary or '1' for detailed line plot summary")
@click.option('--summary_types', default=1, help="Choose '1' for Clique based and '2' for Percentage Connections based summary")
@click.option('--recache', default=0, help="If set to 1, removes cache, and reruns the program")
def main(file, cluster, summary, summary_types, recache):
    global infile_name, outdir, CACHE
    if not os.path.isfile(file):
        exit("File not provided or does not exist. Exited.")
    infile_name = file
    outdir = f"./dat_{infile_name.split('.')[0]}/"

    CACHE = not recache

    process()
    if cluster >= 0:
        draw_connecting_graphs(cluster)
    elif summary == 0:
        # Extra parameters because of how we run this
        summary_plots(summary_types, outpath=None, toprint=False)
    else:
        plot_clique_plots(summary_types)


if __name__ == "__main__":
    main()
