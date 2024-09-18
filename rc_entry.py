import click
import os
from Bio import SeqIO
import pickle
from mkclus import head_of_clustering
import time
import random
import re

todaysdate = time.strftime("%b%d") + "_" + str(random.random())[3:6]
all_parameters = {}


def get_files_from_rarefan(rarefan_path, reptypes):

    if os.path.isfile(rarefan_path):
        return rarefan_path

    genomes = []
    for file in list(os.walk(rarefan_path))[0][1]:
        matches = re.finditer(r"(.*)_\d", file, re.MULTILINE)
        matches = [match.group() for match in matches]
        if len(matches) > 0:
            genomes.append(matches[0])

    genomes.sort()
    allrepins = []
    remove_repeats = {}
    for gen in genomes:
        try:
            repins = open(rarefan_path + "/" + gen + "/" + gen + ".ss").read()
        except Exception:
            continue

        repins = [i.replace("\t", "\n").split("\n")
                  for i in repins.split(">") if len(i) > 0]
        repins = [i[1:] for i in repins]
        repins = [[j for j in i if len(j) > 0] for i in repins]

        genname = "_".join(gen.split("_")[:-1])
        repintype = int(gen.split("_")[-1])

        if reptypes is not None:
            if repintype not in reptypes:
                continue

        if genname not in remove_repeats.keys():
            remove_repeats[genname] = {}
        remove_repeats[genname][repintype] = []

        for repin in repins:
            rseq = repin[-1]
            rname = repin[:-1]
            for rep in rname:
                rep = rep.split("_")
                newr = "{} {} {} type{} {}".format(
                    genname, rep[1], rep[2], repintype, rseq)
                
                keep = True
                
                # Check size of sequence to determine if it's a REP or REPIN
                if abs(int(rep[1]) - int(rep[2])) <= all_parameters['MINREPINSIZE']:
                    keep = False
                
                for rtype, val in remove_repeats[genname].items():
                    if f"{rep[1]}_{rep[2]}" in val:
                        keep = False
                remove_repeats[genname][repintype].append(f"{rep[1]}_{rep[2]}")
                if keep:
                    allrepins.append(newr)

    allrepins = "\n".join(allrepins)
    final_filename = rarefan_path + "/sortedrepins.txt"
    with open(final_filename, "w") as f:
        f.write(allrepins)
    return final_filename


def quick_check_files(repin, genomes):
    if not os.path.isfile(repin):
        exit("File containing REPINs does not exist\nExiting......")

    with open(repin, "r") as f:
        existing_in_gens = f.read().split("\n")
        existing_in_gens = [x.split()[0]
                            for x in existing_in_gens if len(x.split()) > 0]

    store_extensions = {}

    if not os.path.isdir(genomes):
        print("Genome directory does not exist")
        exit("Exiting......")
    else:
        gens = next(os.walk(genomes), (None, None, []))[2]
        for gen in gens:
            try:
                list(SeqIO.parse(genomes + "/" + gen, 'fasta'))[0]
            except Exception:
                if ".DS_Store" not in gen:
                    print(f"Ignoring {genomes}/{gen} - Not a fasta file")
                continue
            filesplit = gen.split(".")
            all_parameters["genomes"].append(filesplit[0])
            store_extensions[filesplit[0]] = ".".join(filesplit[1:])

        for gen in existing_in_gens:
            if gen not in all_parameters["genomes"]:
                exit(
                    f"Genome fasta file for {gen} not provided but REPINs from {gen} exist\nExisting Gens: {','.join(existing_in_gens)}")

        extraas = []
        for gen in all_parameters["genomes"]:
            if gen not in existing_in_gens:
                extraas.append(gen)

        all_parameters["genomes"] = [
            x for x in all_parameters["genomes"] if x not in extraas]

        for pos, val in enumerate(all_parameters["genomes"]):
            if len(store_extensions[val]) < 1:
                all_parameters["genomes"][pos] = f"{genomes}/{val}"
            else:
                all_parameters["genomes"][pos] = f"{genomes}/{val}.{store_extensions[val]}"

            store_genomes = list(SeqIO.parse(
                all_parameters["genomes"][pos], "fasta"))
            if len(store_genomes) > 1:
                warning_message = f"""The fasta file {gen} contains more than one sequence - Genome files should have a single fasta entry. Proceeding by selecting only the first entry in the fasta file.
                Note: The program is not designed to work with contigs.
                Genome file is located at {genomes}/{gen}"""
                newpath = all_parameters["bank"] + f"/{gen}_firstentry"
                SeqIO.write(store_genomes[0], newpath, "fasta")
                print(warning_message)



@click.command()
@click.option('--repin', help='Path to file containing repin sequences or RAREFAN Output', default=None)
@click.option('--genomes', help='Path to directory containing genomes', default=None)
@click.option('--visualfile', help='Metafile produced by this program', default=None)
@click.option('--visualtype', help='Cluster number to visualise cluster (or) -1 for Detailed Summary for Clique based (or) -1 for Detailed Summary for Connectivity based', default=None)
@click.option('--out', help="Output file destination", default='./cluster_output')
@click.option('--win', help="Repin flanking window", default=250)
@click.option('--fsize', help="Size of flanking region", default=1000)
@click.option('--pident', help="Percentage sequence similarity", default=90)
@click.option('--coverage', help="Minimum length of alignment", default=90)
@click.option('--reptypes', help="Mention the specific repin types to accept from rarefan output")
@click.option('--minrepinsize', help="To be used if REP singlets and doublets are to be considiered", default=50)
def main(repin, genomes, visualfile, visualtype, out, win, fsize, pident, coverage, reptypes, minrepinsize):
    global all_parameters

    if repin is None and genomes is None and visualfile is None and visualtype is None:
        exit("Appropriate parameters ")

    if reptypes is not None:
        reptypes = [int(x) for x in reptypes.split(",")]

    all_parameters = {
        "repin": os.path.abspath(repin),
        "genomes": [],
        "out": os.path.abspath(out) + f"_{todaysdate}/",
        "win": win,
        "fsize": fsize,
        "pident": pident,
        "coverage": coverage,
        "reptypes": reptypes,
        "MINREPINSIZE": minrepinsize
    }

    # File names cannot contain whitespaces
    if " " in all_parameters['out'] or "\t" in all_parameters['out']:
        exit("Filename / Filepath cannot contain whitespaces. Aborting program...")

    all_parameters['bank'] = all_parameters['out'] + "bank/"

    all_parameters['repin'] = get_files_from_rarefan(
        all_parameters['repin'], all_parameters['reptypes'])

    # Check validity of all files and genome files
    quick_check_files(all_parameters['repin'], os.path.abspath(genomes))

    # Make Temporary files
    if not os.path.isdir(all_parameters['out']):
        os.system("mkdir {}".format(all_parameters['out']))
    os.system("mkdir {}".format(f"{all_parameters['bank']}/"))
    os.system("mkdir {}".format(f"{all_parameters['bank']}/dumpyard"))
    os.system("mkdir {}".format(f"{all_parameters['bank']}/genomes_blastdb"))
    pickle.dump(all_parameters, open(
        f"{all_parameters['bank']}/all_parameters.p", "wb"))

    # Begin the main clustering program
    head_of_clustering.main(all_parameters['bank'])

    # Clear temp files after running program
    os.system("rm -rf {}".format(f"{all_parameters['bank']}/"))


if __name__ == '__main__':
    try:
        main()
        print("Program Completed.")
    except Exception as e:
        os.system(f"rm -rf {all_parameters['out']}")
        exit(f"reporth encountered an error:\n{e}\nExiting...")
