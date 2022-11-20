import click
import os
from Bio import SeqIO
import pickle
from mkclus import head_of_clustering

all_parameters = {}


def get_files_from_rarefan(rarefan_path):
    genomes = []
    for file in list(os.walk(rarefan_path))[0][1]:
        if "_" in file and "_0" not in file: 
            genomes.append(file)

    allrepins = []
    for gen in genomes:
        try:
            repins = open(rarefan_path + "/" + gen + "/" + gen + ".ss").read()
        except Exception:
            continue
        repins = [i.replace("\t", "\n").split("\n")
                  for i in repins.split(">") if len(i) > 0]
        repins = [i[1:] for i in repins]
        repins = [[j for j in i if len(j) > 0] for i in repins]

        repintype = gen.split("_")[-1]

        for repin in repins:
            rseq = repin[-1]
            rname = repin[:-1]
            for rep in rname:
                rep = repin[0].split("_")
                newr = "{} {} {} type{} {}".format(
                    gen.split("_")[0], rep[1], rep[2], repintype, rseq)
                allrepins.append(newr)

    allrepins = "\n".join(allrepins)
    open(rarefan_path + "/sortedrepins.txt", "w").write(allrepins)

def quick_check_files(repin, genomes):
    if not os.path.isfile(repin):
        print("File containing REPINs does not exist")
        print("If you have RAREFAN output use the tag --withrarefan 1")
        print("Exiting......")
        exit()
    if not os.path.isdir(genomes):
        print("Genome directory does not exist")
        print("Exiting......")
        exit()
    else:
        gens = next(os.walk(genomes), (None, None, []))[2]
        for gen in gens:
            try:
                list(SeqIO.parse(genomes + "/" + gen, 'fasta'))[0]
                all_parameters["genomes"].append(genomes + "/" + gen)
            except Exception:
                print(f"Ignoring {genomes}/{gen} - Not a fasta file")


@click.command()
@click.option('--repin', prompt="Repin File or RAREFAN Dir", help='Path to file containing repin sequences or RAREFAN Output')
@click.option('--withrarefan', help='Using RAREFAN output? 1 for yes, 0 for no', default=0)
@click.option('--genomes', prompt="Genomes Directory", help='Path to directory containing genomes')
@click.option('--out', help="Output file destination", default='./cluster_output/')
@click.option('--win', help="Repin flanking window", default=250)
@click.option('--fsize', help="Size of flanking region", default=1000)
@click.option('--pident', help="Percentage sequence similarity", default=90)
@click.option('--coverage', help="Minimum length of alignment", default=90)
def main(repin, withrarefan, genomes, out, win, fsize, pident, coverage):
    global all_parameters

    all_parameters = {
        "repin": os.path.abspath(repin),
        "genomes": [],
        "out": os.path.abspath(out),
        "win": win,
        "fsize": fsize,
        "pident": pident,
        "coverage": coverage
    }

    if withrarefan == 1:
        get_files_from_rarefan(all_parameters['repin'])
        all_parameters['repin'] += "/sortedrepins.txt"

    quick_check_files(all_parameters['repin'], os.path.abspath(genomes))
    os.system("mkdir {}".format("./bank/"))
    os.system("mkdir {}".format("./bank/dumpyard"))
    os.system("mkdir {}".format("./bank/genomes_blastdb"))

    pickle.dump(all_parameters, open(
        "./bank/all_parameters.p", "wb"))

    head_of_clustering.main()

    os.system("rm -rf {}".format("./bank/"))


if __name__ == '__main__':
    try:
        main()
        print("Program Completed.")
    except Exception as e:
        os.system("rm -rf ./bank")
        exit(f"repinclusterer encountered an error:\n{e}\nExiting...")
