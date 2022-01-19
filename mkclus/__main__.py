import click
import os
from Bio import SeqIO
import pickle
file_path = os.path.realpath(__file__)
file_path = os.path.relpath(file_path, os.getcwd())
file_path = "./" + file_path.split("/")[0]
os.chdir(file_path)
import head_of_clustering


def quick_check_files(repin, genomes):
    if not os.path.isfile(repin):
        print("File containing REPINs does not exist")
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
            except Exception:
                print(f"Error in genome fasta files in {genomes}/{gen}")
                print("Exiting......")
                exit()


@click.command()
@click.option('--repin', prompt="Repin File", help='Path to file containing repin sequences')
@click.option('--genomes', prompt="Genomes Directory", help='Path to directory containing genomes')
@click.option('--out', help="Output file destination", default='./cluster_output/')
@click.option('--win', help="Repin flanking window", default=250)
@click.option('--fsize', help="Size of flanking region", default=1000)
@click.option('--pident', help="Percentage sequence similarity", default=90)
@click.option('--coverage', help="Minimum length of alignment", default=90)
def main(repin, genomes, out, win, fsize, pident, coverage):
    all_parameters = {
        "repin": os.path.abspath(repin),
        "genomes": os.path.abspath(genomes),
        "out": os.path.abspath(out),
        "win": win,
        "fsize": fsize,
        "pident": pident,
        "coverage": coverage
    }
    quick_check_files(all_parameters['repin'], all_parameters['genomes'])
    pickle.dump(all_parameters, open(
        "./bank/all_parameters.p", "wb"))
    head_of_clustering.main()


if __name__ == '__main__':
    main()
