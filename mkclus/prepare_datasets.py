# -------------------------------------------------------------------------------------------
"""
Input Files Required:
1.  List of genomes
    Path = ./input/genomes/*.fas
2.  Repins File
    Path = ./input/repins.p

Output Files:
1.  blastgenomedb
    Path = ./bank/genomes_blastdb/
2.  repin_with_1k_flanks
    Path = ./bank/repin_with_1k_flanks.p

setup_blastd()
    Creates the output files blastgenomedb and repin_with_1k_flanks

search_blastdb(sequence)
    Input = string - DNA Sequence
    Output = Top blast hits from all the genomes provided (one per genome)
    Blast Parameters
        Percentage Identity (pident) is the Levenshtein distance between the query and search sequence.
        length is the length of overlapping sequence
            flank_gene_param = {'pident': 90, 'lengthmatch': 90}
            This means that pident has to be over 90% and the length of the seqeunce that matches should
            be over 90% of the query
        We only return the best hit for each genome

"""
# -------------------------------------------------------------------------------------------
import os
import time
from Bio import SeqIO
import pickle
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

blast_path = "./bank/genomes_blastdb/"
temp_files = "./bank/dumpyard/"
all_parameters_path = "./bank/all_parameters.p"


def setup_blastdb():
    all_parameters = pickle.load(open(all_parameters_path, "rb"))
    genomes_path = all_parameters['genomes']
    genomes_list = next(
        os.walk(genomes_path), (None, None, []))[2]
    repins = open(all_parameters['repin'], "r").read().split("\n")

    flank_gene_range = {
        'window': all_parameters['win'], 'flanklength': all_parameters['fsize']}
    fgr = flank_gene_range['window'] + flank_gene_range['flanklength']

    genome_sequences = []
    repin_with_1k_flanks = {}
    repin_per_genome = {genome[:-4]: [] for genome in genomes_list}

    for rep in repins:
        gen = rep.split(" ")[0]
        posa = int(rep.split(" ")[1])
        posb = int(rep.split(" ")[2])
        col = rep.split(" ")[3]
        seq = rep.split(" ")[4]
        rep = [gen, posa, posb, col, seq]
        repin_per_genome[gen].append(rep)

    for genome in genomes_list:
        sequence = str(SeqIO.read(genomes_path + "/" + genome, "fasta").seq)
        for rep in repin_per_genome[genome[:-4]]:
            left_flank = sequence[rep[1] -
                                  fgr: rep[1] - flank_gene_range['window']]
            right_flank = sequence[rep[2] +
                                   flank_gene_range['window']: rep[2] + fgr]
            repname = rep[0] + " " + str(rep[1]) + " " + str(rep[2])
            repin_with_1k_flanks[repname] = [
                repname, rep[3], left_flank, rep[4], right_flank]

        genome_sequences.append(">{}\n{}".format(
            genome.split(".")[0], sequence))

    pickle.dump(repin_with_1k_flanks, open(
        "./bank/repin_with_1k_flanks.p", "wb"))

    genome_sequences = "\n".join(genome_sequences)
    open(blast_path + "allgenomes.fas", "w").write(genome_sequences)
    cmd = f"makeblastdb -in {blast_path}allgenomes.fas -out {blast_path}allgenomes -parse_seqids -dbtype nucl"

    cline = NcbimakeblastdbCommandline(
        dbtype="nucl", input_file=blast_path + "allgenomes.fas", out=blast_path + "allgenomes", parse_seqids=True)

    # Using BLAST CLI
    # os.system(cmd)
    # Using Biopython
    cline()


def search_blastdb(sequence, flank_gene_param):
    infile = temp_files + "test1_in.fas"
    outfile = temp_files + "test1_out.fas"
    open(infile, "w").write(">query_seq\n{}".format(sequence))

    cmd = f"blastn -query {infile} -db {blast_path+'allgenomes'} -out {outfile}"
    cmd_format = "-outfmt '6 qseqid sseqid pident length sstart send'"

    cline = NcbiblastnCommandline(query=infile, db=blast_path +
                                  "allgenomes", out=outfile, outfmt='6 qseqid sseqid pident length sstart send')

    # Using BLAST CLI
    # os.system(cmd + " " + cmd_format)
    # Using Biopython
    cline()

    outfile = [i.split("\t") for i in open(
        outfile, "r").read().split("\n") if len(i) > 0]
    good_output = []
    for i in range(len(outfile)):
        outfile[i] = [outfile[i][0], outfile[i][1]] + \
            [int(float(x)) for x in outfile[i][2:]]
        lengthmatch = int(
            100 * (abs(outfile[i][5] - outfile[i][4]) / len(sequence)))
        pident = outfile[i][2]
        if lengthmatch >= flank_gene_param['lengthmatch'] and pident >= flank_gene_param['pident']:
            good_output.append(outfile[i])

    # Makes sure that only one hit is recorded per genome and this hit is the highest hit

    seen_gens = {}
    to_keep = []
    for i in range(len(good_output)):
        gen = good_output[i][1]
        if gen not in seen_gens.keys():
            seen_gens[gen] = [0, 0]
        if good_output[i][2] > seen_gens[gen][0] and good_output[i][3] > seen_gens[gen][1]:
            to_keep.append(i)
            seen_gens[gen] = [good_output[i][2], good_output[i][3]]

    good_output = [good_output[i]
                   for i in range(len(good_output)) if i in to_keep]

    return good_output


def main():
    # setup_blastdb()
    test_seq = "TAGGAGCGAGCTTGCTCGCGATGGTCGTCAACGATAACGCGCCCAACCTGGGTGAATGCGCTGTCTTGACGTTTTTCGCGAGCAAGCTCGCTCCTA"
    blast_out = search_blastdb(test_seq)
    print(blast_out)


if __name__ == "__main__":
    st = time.time()
    main()
    print("Runtime: {:.2}s".format(time.time() - st))
