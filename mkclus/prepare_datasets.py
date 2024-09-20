# -------------------------------------------------------------------------------------------
"""
Input Files Required:
1.  List of genomes
    Path = ./input/genomes/*.fas
2.  Repins File
    Path = ./input/repins.p

Output Files:
1.  blastgenomedb
    Path = bank_path/genomes_blastdb/
2.  repin_with_1k_flanks
    Path = bank_path/repin_with_1k_flanks.p

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

blast_path = "/genomes_blastdb/"
temp_files = "/dumpyard/"
all_parameters_path = "/all_parameters.p"


def setup_blastdb(bank_path):
    global blast_path, temp_files, all_parameters_path

    blast_path = bank_path + blast_path
    temp_files = bank_path + temp_files
    all_parameters_path = bank_path + all_parameters_path

    all_parameters = pickle.load(open(all_parameters_path, "rb"))
    genomes_list = all_parameters['genomes']

    repins = open(all_parameters['repin'], "r").read().split("\n")

    flank_gene_range = {
        'window': all_parameters['win'], 'flanklength': all_parameters['fsize']}
    fgr = flank_gene_range['window'] + flank_gene_range['flanklength']

    genome_sequences = []
    repin_with_1k_flanks = {}
    repin_per_genome = {genome.split("/")[-1].split(".")[0]: []
                        for genome in genomes_list}

    for rep in repins:
        if len(rep) < 1:
            continue
        splitrep = rep.split()
        splitrep[1] = int(splitrep[1])
        splitrep[2] = int(splitrep[2])
        repin_per_genome[splitrep[0]].append(splitrep)

    for genome in genomes_list:
        sequence = str(SeqIO.read(genome, "fasta").seq)
        genname = genome.split("/")[-1].split(".")[0]
        for rep in repin_per_genome[genname]:
            left_flank = sequence[rep[1] -
                                  fgr: rep[1] - flank_gene_range['window']]
            right_flank = sequence[rep[2] +
                                   flank_gene_range['window']: rep[2] + fgr]
            repname = rep[0] + " " + str(rep[1]) + " " + str(rep[2])
            if len(left_flank) == 0:
                left_flank = "N" * 1000
            if len(right_flank) == 0:
                right_flank = "N" * 1000
            repin_with_1k_flanks[repname] = [
                repname, rep[3], left_flank, rep[4], right_flank]

        genome_sequences.append(">{}\n{}".format(
            genname, sequence))

    pickle.dump(repin_with_1k_flanks, open(
        f"{bank_path}/repin_with_1k_flanks.p", "wb"))

    genome_sequences = "\n".join(genome_sequences)
    open(blast_path + "allgenomes.fas", "w").write(genome_sequences)

    # Using BLAST CLI
    cmd = f"makeblastdb -in {blast_path}allgenomes.fas -out {blast_path}allgenomes -dbtype nucl"
    os.system(cmd)
    # DEPRACATED - Using Biopython
    # cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=blast_path + "allgenomes.fas", out=blast_path + "allgenomes")
    # cline()


def search_blastdb(bank_path, sequence_list, flank_gene_param):
    infile = temp_files + "test1_in.fas"
    outfile = temp_files + "test1_out.fas"
    with open(infile, "w+") as f:
        for key, val in sequence_list.items():
            f.write(f">{key}\n{val}\n")

    # Using BLAST CLI
    cmd_format = "6 qseqid sseqid pident length sstart send"
    cmd = f"blastn -query {infile} -db {blast_path+'allgenomes'} -out {outfile} -outfmt '{cmd_format}'"
    os.system(cmd)
    # DEPRACATED - Using Biopython
    # cline = NcbiblastnCommandline(query=infile, db=blast_path "allgenomes", out=outfile, outfmt=cmd_format)
    # cline()

    outfile = [i.split("\t") for i in open(
        outfile, "r").read().split("\n") if len(i) > 0]
    good_output = []
    for i in range(len(outfile)):
        seq_name = outfile[i][0]
        outfile[i] = [seq_name, outfile[i][1]] + \
            [int(float(x)) for x in outfile[i][2:]]
        lengthmatch = int(
            100 * (abs(outfile[i][5] - outfile[i][4]) / len(sequence_list[seq_name])))
        pident = outfile[i][2]
        if lengthmatch >= flank_gene_param['lengthmatch'] and pident >= flank_gene_param['pident']:
            good_output.append(outfile[i])

    # Reformat output so that we have a key-value pair for each query and each hit
    refmt = {}
    for i in range(len(good_output)):
        query_name = good_output[i][0]
        if query_name not in refmt.keys():
            refmt[query_name] = []
        refmt[query_name].append(good_output[i])

    # Makes sure that only one hit is recorded per genome and this hit is the highest hit
    seen_gens = {}
    to_keep = []
    for key, gop in refmt.items():
        for i in range(len(gop)):
            gen = gop[i][1]
            if gen not in seen_gens.keys():
                seen_gens[gen] = [0, 0]
            if gop[i][2] > seen_gens[gen][0] and gop[i][3] > seen_gens[gen][1]:
                to_keep.append(i)
                seen_gens[gen] = [gop[i][2], gop[i][3]]

        gop = [gop[i] for i in range(len(gop)) if i in to_keep]
        refmt[key] = gop

    return refmt


def main():
    print("Why is this being run?")


if __name__ == "__main__":
    st = time.time()
    main()
    print("Runtime: {:.2}s".format(time.time() - st))
