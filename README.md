# REPORTH

## List of Contents

- [Installation](#installation)
- [Introduction](#introduction)
- [Using this tool](#using-this-tool)
- [Example Dataset](#example-dataset)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [Python packages used](#python-packages-used)

## Installation
The tool can be installed with the command:

    pip install reporth

Note: `makeblastn` and `blastn` tools should be installed and in your $PATH

## Introduction
repin_clusterer is a tool developed for the clustering of REPIN sequences based on their position in the genome. REPIN Sequences that are present in orthologous regions will be placed in the same cluster. The orthology is determined based on the similarity of the nucleotide sequence that flanks either side of the REPIN.
<br>See also: rarefan.evolbio.mpg.de

![Process of REPIN Clustering](./readme_images/repin_process.heic)

## Using this tool
The tool can be run with the command:

    reporth --repin all_repins.txt --genomes genomes/ --out output_dir --reptypes 1,2
  | Tag      |      Function      |      Default      |
|:------------:|:-------------:|:-------------:|
  |repin|path to the text file containing REPIN sequences or RAREFAN Dir| None|
  |genomes|path to the directory containing all the genome sequences| None|
  |out|path to where the output file should be stored| ./cluster\_output\_[date]_[run-id]/|
  |reptypes|Mention the specific repin types to accept from rarefan output | None|
  |win|begin flanking region after. This means that if the REPIN starts at 'x' and ends at 'y', the flanking region will begin at 'x-250' and from 'y+250'| 250 |
  |fsize|length of the flanking region to consider|1000|
  |pident|percentage sequence similarity that needs to be met|90(%)|
  |coverage|minimum length of sequence that has to align/match|90(%)|
  <br>
  **Note**:
  1. Each tag in the command begins with two '-' dash characters followed by a space and then the argument (see example above)
  2. For --reptypes, the numbers should be written separated by a comma without a space between the values

![Clustering Parameters](./readme_images/repin_flank.png)

## Example Dataset
We have created an example dataset for you to understand how to use the software and for a demonstration.
<br>Download the sample input dataset from the link below.
<br>Install `reporth` and run using the command as mentioned above.
<br>The output will be stored in the `cluster_output_[date]_[run-id]` directory in the format as mentioned below.
<br>[Download Test Dataset](https://github.com/blackthorne18/reporth_cli/files/8309764/test_data.zip)

## Input Format
### RAREFAN Output
When a RAREFAN output directory is provided as input REPORTH will go through the RAREFAN output folder and find all genomes directories ([genome]\_[x]/) where x is any digit signifying the type of REP. Using the ‘--reptype' tag, the range of x can be specified. Within these folders the locations listed in the files [genome]_[x].ss are parsed. 
> reporth --repin rarefan_output_dir --genomes genomes --reptypes 0,1,2

### List of REPINs as input
The input files that are required are by software:
<ol><li> List of REPIN Sequences with position </li>
<li> List of all genomes whose REPINs are provided. </li>
</ol>
#### Note
<ol>
<li> If there are REPINs whose genome sequence files are not provided, those REPINs will be dropped from the analysis</li>
<li>The file containing REPIN sequences should be formatted such that each line contains (only) the following information:<br>
genome_name repin_start repin_end repin_type repin_sequence
<br>Ex: `TAMOak81 1008421 1008530 type0<br>Pb-St2 1008000 1008123 type0` </li>
<li>It is preferrable to provide the genome sequences in fasta format. </li>
</ol>
> reporth --repin listofrepins.txt --genomes genomes --reptypes 0,1,2

## Output Format
The primary output file `clusters_[date].txt` is very similar to the input file with the addition of a number at the beginning of the link representing the cluster number.
> num genome_name repin_start repin_end repin_type<br>
> 0 chlTAMOak81 1008421 1008530 type0

Implying that this particular REPIN belongs to cluster number 0 and so on.

| Output File               | Description of contents                                                                                                                                                             |
|---------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| clusters_[date].txt       | Listing all repetitive elements with unique cluster ID                                                                                                                              |
| meta_cluster_[date].txt   | Lists the sequence similarity of the flanking sequences of a REPIN with the flanking sequences of each other REPIN in a cluster (in the order they are present in the cluster file) |
| lhs_hits.p/rhs_hits.p     | Pickle(python package) file storing the BLAST results of all flanking sequences                                                                                                     |
| store_nearby_repins.p     | Pickle(python package) file storing the flanking sequences and their corresponding BLAST hits and the presence of REPINs near these flanking sequences                              |
| flanking_pairwise_dists.p | Stores the pairwise distances between all flanking sequences from the BLAST hits                                                                                                    |
| path_making_[date].txt    | File listing whether a REPIN was merged into a cluster because of left flanking sequence/right flanking sequence or both                                                            |

## Python packages used
**Biopython**
<br>Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. _Bioinformatics_, _25_(11), 1422–1423.
<br>**Pickle**
<br>Van Rossum, G. (2020). _The Python Library Reference, release 3.8.2_. Python Software Foundation.
<br>**Networkx**
<br>Hagberg, A., Swart, P., & S Chult, D. (2008). _Exploring network structure, dynamics, and function using NetworkX_.
