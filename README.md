# REPIN CLUSTERER

## Installation
The tool can be installed with the command:

    pip install repinclusterer

## Introduction
repin_clusterer is a tool developed for the clustering of REPIN sequences based on their position in the genome. REPIN Sequences that are present in orthologous regions will be placed in the same cluster. The orthology is determined based on the similarity of the nucleotide sequence that flanks either side of the REPIN.

![rocess of REPIN Clustering](./readme_images/repin_process.jpeg)

## Using this tool
The tool can be run with the command:

    repinclusterer --repin all_repins.txt --genomes genomes/ -out output_dir
  | Tag      |      Function      |
|:------------:|:-------------:|
  |repin|path to the text file containing REPIN sequences|
  |genomes|path to the directory containing all the genome sequences|
  |out|path to where the output file should be stored|
  |win|begin flanking region after, default  = 250. This means that if the REPIN starts at 'x' and ends at 'y', the flanking region will begin at 'x-250' and from 'y+250'|
  |fsize|length of the flanking region to consider, default = 1000|
  |pident|percentage sequence similarity that needs to be met, default= 90%|
  |coverage|minimum length of sequence that has to align/match, default = 90%|
  <br>
  **Note**: Each tag in the command begins with two '-' dash characters followed by a space and then the argument (see example above)

![Clustering Parameters](./readme_images/repin_flank.png)

## Input Format
The input files that are required are by software:
<ol><li> List of REPIN Sequences with position </li>
<li> List of all genomes whose REPINs are provided. </li>
</ol>
**Keep in mind** <br>
<ol>
<li> If there are REPINs whose genome sequence files are not provided, those REPINs will be dropped from the analysis</li>
<li>The file containing REPIN sequences should be formatted such that each line contains (only) the following information:<br>
genome_name repin_start repin_end repin_type repin_sequence
<br>Ex:<br>
chlTAMOak81 1008421 1008530 green AGCTATCGTAC.......GTACGATAGCT </li>
<li>It is preferrable to provide the genome sequences in fasta format. </li>
</ol>

## Output Format
The output file is very similar to the input file with the addition of a number at the beginning of the link representing the cluster number.

> num genome_name repin_start repin_end repin_type repin_sequence<br>
> 0 chlTAMOak81 1008421 1008530 green AGCTATCGTAC.......GTACGATAGCT

Implying that this particular REPIN belongs to cluster number 0 and so on.

## Python packages used
<br>**Biopython**
<br>Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. _Bioinformatics_, _25_(11), 1422–1423.
<br>**Pickle**
<br>Van Rossum, G. (2020). _The Python Library Reference, release 3.8.2_. Python Software Foundation.
<br>**Networkx**
<br>Hagberg, A., Swart, P., & S Chult, D. (2008). _Exploring network structure, dynamics, and function using NetworkX_.
