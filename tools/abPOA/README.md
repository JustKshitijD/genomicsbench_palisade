# abPOA: adaptive banded Partial Order Alignment
[![Latest Release](https://img.shields.io/github/release/yangao07/abPOA.svg?label=Release)](https://github.com/yangao07/abPOA/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/abPOA/total.svg?label=Download)](https://github.com/yangao07/abPOA/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abpoa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/abpoa)
[![PyPI](https://img.shields.io/pypi/dm/pyabpoa.svg?label=pip%20install)](https://pypi.python.org/pypi/pyabpoa)
[![Published in bioRxiv](https://img.shields.io/badge/Preprint-bioRxiv-red.svg)](https://doi.org/10.1101/2020.05.07.083196)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/abPOA.svg?label=Issues)](https://github.com/yangao07/abPOA/issues)
[![Build Status](https://img.shields.io/travis/yangao07/abPOA/master.svg?label=Master)](https://travis-ci.org/yangao07/abPOA)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)
<!-- [![PyPI](https://img.shields.io/pypi/v/pyabpoa.svg?style=flat)](https://pypi.python.org/pypi/pyabpoa) -->
## Updates (v1.0.4)

- Added read ID as head in MSA output: `-A`
- Added GFA output: `-r3`/`-r4`
- Added ambiguous strand mode: `-s`

## Getting started
Download the [latest release](https://github.com/yangao07/abPOA/releases):
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.4/abPOA-v1.0.4.tar.gz
tar -zxvf abPOA-v1.0.4.tar.gz && cd abPOA-v1.0.4
```
Make from source and run with test data:
```
make; ./bin/abpoa ./test_data/seq.fa > cons.fa
```
Or, install via conda and run with test data:
```
conda install -c bioconda abpoa
abpoa ./test_data/seq.fa > cons.fa
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing abPOA via conda](#conda)
  - [Building abPOA from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate consensus sequence](#gen_cons)
  - [To generate row-column multiple sequence alignment](#gen_msa)
  - [To generate graph information in GFA format](#gen_gfa)
  - [To generate a plot of the alignment graph](#gen_plot)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [Consensus sequence](#cons)
  - [Row-column multiple sequence alignment](#msa)
  - [Full graph information](#gfa)
  - [Plot of alignment graph](#plot)
- [For development](#dev)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
abPOA is an extended version of [Partial Order Alignment (POA](10.1093/bioinformatics/18.3.452)) 
that performs adaptive banded dynamic programming (DP) with an SIMD implementation. 
abPOA can perform multiple sequence alignment (MSA) on a set of input sequences and 
generate a consensus sequence by applying the [heaviest bundling algorithm](10.1093/bioinformatics/btg109) 
to the final alignment graph.

abPOA can generate high-quality consensus sequences from error-prone long reads and offer 
significant speed improvement over existing tools.

abPOA supports three alignment modes (global, local, extension) and flexible scoring schemes that allow linear, affine and convex gap penalties. 
It right now supports SSE2/SSE4.1/AVX2/AVX512F/AVX512BW vectorization and more advanced instructions 
will be supported in the future.

For more information, please refer to our [preprint paper](https://doi.org/10.1101/2020.05.07.083196).

## <a name="install"></a>Installation

### <a name="conda"></a>Installing abPOA via conda
On Linux/Unix and Mac OS, abPOA can be installed via
```
conda install -c bioconda abpoa   # install abPOA program
```

### <a name="build"></a>Building abPOA from source files
You can also build abPOA from source files. 
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/abPOA/releases).
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.4/abPOA-v1.0.4.tar.gz
tar -zxvf abPOA-v1.0.4.tar.gz
cd abPOA-v1.0.4; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of abPOA, which might be still under development.
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.4/abPOA-v1.0.4_x64-linux.tar.gz
tar -zxvf abPOA-v1.0.4_x64-linux.tar.gz
```

## <a name="usage"></a>General usage
### <a name="gen_cons"></a>To generate consensus sequence

```
abpoa seq.fa > cons.fa
```

### <a name="gen_cons"></a>To generate row-column multiple sequence alignment in PIR format

```
abpoa seq.fa -r2 > cons.out
```

### <a name="gen_gfa"></a>To generate graph information in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format

```
abpoa seq.fa -r3 > abpoa.gfa
```
To include the generated consensus sequence as a path in the GFA file:
```
abpoa seq.fa -r4 > abpoa.gfa
```

### <a name="gen_plot"></a>To generate a plot of the alignment graph

```
abpoa seq.fa -g poa.png > cons.fa
```
See [Plot of alignment graph](#plot) for more details about the plot file.

## <a name="cmd"></a>Commands and options
```
Usage: abpoa [options] <in.fa/fq> > cons.fa/msa.out/abpoa.gfa

Options:
  Alignment:
    -m --aln-mode INT       alignment mode [0]
                              0: global, 1: local, 2: extension
    -M --match    INT       match score [2]
    -X --mismatch INT       mismatch penalty [4]
    -O --gap-open INT(,INT) gap opening penalty (O1,O2) [4,24]
    -E --gap-ext  INT(,INT) gap extension penalty (E1,E2) [2,1]
                            abPOA provides three gap penalty modes, cost of a g-long gap:
                            - convex (default): min{O1+g*E1, O2+g*E2}
                            - affine (set O2 as 0): O1+g*E1
                            - linear (set O1 as 0): g*E1
    -s --amb-strand         ambiguous strand mode [False]
                            for each input sequence, try the reverse complement if the current
                            alignment score is too low, and pick the strand with a higher score
  Adaptive banded DP:
    -b --extra-b  INT       first adaptive banding parameter [10]
                            set b as < 0 to disable adaptive banded DP
    -f --extra-f  FLOAT     second adaptive banding parameter [0.01]
                            the number of extra bases added on both sites of the band is
                            b+f*L, where L is the length of the aligned sequence
  Input/Output:
    -l --in-list            input file is a list of sequence file names [False]
                            each line is one sequence file containing a set of sequences
                            which will be aligned by abPOA to generate a consensus sequence
    -o --output   FILE      ouput to FILE [stdout]
    -r --result   INT       output result mode [0]
                            - 0: consensus (FASTA format)
                            - 1: MSA (PIR format)
                            - 2: both 0 & 1
                            - 3: graph (GFA format)
                            - 4: graph with consensus path (GFA format)
    -A --msa-header         add read ID as header of each sequence in MSA output [False]
    -g --out-pog  FILE      dump final alignment graph to FILE (.pdf/.png) [Null]

    -h --help               print this help usage information
    -v --version            show version number
```

## <a name="input"></a>Input
abPOA works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is 
expected to contains multiple sequences which will be processed sequentially to perform the iterative 
sequence-to-graph (partial order) alignment.

abPOA can also take a list of file names as input with option `-l`, where each line is the path to one 
file containing multiple sequences. Each sequence file is then individually aligned by abPOA to generate a
consensus sequence.

## <a name="output"></a>Output
### <a name="cons"></a>Consensus sequence 
By default, abPOA only outputs the consensus sequence generated from the final alignment agraph.
It is in FASTA format with the name field set as "Consensus_sequence".
For example:
```
>Consensus_sequence
ACGTGTACACGTTGAC
```
### <a name="msa"></a>Row-column multiple sequence alignment
abPOA can also output the row-column multiple sequence alignment (RC-MSA) of all the aligned sequences in PIR format
with an additional FASTA header `>Multiple_sequence_alignment`. For example:
```
>Multiple_sequence_alignment
ACGTGTACA-GTTGAC
A-G-GTACACGTT-AC
A-GTGT-CACGTTGAC
ACGTGTACA--TTGAC
```
The `-` in the sequence stands for alignment gap. 

### <a name="gfa"></a>Full graph information
abPOA can output the final alignment graph in GFA format.
Each segment line (`S` line) represents one node and each link line (`L` line) represents one edge between two nodes.
The original input sequences and the generated consensus sequence are described as paths in `P` lines.

abPOA outputs two graph-related numbers in the header line (`H` line):
`NS` and `NL`, which denote the total number of nodes and edges in the GFA file, respectively.

Please refer to the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for more details of the GFA format.

### <a name="plot"></a>Plot of alignment graph

abPOA can generate a plot of the final partial order alignment graph with the help of `graphviz dot`. 
For example:

![pog](https://github.com/yangao07/abPOA/blob/master/pog.png)

The numbers inside the nodes are the node IDs. The numbers on the edges are the edge weights.
`S` and `E` are the auxiliary start and end nodes that have no sequence bases.

Make sure you have `dot` installed beforing using abPOA to generate the plot.
For Linux/Unix systems: `sudo apt-get install graphviz`.

## <a name="dev"></a>For development
abPOA is not only a stand-alone tool for MSA and consensus calling, it can also work as a programming library. [example.c](example.c) shows how to use the C APIs of abPOA to take a set of sequences as input and perform MSA and consensus calling. Basically, the library file `libabpoa.a` and two header files [abpoa.h](include/abpoa.h) and [simd_instruction.h](include/simd_instruction.h) are needed to make the abPOA library work in your program.

abPOA also provides Python bindings to all the primary C APIs. Refer to [python/README.md](python/README.md) for more details.

## <a name="contact"></a>Contact
Yan Gao yangao07@hit.edu.cn

Yi Xing XINGYI@email.chop.edu

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/abPOA/issues)
