# LibraryQC

## Description

A script for sgRNA extraction and quality control.

## Getting started

Clone this repo:

```shell
git clone https://github.com/Jianhua-Wang/LibraryQC
```

set up conda environment and activate the environment:

```shell
cd LibraryQC
conda env create -f libraryqccondaenv.yml
conda activate libraryqc
```

## Usage

```shell
(base) ➜  LibraryQC git:(master) ✗ python libraryqc.py -h    
usage: Librayqc using pattern matching to extract the sgRNA from pair-end sequencing data.

python libraryqc.py -i input/file_path.csv -o output -v input/vector.fa -l
input/library.csv

optional arguments:
  -h, --help       show this help message and exit
  -i , --input     csv file contains the sampl name and fq1, fq2 file path of
                   samples
  -o , --output    directory of output file
  -v , --vector    the fasta file of your vector, containing two sequences, 5'
                   and 3'.
  -l , --library   csv file of sgRNA you designed
  -s , --shot      the number of bases upstream and downstream of sgRNA you
                   exptected to match. More shot might reduce the false
                   positive. default=4
```

## Input

### 1. Path of fastq files

The path of fastq files are specified in a **csv** file, for example:

| sample_name | fq_1                       | fq_2                       |
| ----------- | -------------------------- | -------------------------- |
| day0        | ./input/data/day0_1.fq.gz  | ./input/data/day0_1.fq.gz  |
| day14       | ./input/data/day14_1.fq.gz | ./input/data/day14_1.fq.gz |
| vf          | ./input/data/vf_1.fq.gz    | ./input/data/vf_1.fq.gz    |

### 2. sgRNA library

The library you designed or download from other website, for example:

| id   | sequence             | note |
| ---- | -------------------- | ---- |
| 1    | ATAGGCACACATGAAGCGGA |      |
| 2    | TTTGCTGATAACTAGATCTA |      |
| 3    | TTGCAGGCCGCGATCTGTGC |      |
| 4    | GGTGGAGACTCCGAGTGTAG |      |
| 5    | ACTTCAAGTACGAGAACCAG | -    |
| 6    | AATTGCATTCGGTTTCTATC | -    |
| 7    | GAGCTTTTTGGGGTGTGACC | +    |
| 8    | AGCCGGTGTTAGTAAGAAAT | +    |
| ...  | ...                  | ...  |

### 3. Sequence of vector

For non-directional sequencing library, we need to know the sequence vector without sgRNA to determine the strand of extracted sgRNA. The sequence should be split into two end at the insertion site of sgRNA, one end should be named as 5 and another is 3.

```
>5
ccagaGagggcctatttcccatgattccttcatatttgcatatacgatacaaggctgttagagagataattagaattaatttgactgtaaacacaaagatattagtacaaaatacgtgacgtagaaagtaataatttcttgggtagtttgcagttttaaaattatgttttaaaatggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccg
>3
gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcg
```

## Output

Extracted sgRNA for each sample and an HTML report for quality control.

The HTML report was powered by [Multiqc](<https://multiqc.info/>), see the example in `output` directory.

## Note

1. This script only suits pair-end non-directional sequencing data and the insertion site of sgRNA is near 3' end.