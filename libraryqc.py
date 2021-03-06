#! /usr/bin/env python3
# coding: utf-8

import mappy as mp
import pandas as pd
import re, tempfile, shutil, os, argparse, subprocess
from subprocess import call
from multiprocessing import Pool


def print_logo():
    logo = '''
==================================================================
     \033[1;33m/\\\033[0m
    \033[1;33m/__\\\033[0m\033[1;31m\\\033[0m         This is a script for sgRNA Library QC
   \033[1;33m/\033[0m  \033[1;31m---\\\033[0m        Author: Jianhua Wang
  \033[1;33m/\\\033[0m      \033[1;31m\\\033[0m       Date:   06-10-2019
 \033[1;33m/\033[0m\033[1;32m/\\\033[0m\033[1;33m\\\033[0m     \033[1;31m/\\\033[0m
 \033[1;32m/  \   /\033[0m\033[1;31m/__\\\033[0m
\033[1;32m`----`-----\033[0m
==================================================================
    '''
    print(logo)


def parseArguments():
    parser = argparse.ArgumentParser(
        usage="Librayqc using pattern matching to extract the sgRNA from pair-end sequencing data.",
        description="python libraryqc.py -i input/file_path.csv -o output -v input/vector.fa -l input/library.csv", )
    parser.add_argument('-i', '--input', type=str,
                        help='csv file contains the sampl name and fq1, fq2 file path of samples', metavar=''),
    parser.add_argument('-o', '--output', type=str, help='directory of output file', metavar=''),
    parser.add_argument('-v', '--vector', type=str,
                        help='the fasta file of your vector, containing two sequences, 5\' and 3\'.', metavar=''),
    parser.add_argument('-l', '--library', type=str, help='csv file of sgRNA you designed', metavar=''),
    parser.add_argument('-s', '--shot', type=int,
                        help='the number of bases upstream and downstream of sgRNA you exptected to match. More shot might reduce the false positive. default=4',
                        default=4, metavar='')
    args = parser.parse_args()
    return args


lib_distribution_header = '''# id: lib_distribution
# section_name: \'Designed sgRNA library\'
# description: \': the number of designed sgRNA found and NOT found in sequencing data\'
# format: \'tsv\'
# plot_type: \'bargraph\'
# pconfig:
#    id: \'lib_distribution\'
#    title: Desinged sgRNA Bar Plot
#    ylab: \'Number of sgRNA\''''

read_distribution_header = '''# id: read_distribution
# section_name: \'Reads distribution\'
# description: \': the number of read counts in sequencing data\'
# format: \'tsv\'
# plot_type: \'bargraph\'
# pconfig:
#    id: \'read_distribution\'
#    title: Read counts Bar Plot
#    ylab: \'Read count\''''

sgrna_distribution_header = '''# id: sgrna_distribution
# section_name: \'Extracted sgRNA\'
# description: \': the number of sgRNA counts mapped or unmapped to the library\'
# format: \'tsv\'
# plot_type: \'bargraph\'
# pconfig:
#    id: \'sgrna_distribution\'
#    title: sgRNA counts Bar Plot
#    ylab: \'sgRNA count\''''

position_distribution_header = '''# id: position_distribution
# section_name: \'Position of mactched\'
# description: \': the distance between mactch site and 3 end of reads\'
# format: \'tsv\'
# plot_type: \'linegraph\'
# pconfig:
#    id: \'position_distribution\'
#    title: position distribution llne Plot
#    ylab: \'bp from 3 end\''''


def extract_sgrna_from_single_read(read):
    read_out = [0, 0, 0, 0]
    for hit in vector_aligner.map(read[1]):
        match = re.match(r'{}'.format(pattern), mp.revcomp(read[1]))
        if match:
            read_out = [read[0], match.group(2), read[2][len(match.group(3)):-len(match.group(1))][::-1], read[3],
                        len(match.group(3))]
        else:
            read_out = read
    return read_out


def run_pair_end(sample_name, fq1, fq2):
    sgrna_fq = open(f'{temp_dir}/{sample_name}_sgrna.fq', 'w')
    sgrna_txt = open(f'{temp_dir}/{sample_name}_sgrna.txt', 'w')
    print(f'Extract sgRNA from {sample_name}')
    for read1, read2 in zip(mp.fastx_read(fq1, read_comment=True), mp.fastx_read(fq2, read_comment=True)):
        read1, read2 = extract_sgrna_from_single_read(read1), extract_sgrna_from_single_read(read2)
        if len(read1) == 5:
            sgrna_txt.write(f'{read1[1]}\t{read1[4]}\t{1}\n')
            sgrna_fq.write(f'@{read1[0]}\n{read1[1]}\n+\n{read1[2]}\n')
        elif len(read2) == 5:
            sgrna_txt.write(f'{read2[1]}\t{read2[4]}\t{1}\n')
            sgrna_fq.write(f'@{read2[0]}\n{read2[1]}\n+\n{read2[2]}\n')
        elif read1[0] == 0 and read2[0] != 0:
            sgrna_txt.write(f'\t\t{3}\n')
        elif read2[0] == 0 and read1[0] != 0:
            sgrna_txt.write(f'\t\t{3}\n')
        elif read1[0] == 0 and read2[0] == 0:
            sgrna_txt.write(f'\t\t{4}\n')
        else:
            sgrna_txt.write(f'\t\t{2}\n')
    sgrna_fq.close()
    sgrna_txt.close()

    print(f'Count sgRNA quantity {sample_name}')
    sgrna_df = pd.read_csv(f'{temp_dir}/{sample_name}_sgrna.txt', sep='\t', names=['sgrna', 'position', 'lable'])
    sgrna_count = sgrna_df['sgrna'].value_counts()
    library_df = pd.read_csv(library)
    library_df['count'] = library_df['sequence'].map(sgrna_count)

    lib_distribution = pd.Series(data=[len(library_df[library_df['count'].notnull()]),
                                       len(library_df[library_df['count'].isnull()])],
                                 index=['Found', 'Not found'])

    sgrna_distribution = pd.Series(data=[library_df.drop_duplicates(subset=['sequence'])['count'].sum(),
                                         sgrna_count.sum() - library_df.drop_duplicates(subset=['sequence'])[
                                             'count'].sum()],
                                   index=['Mapped', 'Unmapped'])

    read_distribution = sgrna_df['lable'].value_counts()
    read_distribution.index = read_distribution.index.map(
        pd.Series(index=[1, 2, 3, 4], data=['matched', 'empty', 'one end unmapped', 'pair end unmapped']))

    position_distribution = sgrna_df['position'].value_counts()

    sgrna_plot = open(f'{temp_dir}/{sample_name}_sgrna_plot_mqc.txt', 'w')
    sgrna_plot.write(f'{sgrna_distribution_header}\n{sgrna_distribution.to_csv(sep=chr(9), header=False)}')
    sgrna_plot.close()
    lib_plot = open(f'{temp_dir}/{sample_name}_lib_plot_mqc.txt', 'w')
    lib_plot.write(f'{lib_distribution_header}\n{lib_distribution.to_csv(sep=chr(9), header=False)}')
    lib_plot.close()
    read_plot = open(f'{temp_dir}/{sample_name}_read_plot_mqc.txt', 'w')
    read_plot.write(f'{read_distribution_header}\n{read_distribution.to_csv(sep=chr(9), header=False)}')
    read_plot.close()
    read_plot = open(f'{temp_dir}/{sample_name}_position_plot_mqc.txt', 'w')
    read_plot.write(f'{position_distribution_header}\n{position_distribution.to_csv(sep=chr(9), header=False)}')
    read_plot.close()

    print(f'Run fastqc on extracted fastq {sample_name}')
    FNULL = open(os.devnull, 'w')
    call(f'fastqc {temp_dir}/{sample_name}_sgrna.fq -o {temp_dir}', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    shutil.move(f"{temp_dir}/{sample_name}_sgrna.fq", f"{output_dir}/{sample_name}_sgrna.fq")


if __name__ == "__main__":
    print_logo()
    args = parseArguments()

    output_dir = args.output
    vector = args.vector
    library = args.library
    fq_path = args.input
    shot = args.shot

    work_dir = os.getcwd()
    temp_dir = tempfile.mkdtemp()
    py_dir = os.path.split(os.path.realpath(__file__))[0]

    fq_path_df = pd.read_csv(fq_path)

    vector_sequence = mp.Aligner(vector)
    with open(f'{temp_dir}/vector1.fa', 'w') as vector_aligner:
        vector_aligner.write(f">vector\n{vector_sequence.seq('5') + vector_sequence.seq('3')}\n")
    vector_aligner = mp.Aligner(f'{temp_dir}/vector1.fa')
    pattern = '([ATCG]*' + vector_sequence.seq('5')[-shot:] + ')([ACTG]{19,20})(' + vector_sequence.seq('3')[
                                                                                    :shot] + '[ATCG]*)'

    p = Pool(len(fq_path_df))
    for i in fq_path_df.index:
        sample_name, fq1, fq2 = fq_path_df.loc[i].values
        p.apply_async(run_pair_end, (sample_name, fq1, fq2))
    p.close()
    p.join()

    os.chdir(temp_dir)
    print('Generate report')
    FNULL = open(os.devnull, 'w')
    call(f'multiqc -f --no-data-dir . -c {py_dir}/config.yml', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    os.chdir(work_dir)
    shutil.move(f"{temp_dir}/multiqc_report.html", f"{output_dir}/multiqc_report.html")
    shutil.rmtree(temp_dir)
    print('Completed.')
