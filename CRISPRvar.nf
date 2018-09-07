#!/usr/bin/env nextflow

params.fastq1 = '/f/mulinlab/jianhua/CRISPRvar/input/test_1.fq'
params.fastq2 = '/f/mulinlab/jianhua/CRISPRvar/input/test_2.fq'

params.output_dir = '/f/mulinlab/jianhua/CRISPRvar/output/'

params.pattern = '([ATCG]*ACCG)([ACTG]{20})GT{2,4}AGAGC[ATCG]*'
params.library = '/f/mulinlab/jianhua/CRISPRvar/input/library.csv'
params.vector = '/f/mulinlab/jianhua/CRISPRvar/input/vector.txt'

fq1 = file(params.fastq1)
fq2 = file(params.fastq2)
pattern_str = params.pattern
library_csv = file(params.library)
vector_txt = file(params.vector)

process build_vector_index {
    input:
        file vector from vector_txt
    output:
        file 'vector.*' into vector_ch
    script:
        template 'build_vector_index.py'
}

process build_lib_index {
    input:
        file library from library_csv
    output:
        file 'lib.*' into lib_ch
    script:
        template 'build_lib_index.py'
}
/*
process fastqtosam {
    input:
        file fq1 from fq1
        file fq2 from fq2
        publishDir params.output_dir, pattern: "*.sam", overwrite:true, mode:'copy'
    output:
        file '*.sam' into samwith2fq
    script:
        m = fq1 =~ /(.*)(_1)(.*)/
        base_name = m[0][1]
        """
        picard FastqToSam F1=${fq1} F2=${fq2} O=fastqtosam.sam SM=${base_name} RG=${base_name}
        """
}
*/
process extract_sg_from_sam {
    input:
        file fq1 from fq1
        file fq2 from fq2
        publishDir params.output_dir, pattern: "*.{txt,fastq}", overwrite:true, mode:'copy'
    output:
        file 'extract.txt' into extract_sgrna
        file 'fail_pattern_mapping.fastq' into fail_pattern_mapping
        file 'succeed_pattern_mapping.fastq' into succeed_pattern_mapping_fq
        file 'two_in_one.txt' into two_in_one_txt
    script:
        template 'extract_sg_from_sam.py'
}

process maptolib {
    input:
        file lib from lib_ch.collect()
        file fastq from succeed_pattern_mapping_fq
        publishDir params.output_dir, pattern: "*.sam", overwrite:true, mode:'copy'
    output:
        file 'maptolib.sam' into maptolib_sam
    script:
        """
        bowtie -q -S -p 8 -l 20 -n 3 lib ${fastq} maptolib.sam --norc
        """
}

process count_mismatch {
    input:
        file sam from maptolib_sam
        publishDir params.output_dir, pattern: "*.txt", overwrite:true, mode:'copy'
    output:
        file 'mismatch.txt' into mismatch_txt
        file 'pbc_qc.txt' into pbc_qc
    script:
        template 'count_mismatch.py'
}

process maptovector {
    input:
        file vector from vector_ch.collect()
        file fastq from fail_pattern_mapping
        publishDir params.output_dir, pattern: "*.{txt,sam}", overwrite:true, mode:'copy'
    output:
        file 'maptovector.txt' into maptovector_txt
        file 'maptovector.sam' into maptovector_sam
    script:
        """
        bowtie2 -x vector -p 8 -U ${fastq} -S maptovector.sam
        cat .command.log > maptovector.txt
        """
}

process calculate_composition {
    input:
        file mismatch from mismatch_txt
        file fastq from fq1
        file extract from extract_sgrna
        file distance_58 from succeed_pattern_mapping_fq
        file maptovector from maptovector_sam
        file lib_csv from library_csv
        publishDir params.output_dir, pattern: "*.txt", overwrite:true, mode:'copy'
    output:
        file 'ngs_composition.txt' into ngs_composition
        file 'lib_composition.txt' into lib_composition
    script:
        template 'calculate_composition.py'
}

process mismatch_analysis {
    input:
        file txt from mismatch_txt
        publishDir params.output_dir, pattern: "*.png", overwrite:true, mode:'copy'
    output:
        file '01pie.png' into pie_png
        file '02his.png' into his_png
        file '03distribution.png' into distribution_png
        file '04norm.png' into norm_png
    script:
        template 'mismatch_analysis.R'
}

process generate_report {
    input:
        file ngs from ngs_composition
        file lib from lib_composition
        file pbc from pbc_qc
        file pie from pie_png
        file his from his_png
        file distribution from distribution_png
        file norm from norm_png
        file fastq1 from fq1
        file fastq2 from fq2
        file library from library_csv
        file vector from vector_txt
        publishDir params.output_dir, pattern: "*.html", overwrite:true, mode:'copy'
    output:
        file 'LibraryQC_Report.html'
    script:
        template 'generate_report.py'
}