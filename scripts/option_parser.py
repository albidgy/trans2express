import argparse
import os
from datetime import datetime


def parser_arguments():
    parser = argparse.ArgumentParser(description='Trans2express: de novo hybrid transcriptome assembly tool')
    general_args = parser.add_argument_group(description='General arguments')
    general_args.add_argument('-1',
                              '--short_reads1',
                              required=True,
                              help='Short reads 1 in fastq or fastq.gz format',
                              )
    general_args.add_argument('-2',
                              '--short_reads2',
                              default=None,
                              help='Short reads 2 in fastq or fastq.gz format. If you have single-end reads do not specify this parameter',
                              )
    general_args.add_argument('--long_reads',
                              required=True,
                              help='Long reads in fastq or fastq.gz format',
                              )
    general_args.add_argument('-t',
                              '--threads',
                              default=1,
                              type=int,
                              help='Number of threads [default: 1]',
                              )
    general_args.add_argument('-o',
                              '--output_dir',
                              default='../res_trans2express_' + datetime.now().strftime('%Y_%m_%d_%H_%M_%S') + '/',
                              help='Output directory [default: ../res_trans2express_YEAR_MONTH_DAY_HOUR_MINUTE_SECOND]',
                              )
    general_args.add_argument('-c',
                              '--config_file',
                              default='CONFIGURATIONS.txt',
                              help='Path to configurations file',
                              )
    general_args.add_argument('--trim_short_reads',
                              default=False,
                              action='store_true',
                              help='Add trimming reads step by using fastp tool [default: False]',
                              )
    general_args.add_argument('--trim_long_reads',
                              default=False,
                              action='store_true',
                              help='Add trimming long reads step by using Porechop tool [default: False]',
                              )
    general_args.add_argument('-m', '--memory_lim',
                              default=128,
                              type=int,
                              help='Memory limit in Gb [default: 128]',
                              )

    optional_args = parser.add_argument_group(description='Optional arguments')
    optional_args.add_argument('--diamond_db',
                              default='db/nr.dmnd',
                              help='DIAMOND protein database nr.dmnd for finding homologous proteins for prediction CDS by TransDecoder and for removing foreign transcripts. The database is downloaded when you run the install.sh script or you can create your own .dmnd db [default: db/nr.dmnd]',
                              )
    optional_args.add_argument('--diamond_taxonomic_id',
                              default='db/taxonomic_id_to_full_taxonomy.txt',
                              help='File with taxonomic ids list by nr.dmnd database for removal foreign rna. The file is downloaded when you run the install.sh script [default: db/taxonomic_id_to_full_taxonomy.txt]',
                              )
    optional_args.add_argument('--go_tree',
                              default='db/goTree.txt',
                              help='File with broad GO terms. The file is downloaded when you run the install.sh script or you can create your own goTree file [default: db/goTree.txt]',
                              )
    optional_args.add_argument('--min_short_read_length',
                               default=50,
                               type=str,
                               help='Minimum length of short reads for fastp tool [default: 50]',
                               )
    optional_args.add_argument('--seq_idy_threshold',
                               default=0.98,
                               type=str,
                               help='Sequence identity threshold for CD-HIT-EST tool [default: 0.98]',
                               )
    optional_args.add_argument('--alignment_type',
                               default=0,
                               type=str,
                               help='Select type of alignment for CD-HIT-EST tool. 0 - local alignment, 1 - global alignment [default: 0]',
                               )
    optional_args.add_argument('--subseq_len_matching_cov',
                               default=0.6,
                               type=str,
                               help='Alignment coverage for the shorter sequence for CD-HIT-EST tool [default 0.6]'
                               )
    return parser.parse_args()


def read_configuration_file(config_file):
    with open(config_file) as inf:
        lines = inf.readlines()
        fpath_fastp = lines[4].strip('\n')
        fpath_porechop = lines[6].strip('\n')
        fpath_rnaspades = lines[8].strip('\n')
        fpath_transdecoder = lines[10].strip('\n')
        fpath_cdhit = lines[12].strip('\n')
        fpath_diamond = lines[14].strip('\n')
    return [fpath_fastp, fpath_porechop, fpath_rnaspades,
            fpath_transdecoder, fpath_cdhit, fpath_diamond]


args_for_searching_homologous_prots = {'evalue': '0.001',
                                       'max_target_seqs': 1,
                                       'outfmt': '6',
                                       }
args_for_removal_foreign_rna = {'evalue': '1e-5',
                                'max_target_seqs': 8,
                                'outfmt': '6 qseqid sseqid staxids',
                                }
