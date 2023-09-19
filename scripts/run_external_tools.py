import os

from scripts import option_parser
from scripts.utils import run_external_tool


def run_fastp(l_of_short_reads, common_output_dir, num_threads, min_length):
    fastp_output_dir = common_output_dir + 'fastp_output/'
    os.mkdir(fastp_output_dir)
    cmd_fastp = [option_parser.fpath_fastp]
    if num_threads > 16:
        n_threads = '16'
    else:
        n_threads = str(num_threads)

    for elem_idx in range(len(l_of_short_reads)):
        out_name_short_reads = fastp_output_dir + 'trimmed.' + os.path.basename(l_of_short_reads[elem_idx])
        cmd_fastp.extend([f'--in{str(elem_idx + 1)}', l_of_short_reads[elem_idx], f'--out{str(elem_idx + 1)}', out_name_short_reads])

    run_external_tool(' '.join(cmd_fastp) + f' --cut_front --cut_right -l {min_length} --thread {n_threads} -h {fastp_output_dir}info.html -j {fastp_output_dir}info.json', common_output_dir)
    return fastp_output_dir


def run_porechop(long_reads, common_output_dir, num_threads):
    porechop_output_dir = common_output_dir + 'porechop_output/'
    os.mkdir(porechop_output_dir)
    cmd_porechop = f'{option_parser.fpath_porechop} -t {str(num_threads)} -i {long_reads} -o {porechop_output_dir}trimmed.{os.path.basename(long_reads)}'
    run_external_tool(cmd_porechop, common_output_dir)
    return porechop_output_dir


def run_rnaspades(l_of_short_reads, long_reads, memory_lim, num_threads, common_output_dir):
    rnaspades_output_dir = common_output_dir + 'rnaspades_output/'
    cmd_rna_spades = ['python', option_parser.fpath_rnaspades, '-m', str(memory_lim), '-t', str(num_threads), '--nanopore', long_reads, '-o', rnaspades_output_dir]
    if len(l_of_short_reads) == 1:
        flag_prefix = '-s'
    else:
        flag_prefix = ''

    for elem_idx in range(len(l_of_short_reads)):
        cmd_rna_spades.extend([f'-{flag_prefix}{str(elem_idx + 1)}', os.path.abspath(l_of_short_reads[elem_idx])])

    run_external_tool(' '.join(cmd_rna_spades), common_output_dir)
    return rnaspades_output_dir


def run_diamond_blastp(fpath_db_diamond, common_output_dir, output_postfix_dir, d_of_optional_args,
                       query_fpath, num_threads, output_filename, taxon_list=None):
    if common_output_dir[-1] != '/':
        common_output_dir += '/'

    diamond_output_dir = common_output_dir + output_postfix_dir
    os.mkdir(diamond_output_dir)
    evalue, max_target_seqs, outfmt = d_of_optional_args['evalue'], d_of_optional_args['max_target_seqs'], d_of_optional_args['outfmt']

    cmd_diamond_blastp = f'{option_parser.fpath_diamond} blastp --db {fpath_db_diamond} --query {query_fpath} --evalue {evalue} --max-target-seqs {str(max_target_seqs)} --threads {str(num_threads)} --out {diamond_output_dir}{output_filename} --outfmt {outfmt}'
    if taxon_list is not None:
        cmd_diamond_blastp += f' --taxonlist {taxon_list}'
    run_external_tool(cmd_diamond_blastp, common_output_dir)

    return diamond_output_dir + output_filename


def run_transdecoder_longorfs(common_output_dir, rnaspades_output_dir):
    transdecoder_output_dir = common_output_dir + 'transdecoder_output'
    cmd_transdecoder_longorfs = f'{option_parser.fpath_transdecoder}TransDecoder.LongOrfs -t {rnaspades_output_dir}transcripts.fasta --output_dir {transdecoder_output_dir}'
    run_external_tool(cmd_transdecoder_longorfs, common_output_dir)
    return transdecoder_output_dir


def run_transdecoder_predict(common_output_dir, rnaspades_output_dir,
                             transdecoder_output_dir, diamond_file_for_searching_homologous_prots):
    cmd_transdecoder_predict = f'{option_parser.fpath_transdecoder}TransDecoder.Predict -t {rnaspades_output_dir}transcripts.fasta --retain_blastp_hits {diamond_file_for_searching_homologous_prots} --output_dir {transdecoder_output_dir}'
    run_external_tool(cmd_transdecoder_predict, common_output_dir)


def run_cdhit(memory_lim, num_threads, common_output_dir, extraction_after_transdecoder_output_dir,
              prop_idy_seqs, type_alignment, subseq_len_matching_frac):
    cdhit_output_dir = common_output_dir + 'clustered_cdhit/'
    os.mkdir(cdhit_output_dir)
    cmd_cdhit_est = f'{option_parser.fpath_cdhit} -M {str(memory_lim * 1000)} -c {prop_idy_seqs} -G {type_alignment} -aS {subseq_len_matching_frac} -p 1 -g 1 -T {str(num_threads)} -i {extraction_after_transdecoder_output_dir}cds_longest_iso.fasta -o {cdhit_output_dir}clust_cds_longest_iso.fasta'
    run_external_tool(cmd_cdhit_est, common_output_dir)
    return cdhit_output_dir
