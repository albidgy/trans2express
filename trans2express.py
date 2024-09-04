import os
import shutil
import sys
from datetime import datetime

from scripts import extract_longest_iso, go_annotation, run_external_tools, option_parser
from scripts.utils import Logger


def run_pipeline():
    arguments = option_parser.parser_arguments()
    if os.path.exists(arguments.output_dir):
        shutil.rmtree(arguments.output_dir)
    os.mkdir(arguments.output_dir)
    if arguments.output_dir[-1] != '/':
        arguments.output_dir = arguments.output_dir + '/'
    sys.stdout = Logger(arguments.output_dir)

    if arguments.short_reads2 is None:
        l_of_short_reads = [arguments.short_reads1]
    else:
        l_of_short_reads = [arguments.short_reads1, arguments.short_reads2]
    long_reads = arguments.long_reads

    # run trimming illumina reads
    if arguments.trim_short_reads:
        print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start trimming short reads by using fastp...')
        fastp_output_dir = run_external_tools.run_fastp(l_of_short_reads,
                                                        arguments.output_dir,
                                                        arguments.threads,
                                                        arguments.min_short_read_length)
        l_of_short_reads = [fastp_output_dir + file for file in os.listdir(fastp_output_dir) if file.startswith('trimmed.')]
        print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run trimming long reads
    if arguments.trim_long_reads:
        print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start trimming long reads by using Porechop...')
        porechop_output_dir = run_external_tools.run_porechop(long_reads,
                                                              arguments.output_dir,
                                                              arguments.threads)
        long_reads = f'{porechop_output_dir}trimmed.{os.path.basename(long_reads)}'
        print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run rnaSPAdes
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start transcriptome assembly by using rnaSPAdes...')
    rnaspades_output_dir = run_external_tools.run_rnaspades(l_of_short_reads,
                                                            long_reads,
                                                            arguments.memory_lim,
                                                            arguments.threads,
                                                            arguments.output_dir)
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run TransDecoder + DIAMOND
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start TransDecoder...')

    # 1st step
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tRun prediction of ORFs')
    transdecoder_output_dir = run_external_tools.run_transdecoder_longorfs(arguments.output_dir,
                                                                           rnaspades_output_dir)
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tORFs predicted')

    # 2nd step
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tStart searching homologous proteins by using DIAMOND')
    diamond_file_for_searching_homologous_prots = run_external_tools.run_diamond_blastp(os.path.abspath(arguments.diamond_db),
                                                                                        arguments.output_dir,
                                                                                        'diamond_output_for_searching_homologous_prots/',
                                                                                        option_parser.args_for_searching_homologous_prots,
                                                                                        f'{transdecoder_output_dir}/transcripts.fasta.transdecoder_dir/longest_orfs.pep',
                                                                                        arguments.threads,
                                                                                        'homologous_compairing.outfmt6',
                                                                                        taxon_list='33090')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tCompleted')

    # 3rd step
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tRun prediction of CDS and translated CDS')
    run_external_tools.run_transdecoder_predict(arguments.output_dir,
                                                rnaspades_output_dir,
                                                transdecoder_output_dir,
                                                diamond_file_for_searching_homologous_prots)
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' \tCDS and translated CDS predicted')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run extraction longest iso after TransDecoder
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Extract longest isoforms...')
    extraction_after_transdecoder_output_dir = arguments.output_dir + 'longest_iso_after_transdecoder/'
    os.mkdir(extraction_after_transdecoder_output_dir)
    extractor_1st = extract_longest_iso.LongestIsoExtractor(cds_fpath=f'{transdecoder_output_dir}/transcripts.fasta.transdecoder.cds',
                                                            ouf_path=extraction_after_transdecoder_output_dir,
                                                            is_first_filtration=True)
    extractor_1st.filter_cds_after_transdecoder()
    extractor_1st.filter_annot_file(f'{transdecoder_output_dir}/transcripts.fasta.transdecoder.gff3')
    extractor_1st.filter_other_fasta(f'{transdecoder_output_dir}/transcripts.fasta.transdecoder.pep')
    extractor_1st.filter_other_fasta(f'{rnaspades_output_dir}transcripts.fasta')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run CD-HIT-EST
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start clustering CDS by using CD-HIT...')
    cdhit_output_dir = run_external_tools.run_cdhit(arguments.memory_lim,
                                                    arguments.threads,
                                                    arguments.output_dir,
                                                    extraction_after_transdecoder_output_dir,
                                                    arguments.seq_idy_threshold,
                                                    arguments.alignment_type,
                                                    arguments.subseq_len_matching_cov)
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run extraction longest iso from protein fasta and annotation files by clustered cds fasta
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Extract longest iso from protein fasta and annotation files after CD-HIT...')
    extractor_2nd = extract_longest_iso.LongestIsoExtractor(cds_fpath=f'{cdhit_output_dir}clust_cds_longest_iso.fasta',
                                                            ouf_path=cdhit_output_dir,
                                                            prefix='clust_')
    extractor_2nd.get_l_of_longest_seqids_from_clust_cds()
    extractor_2nd.filter_annot_file(f'{extraction_after_transdecoder_output_dir}annotation_longest_iso.gff3')
    extractor_2nd.filter_other_fasta(f'{extraction_after_transdecoder_output_dir}proteins_longest_iso.fasta')
    extractor_2nd.filter_other_fasta(f'{extraction_after_transdecoder_output_dir}transcripts_longest_iso.fasta')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run DIAMOND
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Compare translated CDS with nr db by using DIAMOND...')
    diamond_file_for_removal_foreign_rna = run_external_tools.run_diamond_blastp(os.path.abspath(arguments.diamond_db),
                                                                                 arguments.output_dir,
                                                                                 'diamond_output_for_removal_foreign_rna/',
                                                                                 option_parser.args_for_removal_foreign_rna,
                                                                                 f'{cdhit_output_dir}/clust_proteins_longest_iso.fasta',
                                                                                 arguments.threads,
                                                                                 'diamond_compairing_res.tsv')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # get final transcriptome assembly
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Get final transcriptome assemby...')
    final_assembly_output_dir = arguments.output_dir + 'final_assembly/'
    os.mkdir(final_assembly_output_dir)
    contaminating_transcripts = extract_longest_iso.is_contaminating_transcript(os.path.abspath(arguments.diamond_taxonomic_id),
                                                                                diamond_file_for_removal_foreign_rna)
    longest_iso_cds = extract_longest_iso.make_final_cds_file(f'{cdhit_output_dir}clust_cds_longest_iso.fasta',
                                                              contaminating_transcripts,
                                                              final_assembly_output_dir)
    extractor_3rd = extract_longest_iso.LongestIsoExtractor(cds_fpath=f'{final_assembly_output_dir}final.clust_cds_longest_iso.fasta',
                                                            ouf_path=final_assembly_output_dir,
                                                            prefix='final.',
                                                            l_of_longest_seq_ids=longest_iso_cds)
    extractor_3rd.filter_annot_file(f'{cdhit_output_dir}clust_annotation_longest_iso.gff3')
    extractor_3rd.filter_other_fasta(f'{cdhit_output_dir}clust_proteins_longest_iso.fasta')
    extractor_3rd.filter_other_fasta(f'{cdhit_output_dir}clust_transcripts_longest_iso.fasta')
    extractor_3rd.get_correct_orientation(f'{final_assembly_output_dir}final.clust_transcripts_longest_iso.fasta')
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')

    # run PANNZER2
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Start GO annotation...')
    go_annotation.run_go_annotation('final.clust_proteins_longest_iso.fasta',
                                    final_assembly_output_dir,
                                    os.path.abspath(arguments.go_tree))
    print(datetime.now().strftime('%Y.%m.%d %H:%M:%S') + ' Done\n')
    print('Thank you for using Trans2express!\n')


if __name__ == '__main__':
    run_pipeline()
