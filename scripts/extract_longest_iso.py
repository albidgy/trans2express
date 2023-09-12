import os
import pandas as pd
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser


class LongestIsoExtractor:
    def __init__(self, cds_fpath, ouf_path, is_first_filtration=False, prefix=None, l_of_longest_seq_ids=None):
        self.cds_fpath = cds_fpath
        self.ouf_path = ouf_path
        self.is_first_filtration = is_first_filtration
        self.prefix = prefix
        self.l_of_longest_seq_ids = l_of_longest_seq_ids
        self._D_OF_INIT_NAMES = {'cds': 'cds_longest_iso.fasta',
                                 'pep': 'proteins_longest_iso.fasta',
                                 'gff3': 'annotation_longest_iso.gff3',
                                 'fasta': 'transcripts_longest_iso.fasta',
                                 }
        self._d_of_cds = {}

    def filter_cds_after_transdecoder(self):
        self.l_of_longest_seq_ids = []
        output_filename = self.ouf_path + self._D_OF_INIT_NAMES[self.cds_fpath.split('.')[-1]]
        with open(self.cds_fpath) as f:
            d = dict()
            for seq_id, seq in SimpleFastaParser(f):
                d[seq_id] = seq
            sequences = pd.Series(d)
        lengths = sequences.map(len)

        index = lengths.groupby(lambda x: x.split('_')[6]).idxmax()
        longest_isoforms = sequences[index]

        with open(output_filename, 'w') as ouf:
            for seq_id, seq in longest_isoforms.items():
                print('>{}\n{}'.format(seq_id.split(' ')[0], seq), file=ouf)
                self.l_of_longest_seq_ids.append(seq_id.split(' ')[0])

    def get_l_of_longest_seqids_from_clust_cds(self):
        self.l_of_longest_seq_ids = []
        with open(self.cds_fpath) as inf:
            for line in inf:
                if line[0] == '>':
                    transcript_id = line.strip('>').strip('\n')
                    self.l_of_longest_seq_ids.append(transcript_id)

    @staticmethod
    def _parse_gff_file(line):
        line = re.sub(r'.p\d+', '', line)
        line = re.sub(r'~~.+?;', ';', line)
        line = re.sub(r'GENE.', '', line)
        return line

    def filter_annot_file(self, gff_file):
        is_longest_transcript = False
        if self.is_first_filtration:
            output_filename = self.ouf_path + self._D_OF_INIT_NAMES[gff_file.split('.')[-1]]
        else:
            output_filename = self.ouf_path + self.prefix + os.path.basename(gff_file)

        with open(gff_file) as inf:
            for line in inf:
                if line.startswith('NODE_'):
                    if line.split('\t')[2] == 'gene':
                        transcript_id = re.findall(r'~~(.+)?;Name', line)[0]
                        if transcript_id in self.l_of_longest_seq_ids:
                            is_longest_transcript = True
                            with open(output_filename, 'a') as ouf:
                                if self.prefix == 'final.':
                                    line = self._parse_gff_file(line)
                                ouf.write(line)
                        else:
                            is_longest_transcript = False
                    else:
                        if is_longest_transcript:
                            with open(output_filename, 'a') as ouf:
                                if self.prefix == 'final.':
                                    line = self._parse_gff_file(line)
                                ouf.write(line)

    def filter_other_fasta(self, other_fasta):
        is_good_transcript = False
        is_first_line = True
        if self.is_first_filtration:
            output_filename = self.ouf_path + self._D_OF_INIT_NAMES[other_fasta.split('.')[-1]]
        else:
            output_filename = self.ouf_path + self.prefix + os.path.basename(other_fasta)

        with open(other_fasta) as inf:
            for line in inf:
                transcript_id = line.strip('\n').split(' ')[0]
                if '.p' not in transcript_id and is_first_line:
                    self.l_of_longest_seq_ids = [seqid.split('.p')[0] for seqid in self.l_of_longest_seq_ids]
                is_first_line = False

                if line[0] == '>':
                    if transcript_id[1:] in self.l_of_longest_seq_ids:
                        is_good_transcript = True
                        with open(output_filename, 'a') as ouf:
                            ouf.write(transcript_id.split('.p')[0] + '\n')
                    else:
                        is_good_transcript = False
                else:
                    if is_good_transcript:
                        with open(output_filename, 'a') as ouf:
                            ouf.write(line)

    @staticmethod
    def _complement(seq):
        dna_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        complement_l = []
        for nucl in seq:
            complement_l.append(dna_dict[nucl])
        return ''.join(complement_l[::-1])

    @staticmethod
    def _split_fasta_by_n_nucl(seq):
        seq_splitted = []
        prev_idx = 0
        for idx in range(60, len(seq), 60):
            seq_splitted.append(seq[prev_idx:idx])
            prev_idx = idx
        seq_splitted.append(seq[prev_idx:len(seq)])
        return '\n'.join(seq_splitted)

    def read_cds(self):
        with open(self.cds_fpath) as cds:
            for line in cds:
                line = line.strip('\n')
                if line[0] == '>':
                    transcript_id = line
                    self._d_of_cds[transcript_id] = []
                else:
                    self._d_of_cds[transcript_id] += list(line)

        with open(self.cds_fpath, 'w') as ouf:
            for key, val in self._d_of_cds.items():
                splitted_val = self._split_fasta_by_n_nucl(''.join(val))
                ouf.write(f'{key}\n{splitted_val}\n')

    def get_correct_orientation(self, transcript_fpath):
        d_of_transcripts = {}
        self.read_cds()

        with open(transcript_fpath) as transcript:
            for line in transcript:
                line = line.strip('\n')
                if line[0] == '>':
                    transcript_id = line
                    d_of_transcripts[transcript_id] = []
                else:
                    d_of_transcripts[transcript_id] += list(line)

        with open(transcript_fpath, 'w') as ouf:
            for key, val in d_of_transcripts.items():
                transcript_seq = ''.join(val)
                cds_seq = ''.join(self._d_of_cds[key])
                if cds_seq in transcript_seq:
                    splitted_seq = self._split_fasta_by_n_nucl(transcript_seq)
                    ouf.write(f'{key}\n{splitted_seq}\n')
                elif cds_seq in self._complement(transcript_seq):
                    splitted_seq = self._split_fasta_by_n_nucl(self._complement(transcript_seq))
                    ouf.write(f'{key}\n{splitted_seq}\n')


def make_dict_by_diamond_db(path_to_taxonomic_id):
    d_of_taxid_and_clade = {}
    with open(path_to_taxonomic_id) as taxonomy:
        for line in taxonomy:
            taxid = line.split(' ')[0]
            clade = re.findall(f'clade=(.+?);', line)
            if clade:
                if taxid not in d_of_taxid_and_clade.keys():
                    d_of_taxid_and_clade[taxid] = []
                d_of_taxid_and_clade[taxid] += clade
            else:
                d_of_taxid_and_clade[taxid] = ['Other taxon']
    return d_of_taxid_and_clade


def get_compairing_taxon(path_to_diamond_res_file):
    d_of_protein_pairs = {}

    with open(path_to_diamond_res_file) as diamond:
        for line in diamond:
            line_l = line.strip().split('\t')
            transcript_id = line_l[0]
            if transcript_id not in d_of_protein_pairs.keys():
                d_of_protein_pairs[transcript_id] = []
            if len(line_l) == 3:
                taxon_id = line_l[2].split(';')[0]
                if len(d_of_protein_pairs[transcript_id]) < 5:
                    d_of_protein_pairs[transcript_id].append(taxon_id)
    return d_of_protein_pairs


def is_contaminating_transcript(path_to_taxonomic_id, path_to_diamond_res_file):
    contaminating_transcripts = []
    d_of_taxid_and_clade = make_dict_by_diamond_db(path_to_taxonomic_id)
    d_of_protein_pairs = get_compairing_taxon(path_to_diamond_res_file)

    for key, val_l in d_of_protein_pairs.items():
        counter_contaminating_rna = 0
        for taxon_id in val_l:
            if taxon_id not in d_of_taxid_and_clade.keys():
                continue
            if 'Embryophyta' not in d_of_taxid_and_clade[taxon_id]:
                counter_contaminating_rna += 1
        if len(val_l) == 5 and counter_contaminating_rna >= 3:
            contaminating_transcripts.append(key.split('.p')[0])
        elif len(val_l) > 0 and counter_contaminating_rna / len(val_l) >= 0.5:
            contaminating_transcripts.append(key.split('.p')[0])
    return contaminating_transcripts


def make_final_cds_file(cds_file, contaminating_transcripts, output_dir, prefix='final.'):
    d_of_cds_seq = {}
    l_of_good_cds = []
    file_out = output_dir + prefix + os.path.basename(cds_file)

    with open(cds_file) as inf:
        for line in inf:
            if line[0] == '>':
                transcript_id_full_name = line.strip('\n')
                d_of_cds_seq[transcript_id_full_name] = []
            else:
                d_of_cds_seq[transcript_id_full_name].append(line.strip('\n'))

    with open(file_out, 'w') as ouf:
        for key, val in d_of_cds_seq.items():
            transcript_id = key.split('.p')[0]
            if transcript_id[1:] not in contaminating_transcripts and 'N' not in val[0]:
                l_of_good_cds.append(key.strip('\n').strip('>'))
                seq_res = ''.join(val)
                ouf.write(f'{transcript_id}\n{seq_res}\n')
    return l_of_good_cds
