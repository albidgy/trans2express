import os
import re
import requests
import time
from bs4 import BeautifulSoup


URL = 'http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/dump_panz.cgi'
data = {'method': 'Pannzer',
        'qcov': '0.6',
        'coverageoperator': 'AND',
        'scov': '0.6',
        'minlali': '30',
        'minpide': '0.4',
        'maxpide': '1.0',
        'bestcluster_DE': 'bestcluster_DE',
        'predictor': 'ARGOT',
        'filter_output': 'filter_output',
        'b2g_thresh': '55',
        'goslim': 'None',
        'maxquery': '100000',
        }


def prepare_prot_fasta(fasta_file):
    d_of_prot_seq = {}
    input_l_of_seqs = []
    with open(fasta_file) as prot_fasta:
        for line in prot_fasta:
            if line[0] == '>':
                transcript_id = line
                d_of_prot_seq[transcript_id] = []
            else:
                if line.endswith('\n'):
                    line = line.strip('\n')
                if line.endswith('*'):
                    line = line.strip('*')
                d_of_prot_seq[transcript_id].append(line)

    for key, val in d_of_prot_seq.items():
        val_str = ''.join(val) + '\n'
        input_l_of_seqs.extend([key, val_str])

    return ''.join(input_l_of_seqs)


def send_requests(fasta_file):
    is_ready = False
    query_seqs = prepare_prot_fasta(fasta_file)
    data['query'] = query_seqs

    response_sent = requests.post(URL, data=data)
    print(response_sent.request.url)  # TODO: remove this line
    while not is_ready:
        response_get = requests.get(response_sent.request.url)
        soup = BeautifulSoup(response_get.content, 'lxml')
        if soup.find('h1', string='Job status: Finished') is not None:
            is_ready = True
        else:
            is_ready = False
            time.sleep(600)
    return soup, response_sent


def make_full_go_lists(d_of_go_terms, fpath_to_go_tree):
    d_of_go_tree = {}
    d_full_go_lists = {}
    with open(fpath_to_go_tree) as go_tree:
        for line in go_tree:
            go_l = line.strip('\n').split(' ')
            if go_l[-1] not in d_of_go_tree.keys():
                d_of_go_tree[go_l[-1]] = []
            d_of_go_tree[go_l[-1]].extend(go_l)

    for key, val in d_of_go_terms.items():
        tmp_go_l = []
        for go_term in val:
            if go_term in d_of_go_tree.keys():
                tmp_go_l += d_of_go_tree[go_term]
        d_full_go_lists[key] = list(set(tmp_go_l))
    return d_full_go_lists


def run_go_annotation(prot_fasta_file, final_assembly_output_dir, fpath_to_go_tree):
    prot_fasta_file = final_assembly_output_dir + prot_fasta_file
    soup, response_sent = send_requests(prot_fasta_file)
    d_of_go_terms = {}
    html_files = soup.find('li').find_all('a')
    for html_idx in range(len(html_files)):
        file_postfix = soup.find('li').find_all('a')[html_idx]['href']
        abspath = os.path.dirname(response_sent.request.url) + '/' + file_postfix
        cur_file = BeautifulSoup(requests.get(abspath).content, 'lxml')
        main_text = cur_file.find_all('tr')

        for idx in range(len(main_text)):
            transcript_id = re.findall(r'<td>(NODE_.+?)<br/>', str(main_text[idx]))
            if transcript_id:
                go_l = re.findall(r'(GO:.+?)<', str(main_text[idx]))
                d_of_go_terms[transcript_id[0]] = go_l
    d_full_go_lists = make_full_go_lists(d_of_go_terms, fpath_to_go_tree)

    with open(f'{final_assembly_output_dir}/GO_annotation.txt', 'w') as ouf:
        for key, val in d_full_go_lists.items():
            val_str = ', '.join(val)
            ouf.write(f'{key}\t{val_str}\n')
