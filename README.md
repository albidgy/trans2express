# Trans2express: de novo hybrid transcriptome assembly tool

This pipeline allows transcriptome assembly, aiming to leave one transcript per gene. 64-bit Linux and macOS are supported.

<div style='justify-content: center'>
<img src="img/Fig1.png" align='center', width="50%">
</div>

### Installation
```bash
git clone https://github.com/albidgy/trans2express
cd trans2express/
```
install python libraries and __required__ databases:
```bash
./install.sh
```
Set pathways to the following tools on your system in the __CONFIGURATIONS.txt__ file. By default, they are all supposed to be in PATH.

### Command line options
```bash
python3 trans2express.py [options]
```

#### General options

`-1 / --short_reads1 ` [__Requied__] Forward short reads in fastq or fastq.gz format

`-2 / --short_reads2 `Reversed short reads in fastq or fastq.gz format. If you have single-end reads do not specify this parameter.

`--long_reads ` [__Requied__] Nanopore long reads in fastq or fastq.gz format.

`-o / --output_dir ` Output directory. By default, output directory is ../res_trans2express_YEAR_MONTH_DAY_HOUR_MINUTE_SECOND.

`-c/ --config_file ` Path to configurations file. By default, configuration file is CONFIGURATIONS.txt.

`-t / --threads ` Number of threads. By default, is 1.

`-m / --memory_lim ` Memory limit in Gb. By default, is 10 Gb.

`-h / --help ` See more information.

#### Optional arguments

`--diamond_db ` DIAMOND nr database (nr.dmnd) for finding homologous proteins for prediction CDS by TransDecoder and for removing foreign transcripts. The database is downloaded when you run the install.sh script or you can create your own .dmnd db. By default, is db/nr.dmnd.

`--diamond_taxonomic_id ` File with taxonomic ids list by nr.dmnd database for removal foreign rna. The file is downloaded when you run the install.sh script. By default, is db/taxonomic_id_to_full_taxonomy.txt.

`--go_tree ` File with broad GO terms. The file is downloaded when you run the install.sh script, or you can create your own goTree file. By default, is db/goTree.txt.

`--min_short_read_length ` Minimum length of short reads for fastp tool. By default, is 50.

`--seq_idy_threshold ` Sequence identity threshold for CD-HIT-EST tool. By default, is 0.98.

`--alignment_type ` Select type of alignment for CD-HIT-EST tool. 0 - local alignment, 1 - global alignment. By default, is 0.

`--subseq_len_matching_cov ` Alignment coverage for the shorter sequence for CD-HIT-EST tool. By default, is 0.6.

### Output data

As a result of the pipeline's work, the __final_assemly__ main directory is created, in which the following files are located:
 
- `final.clust_transcripts_longest_iso.fasta` - assemblied transcriptome fasta file;
- `final.clust_annotation_longest_iso.gff3` - annotation file in gff format;
- `final.clust_proteins_longest_iso.fasta` - proteins fasta file;
- `final.clust_cds_longest_iso.fasta` - CDS fasta file;
- `GO_annotation.txt` - file with GO terms.

### Citations

If you use Trans2express in your research, please cite [link](https://doi.org/10.1101/2024.01.11.575187).
