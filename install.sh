#!/bin/bash

# install external tools and python packages
conda env create -f trans2express_env.yml
conda activate trans2express_env

# download databases
mkdir db
wget http://gcf.fbb.msu.ru/shelkmike/Diff/Trans2express_files/Ver1/goTree.txt -P db
wget http://gcf.fbb.msu.ru/shelkmike/Diff/Trans2express_files/Ver1/nr.dmnd -P db
wget http://gcf.fbb.msu.ru/shelkmike/Diff/Trans2express_files/Ver1/taxonomic_id_to_full_taxonomy.txt -P db
