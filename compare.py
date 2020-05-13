#!/usr/bin/env python
# coding: utf-8

# # Find reads that match in K2G but not in the non K2G proteome
# 1. Inputs
#     a. Predicted protein sequence from kmermaid + non K2G fasta 12, one for each sample
#     b. Predicted protein sequence from kmermaid + K2G fasta 12, one for each sample
# 2. Tasks
#     a. Use some scripting (python, bioawk, etc) to identify read IDs from the same sample that appear ONLY in the K2G fasta
# 3. Outputs
#     a. Peptide fasta file containing subset of K2G reads that appears only when using K2G peptide database 12, one for each sample
# 
import glob
import os
import screed
import pandas as pd
import os

K2G_FASTA_PATH = "/data/pranathi/kmer-hashing/golumbeanu2018/kmermaid/k2g_uniprot/extract_coding/"
NO_K2G_FASTA_PATH = "/data/olga/kmer-hashing/golumbeanu2018/kmermaid/no_k2g_uniprot/extract_coding/"
K2G_INTERSECT = "/data/pranathi/kmer-hashing/golumbeanu2018/kmermaid/k2g_uniprot/analysis_april26/"

if not os.path.exists(K2G_INTERSECT):
    os.makedirs(K2G_INTERSECT)
else:
    pass

k2g_fastas = glob.glob(os.path.join(K2G_FASTA_PATH, "*_peptides.fasta"))
no_k2g_fastas = glob.glob(os.path.join(NO_K2G_FASTA_PATH, "*_peptides.fasta"))


def save_reads_only_in_k2g(k2g, no_k2g):
    record_names_no_k2g = []
    for record in screed.open(no_k2g):
        record_names_no_k2g.append(record['name'])
    filename = os.path.join(K2G_INTERSECT, os.path.basename(k2g).replace(".fasta", "_intersect.fasta"))
    result_fasta = open(filename, "a")
    
    for record in screed.open(k2g):
        if record['name'] not in record_names_no_k2g:
            continue
        result_fasta.write(">{}\n{}\n".format(record['name'], record['sequence']))  


for k2g_fasta_peptide_sample in k2g_fastas:
    no_k2g_fasta_peptide_sample = os.path.join(NO_K2G_FASTA_PATH, os.path.basename(k2g_fasta_peptide_sample))
    save_reads_only_in_k2g(k2g_fasta_peptide_sample, no_k2g_fasta_peptide_sample)
