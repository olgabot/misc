import glob
import screed
import os
from multiprocessing import Pool


K2G_FASTA_PATH = "/mnt/ibm_lg/pranathi/kmer-hashing/golumbeanu2018/kmermaid/k2g_uniprot/extract_coding/"
NO_K2G_FASTA_PATH = "/mnt/ibm_lg/pranathi/kmer-hashing/golumbeanu2018/kmermaid/no_k2g_uniprot/extract_coding/"
K2G_INTERSECT = "/mnt/ibm_lg/pranathi/kmer-hashing/golumbeanu2018/kmermaid/k2g_uniprot/analysis_july9/"


if not os.path.exists(K2G_INTERSECT):
    os.makedirs(K2G_INTERSECT)
else:
    pass

k2g_fastas = glob.glob(os.path.join(K2G_FASTA_PATH, "*_peptides_clean.fasta"))
no_k2g_fastas = glob.glob(
    os.path.join(NO_K2G_FASTA_PATH, "*_peptides_clean.fasta"))
print(len(k2g_fastas))
print(len(no_k2g_fastas))


def get_record_names(path):
    print(path)
    db = screed.read_fasta_sequences(path)
    names = db.keys()
    db.close()
    return names


def save_reads_only_in_k2g(index):
    k2g = k2g_fastas[index]
    no_k2g = os.path.join(NO_K2G_FASTA_PATH, os.path.basename(k2g))
    record_names_no_k2g = get_record_names(no_k2g)
    filename = os.path.join(K2G_INTERSECT, os.path.basename(k2g).replace(
        ".fasta", "_intersect.fasta"))
    result_fasta = open(filename, "a")
    print(record_names_no_k2g[:5])
    print(k2g)
    db_k2g = screed.read_fasta_sequences(k2g)
    for name in db_k2g:
        if name not in record_names_no_k2g:
            continue
        result_fasta.write(">{}\n{}\n".format(name, db_k2g[name].sequence))


with Pool(12) as p:
    p.map(save_reads_only_in_k2g, range(len(k2g_fastas)))
