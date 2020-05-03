#!/usr/bin/env python
# coding: utf-8

# # Definitions

# Annotate each gene as being:
# 1. 1:1 across all species → label as: “ortholog_one2one_allspecies”
#     6,725 genes
# 2. 1:1, 1:many, many:many across some species → label as: “ortholog_some_species”
# 3. Not orthologous in any other species → label as: “no_known_ortholog”
# 

# ## Imports

# In[26]:


import csv
import gzip
import gtfparse
import os

import pandas as pd


pd.options.display.max_colwidth = 500


# ## Outdir

# In[41]:


annotation_gtfs = {
    "human": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/gencode.v30.annotation.gtf",
    "chimpanzee": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/GCA_002880755.3_Clint_PTRv2_genomic.gtf.gz",
    "bonobo": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Pan_paniscus.panpan1.1.97.gtf",
    "gorilla": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Gorilla_gorilla.gorGor4.97.gtf.gz",
    "orangutan": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Pongo_abelii.PPYG2.97.gtf.gz",
    "macaque": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz",
    "mouse": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/gencode.vM21.annotation.gtf",
    "opossum": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Monodelphis_domestica.monDom5.97.gtf.gz",
    "platypus": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Ornithorhynchus_anatinus.OANA5.97.gtf.gz",
    "chicken": "/Users/pranathivemuri/Downloads/brawand2011/gtfs/Gallus_gallus.GRCg6a.97.gtf.gz"
 }


# In[40]:


orthology_csvs = {
    "chimpanzee": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_chimpanzee.csv",
    "bonobo": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_bonobo.csv",
    "gorilla": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_gorilla.csv",
    "orangutan": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_orangutan.csv",
    "macaque": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_macaque.csv",
    "mouse": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_mouse.csv",
    "opossum": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_opossum.csv",
    "platypus": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_platypus.csv",
    "chicken": "/Users/pranathivemuri/Downloads/brawand2011/csvs/human_chicken.csv"
 }


# In[13]:


homology_types = [
    'ortholog_one2one', 'ortholog_one2many', 'ortholog_many2many', "no_known_ortholog"]


# # Add orthology as gtf attribute

# In[27]:


def add_gtf_column(species):
    species_gene_id = "{} gene stable ID".format(species.capitalize())
    homology_type_column = "{} homology type".format(species.capitalize())
    gtf_file = annotation_gtfs[species]
    orthology_df = pd.read_csv(orthology_csvs[species])
    filtered_df = orthology_df.loc[
        orthology_df[homology_type_column].isin(homology_types)]

    gtf_df = gtfparse.read_gtf(gtf_file)
    groupby = 'orthology_type'
    gtf_df[groupby] = None
    species_gene_name = "{} gene name".format(species.capitalize())

    for orthology_type, df in filtered_df.groupby(homology_type_column):
        print("--- {} of {} ---".format(orthology_type, species))
        gene_names = set(df[species_gene_name].drop_duplicates().values)
        print(len(gene_names))
        rows = gtf_df.gene_name.isin(gene_names)
        gtf_df.loc[rows, groupby] = orthology_type
    gtf_df[groupby].value_counts(dropna=False)
    gtf_df[groupby] = gtf_df[groupby].fillna('unknown_orthology')
    gtf_df.to_csv(os.path.basename(gtf_file))


# In[31]:


for species, orthology_csv_file in orthology_csvs.items():
    if species in ["human", "chimpanzee", "mouse"]:
        continue
    add_gtf_column(species)


# In[ ]:


homology_types = [
    'ortholog_one2one', 'ortholog_one2many', 'ortholog_many2many', "no_known_ortholog"]


# In[67]:


def add_gtf_column_all_species(species):
    print("species is {}".format(species))
    gtf_file = annotation_gtfs[species]
    gtf_df = gtfparse.read_gtf(gtf_file)
    groupby = 'orthology_type'
    gtf_df[groupby] = None

    all_one_one_genes_species = orthology_table_one_one_all_species_df[species].values.tolist()
    print(len(all_one_one_genes_species))
    rows = gtf_df.gene_id.isin(all_one_one_genes_species)
    print(len(rows))
    gtf_df.loc[rows, groupby] = "ortholog_one2one_allspecies"

    genes = []
    for index, row in enumerate(rows):
        if row:
            genes.append(gtf_df.iloc[index]["gene_id"])

    assert len(set(genes)) == 6725


    filtered_df = orthology_table_with_all_species_df.loc[
        ~orthology_table_with_all_species_df["other_species_id"].isin(all_one_one_genes_species)]

    for orthology_type, df in filtered_df.groupby("orthology_type"):
        gene_ids = set(df["other_species_id"].drop_duplicates().values)
        rows = gtf_df.gene_id.isin(gene_ids)
        gtf_df.loc[rows, groupby] = "ortholog_some_species"
    gtf_df[groupby] = gtf_df[groupby].fillna('no_known_ortholog')
    print(gtf_df[groupby].value_counts(dropna=False))
    gtf_df.to_csv(os.path.basename(gtf_file))


# In[35]:


orthology_table_with_all_species_df = pd.read_csv("/Users/pranathivemuri/brawand_2011_species_orthology_types_csv/tidied_df_all_species_all_homology_types.csv")
orthology_table_one_one_all_species_df = pd.read_csv("/Users/pranathivemuri/brawand_2011_species_orthology_types_csv/all_species_by_name_one_one.csv")


# In[71]:


for species, orthology_csv_file in orthology_csvs.items():
    if species in ["human", "bonobo", "mouse", "gorilla", "macaque", "orangutan", "chimpanzee"]:
        continue
    add_gtf_column_all_species(species)


# In[ ]:




