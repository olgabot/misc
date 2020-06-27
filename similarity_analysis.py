import glob
import os
import pandas as pd
import umap
import seaborn as sns
from os.path import expanduser
import matplotlib.pyplot as plt

# Olga's point - Can you combine both treated and untreated in the same umap,
# colored by treatment?
# they are on different scales, i can check if I can do that
# ksize 93 logsize 14, molecule protein

home = expanduser("~")

K2G_COMPARE_PATH = os.path.join(home, "k2g_results/")
NO_K2G_COMPARE_PATH = os.path.join(home, "no_k2g_results/")
K2G_CSV_COMPARE = os.path.join(home, "umap_analysis_june27/")

if not os.path.exists(K2G_CSV_COMPARE):
    os.makedirs(K2G_CSV_COMPARE)
else:
    pass

k2g_csvs = glob.glob(os.path.join(K2G_COMPARE_PATH, "*.csv"))
no_k2g_csvs = glob.glob(os.path.join(NO_K2G_COMPARE_PATH, "*csv"))
print(len(k2g_csvs))
print(len(no_k2g_csvs))


def umap_similarities(csv1, csv2):
    similarities = pd.DataFrame(pd.read_csv(csv1))
    sample_ids = similarities.columns
    similarities.columns = similarities.index = sample_ids
    embedding = umap.UMAP(
        n_neighbors=5, verbose=True).fit_transform(similarities)
    x = embedding[:, 0]
    y = embedding[:, 1]
    df_treated = pd.DataFrame.from_dict(
        {"x": x, "y": y, "treatment": "treated"})

    similarities = pd.DataFrame(pd.read_csv(csv2))
    sample_ids = similarities.columns
    similarities.columns = similarities.index = sample_ids
    embedding = umap.UMAP(
        n_neighbors=5, verbose=True).fit_transform(similarities)
    x = embedding[:, 0]
    y = embedding[:, 1]
    df_untreated = pd.DataFrame.from_dict(
        {"x": x, "y": y, "treatment": "untreated"})
    print(df_treated.head())
    print(df_untreated.head())
    concatenated = pd.concat(
        [df_treated.assign(dataset='treated'),
         df_untreated.assign(dataset='untreated')])
    print(concatenated.head())
    fig = plt.figure()
    basename = os.path.basename(csv1)
    figure_prefix = basename.replace(".csv", ".png")
    sns.scatterplot(x, y).set_title(basename)
    sns.scatterplot(
        x='x', y='y', data=concatenated, hue="treatment")
    fig.savefig(os.path.join(K2G_CSV_COMPARE, figure_prefix))


for csv1 in k2g_csvs:
    csv2 = os.path.join(NO_K2G_COMPARE_PATH, os.path.basename(csv1))
    umap_similarities(csv1, csv2)
