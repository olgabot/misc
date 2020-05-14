import glob
import os
import pandas as pd
import umap
import seaborn as sns
from os.path import expanduser
import matplotlib.pyplot as plt


home = expanduser("~")

K2G_COMPARE_PATH = os.path.join(home, "k2g_results/")
NO_K2G_COMPARE_PATH = os.path.join(home, "no_k2g_results/")
K2G_CSV_COMPARE = os.path.join(home, "umap_analysis_april26/")

if not os.path.exists(K2G_CSV_COMPARE):
    os.makedirs(K2G_CSV_COMPARE)
else:
    pass

k2g_csvs = glob.glob(os.path.join(K2G_COMPARE_PATH, "*.csv"))
no_k2g_csvs = glob.glob(os.path.join(NO_K2G_COMPARE_PATH, "*csv"))
print(len(k2g_csvs))
print(len(no_k2g_csvs))


def umap_similarities(csv, title):
    similarities = pd.DataFrame(pd.read_csv(csv))
    sample_ids = similarities.columns
    similarities.columns = similarities.index = sample_ids
    embedding = umap.UMAP(
        n_neighbors=5, verbose=True).fit_transform(similarities)
    x = embedding[:, 0]
    y = embedding[:, 1]
    fig = plt.figure()
    basename = os.path.basename(csv)
    figure_prefix = basename.replace(".csv", "_{}.png".format(title))
    sns.scatterplot(x, y).set_title(basename)
    fig.savefig(os.path.join(K2G_CSV_COMPARE, figure_prefix))


for csv1 in k2g_csvs:
    umap_similarities(csv1, "treated")

for csv1 in no_k2g_csvs:
    umap_similarities(csv1, "untreated")
