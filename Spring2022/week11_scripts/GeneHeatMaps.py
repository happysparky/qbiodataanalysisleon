import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def intersectDataFrames(df1, df2):
    shared_samples = np.intersect1d(df1.index, df2.index)
    shared_genes = np.intersect1d(df1.columns, df2.columns)
    return df1.loc[shared_samples, shared_genes], df2.loc[shared_samples, shared_genes]


def createCorrelationMatrix(genes, df):
    ncomparisons = len(genes)
    corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                           index=genes,
                           columns=genes)
    for g1 in genes:
        for g2 in genes:
            rho, pval = stats.spearmanr(df[g1], df[g2], nan_policy="omit")
            corr_df.loc[g1, g2] = rho

    return corr_df


def createCorrelationHeatmap(path, names, df1, df2):
    fig, ax = plt.subplots(1, 2, constrained_layout=True)

    sns.heatmap(
        df1,
        cmap='mako',
        ax=ax[0]
    )

    sns.heatmap(
        df2,
        cmap='mako',
        ax=ax[1]
    )

    # some practice using enumerate
    for i, name in enumerate(names):
        ax[i].set(title=name)

    plt.savefig(path)


def main():
    genes = ["KRAS", "BRAF", "TP53", "CD3E", "TPI1", "GAPDH", "PDHA1"]

    protein_data = pd.read_csv("../analysis_data/cptac/cptac_protein.csv")
    rna_data = pd.read_csv("../analysis_data/cptac/cptac_rna.csv")

    protein_shared, rna_shared = intersectDataFrames(protein_data, rna_data)

    protein_corrs = createCorrelationMatrix(genes, protein_data)
    rna_corrs = createCorrelationMatrix(genes, rna_data)

    createCorrelationHeatmap(
        "correlation_heatmaps.png",
        ["Protein", "RNA"],
        protein_corrs,
        rna_corrs)


if __name__ == "__main__":
    main()