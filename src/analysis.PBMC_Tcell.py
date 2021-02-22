#!/usr/bin/env python


"""
Formerly "src/20200327.analysis.py"
"""

import os
import time
from os.path import join as pjoin
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc

from functools import lru_cache
import time

import numpy as np
import pandas as pd
from anndata import AnnData

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
from ngs_toolkit.general import enrichr, query_biomart


class Params:
    pass


@lru_cache(maxsize=None)
def query(species, ensembl_version):
    return query_biomart(
        attributes=["ensembl_gene_id", "external_gene_name"],
        species=species,
        ensembl_version=ensembl_version,
    )


def write_gene_expression_matrix(
    expr, output_file=None, file_format="h5ad", annotation=None
):
    from scipy.sparse import csr_matrix
    import anndata as an

    if file_format not in ["h5ad", "loom"]:
        raise ValueError("Output format must be one of 'h5ad' or 'loom'.")

    n = expr["umi"].isnull().sum()
    if n > 0:
        print(
            f"# {time.asctime()} - Dropping {n} entries with null UMI values."
        )
        expr = expr.dropna(subset=["umi"])

    print(f"# {time.asctime()} - Categorizing cells and genes.")
    if "plate_well" in expr.columns:
        expr = expr.assign(
            cell=expr["plate_well"] + "-" + expr[args.droplet_column]
        )
        expr["cell"] = (
            expr["plate_well"] + "-" + expr[args.droplet_column]
        ).astype("category")
    else:
        expr["cell"] = expr["r1"].astype("category")
    expr["gene"] = expr["gene"].astype("category")

    print(f"# {time.asctime()} - Creating sparse matrix.")
    sparse_matrix = csr_matrix(
        (expr["umi"].values, (expr["cell"].cat.codes, expr["gene"].cat.codes)),
        shape=(
            len(expr["cell"].cat.categories),
            len(expr["gene"].cat.categories),
        ),
        dtype=np.int,
    )

    print(f"# {time.asctime()} - Creating AnnData object.")
    a = an.AnnData(X=sparse_matrix)
    print(f"# {time.asctime()} - Annotating with cells and genes.")
    a.obs.index = expr["cell"].cat.categories
    a.var.index = expr["gene"].cat.categories
    a.obs["plate_well"] = a.obs.index.str.slice(0, 3)
    if annotation is not None:
        print(f"# {time.asctime()} - Adding additional annotation.")
        a.obs = (
            a.obs.reset_index()
            .set_index("plate_well")
            .join(annotation.set_index("plate_well"))
            .set_index("index")
        )

    print(f"# {time.asctime()} - Writing h5ad object to disk.")
    if output_file is None:
        output_file = args.output_prefix + file_format
    if file_format == "h5ad":
        a.write(output_file)
    if file_format == "loom":
        a.write_loom(output_file)
    return a


def differential_expression(
    an, attribute, groups="all", method="t-test_overestim_var", n_genes=250
):
    sc.tl.rank_genes_groups(
        an,
        attribute,
        groups=groups,
        method=method,
        n_genes=n_genes,
        use_raw=False,
    )
    result = an.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    res = list()
    for group in groups:
        mul = -1 if len(groups) == 2 and group == groups[-1] else 1
        res.append(
            pd.DataFrame(
                [
                    result["names"][group],
                    result["pvals"][group],
                    result["pvals_adj"][group],
                    result["logfoldchanges"][group],
                    result["scores"][group],
                ],
                index=[
                    "gene_name",
                    "pvals",
                    "pvals_adj",
                    "logfoldchanges",
                    "scores",
                ],
            ).T.assign(group=group)
        )
    diff = pd.concat(res).set_index("gene_name")
    for col in ["pvals", "pvals_adj", "logfoldchanges", "scores"]:
        diff[col] = diff[col].astype(float)
    diff = diff.sort_values(["group", "scores"])
    try:
        m = an.X.mean(0).A1
    except AttributeError:
        m = an.X.mean(0)
    mean = np.log(1 + pd.Series(m, index=an.var.index, name="mean"))
    return diff.join(an.var).join(mean)


def plot_differential_expression(diff):
    groups = diff["group"].unique()

    fig, axes = plt.subplots(1, len(groups), figsize=(4 * len(groups), 4 * 1))
    for i, group in enumerate(groups):
        d = diff.query(f"group == '{group}'")
        axes[i].set_title(group)
        axes[i].scatter(
            d["mean"],
            d["logfoldchanges"],
            c=-np.log10(d["pvals"]),
            s=2,
            alpha=0.2,
            rasterized=True,
        )
        v = d["logfoldchanges"].abs().max()
        axes[i].set_ylim((-v, v))
    return fig


def differential_enrichment(
    differential,
    groups,
    attribute,
    sort_by="scores",
    alpha=0.05,
    alpha_col="pvals",
    max_n=250,
):
    enr = list()
    gene_set_libraries = (
        ["ARCHS4_Cell-lines"]
        if "cell_line" in attributes
        else [
            "WikiPathways_2019_Human",
            "BioCarta_2016",
            "KEGG_2019_Human",
            "GO_Biological_Process_2018",
            "NCI-Nature_2016",
            "Human_Gene_Atlas",
            "GTEx_Tissue_Sample_Gene_Expression_Profiles_up",
            "GTEx_Tissue_Sample_Gene_Expression_Profiles_down",
        ]
    )
    for group in groups:
        diff_genes = (
            differential.query(f"group == '{group}' and {alpha_col} < {alpha}")
            .sort_values(sort_by)
            .tail(max_n)["external_gene_name"]
        )
        print(diff_genes)
        if diff_genes.empty:
            continue
        try:
            enr.append(
                enrichr(
                    diff_genes.to_frame(name="gene_name"),
                    gene_set_libraries=gene_set_libraries,
                ).assign(group=group)
            )
        except AttributeError:
            pass
    return pd.concat(enr)


def plot_differential_enrichment(
    enrichments, output_prefix, gene_set_libraries=None, ntop_terms=5
):
    if gene_set_libraries is None:
        gene_set_libraries = enrichments["gene_set_library"].unique()

    for gene_set_library in gene_set_libraries:
        for metric in ["combined_score", "p_value"]:
            combs = enrichments.query(
                f"gene_set_library == '{gene_set_library}'"
            ).pivot_table(
                index="description",
                columns="group",
                values=metric,
                fill_value=0 if metric == "combined_score" else 1,
            )
            combs.columns = combs.columns.astype(str).str.replace("_n", "")

            if metric == "p_value":
                combs = -np.log10(combs)

            # Top X terms per group, together
            for n in [1, 2, 5] if "cell_line" in attributes else [ntop_terms]:
                m = combs.reset_index().melt(id_vars="description")
                terms = m.loc[
                    m.groupby("group")["value"]
                    .nlargest(n)
                    .index.get_level_values(1)
                    .tolist(),
                    "description",
                ]
                combs2 = combs.reindex(terms.drop_duplicates())

                vmax = combs.max().max()
                vmax += vmax * 0.15
                grid = sns.clustermap(
                    combs2,
                    metric="correlation",
                    cbar_kws={"label": f"{metric}"},
                    square=True,
                    vmax=vmax,
                    col_cluster=False,
                    row_cluster=False if n < 10 else True,
                    yticklabels=True,
                )
                grid.ax_heatmap.set_xlabel("scifi-RNA-seq")
                grid.ax_heatmap.set_ylabel(f"top {n} terms entries per group")
                grid.savefig(
                    output_prefix
                    + f".{attribute}.differential_enrichment.top{str(n).zfill(2)}.{metric}.{gene_set_library}.clustermap.svg",
                    dpi=300,
                    bbox_inches="tight",
                )
                grid = sns.clustermap(
                    combs2,
                    z_score=1,
                    center=0,
                    cmap="RdBu_r",
                    metric="correlation",
                    cbar_kws={"label": f"{metric} (Z-score)"},
                    square=True,
                    col_cluster=False,
                    row_cluster=False if n < 10 else True,
                    yticklabels=True,
                )
                grid.ax_heatmap.set_xlabel("scifi-RNA-seq")
                grid.ax_heatmap.set_ylabel(f"top {n} terms entries per group")
                grid.savefig(
                    output_prefix
                    + f".{attribute}.differential_enrichment.top{str(n).zfill(2)}.{metric}.{gene_set_library}.clustermap.zscore.svg",
                    dpi=300,
                    bbox_inches="tight",
                )


def add_sex_ratio(a, plot=True):
    try:
        m = pd.read_csv("human.chromosomes.csv")
    except FileNotFoundError:
        m = query_biomart(
            ["chromosome_name", "external_gene_name"], ensembl_version="grch38"
        )
        m.to_csv("human.chromosomes.csv")

    x = m.query('chromosome_name == "X"')["external_gene_name"].unique()
    xx = a.var.index[a.var["external_gene_name"].isin(x)].tolist()
    y = m.query('chromosome_name == "Y"')["external_gene_name"].unique()
    yy = a.var.index[a.var["external_gene_name"].isin(y)].tolist()

    df2 = pd.DataFrame(a.X, columns=a.var.index, index=a.obs.index).T
    df2 += abs(df2.values.min())
    df2 -= df2.min()

    sex_ratio = (df2.reindex(yy).mean() - df2.reindex(xx).mean()).dropna() * (
        df2.reindex(yy).values.mean() - df2.reindex(xx).values.mean()
    )
    sex_ratio *= 1e4
    sex_ratio = (sex_ratio - sex_ratio.mean()) / sex_ratio.std()
    a.obs = a.obs.join(sex_ratio.to_frame(name="sex_ratio"))
    if plot:
        total = df2.mean(axis=0)
        # sex = df2.reindex(yy + xx).mean()

        fig, ax = plt.subplots(1, 1, figsize=(3 * 1, 3 * 1))
        ax.axhline(0, linestyle="--", color="grey")
        ax.set_xlabel("Mean expression")
        ax.set_ylabel("Sex ratio estimate")
        ax.scatter(
            total,
            sex_ratio,
            s=2,
            alpha=0.25,
            cmap="coolwarm",
            c=sex_ratio,
            vmin=-2,
            vmax=2,
            rasterized=True,
        )
        # axes[1].scatter(total, sex_ratio / sex, s=2, alpha=0.25, cmap="coolwarm", c=sex_ratio / sex, vmin=-2, vmax=2)
        # axes[2].scatter(total, sex_ratio / sex / total, s=2, alpha=0.25, cmap="coolwarm", c=sex_ratio / sex / total, vmin=-2, vmax=2)
        # [a.obs['sex_ratio'].abs() >= 0.5]
        return a, fig
    return a


class Args:
    pass


args = Args()
args.droplet_column = "r2"


base_url = "https://biomedical-sequencing.at/projects/BSA_0383_PD_SCIFI_72203e85f90a4b0d8ad6cee221d59e6d"
samples1 = [
    "PD2XX1_10xscRNA_Human_PBMCs_2S3Qmixed",
    "PD2XX1_10xscRNA_Human_Tcells_2S3Qmixed_Unstimulated",
    "PD2XX1_10xscRNA_Human_Tcells_2S3Qmixed_Stimulated",
    "PD2XX1_10xscRNA_Mouse_LI",
]
samples2 = [
    "PD206_1_scifi_Human_PBMCs_Nuclei1pFA",
    "PD206_4_scifi_Human_Tcells_Nuclei1pFA",
    "PD205_1_scifi_Mouse_LI_CellsMeOH",
    "PD205_2_scifi_Mouse_LI_NucleiFresh",
]

# metrics = dict()
# for sample in samples1:
#     metrics_url = "/".join([base_url, "data", sample, sample]) + ".metrics.csv.gz"
#     m = pd.read_csv(metrics_url).sort_values("umi", ascending=False, ignore_index=True)
#     m.index += 1
#     metrics[sample] = m
# for sample in samples2:
#     metrics_url = "/".join([base_url, "data", sample, sample]) + ".metrics.csv.gz"
#     m = pd.read_csv(metrics_url)
#     if m.shape[1] == 10:
#         m.columns = ['read', 'unique_umi', 'umi', 'gene', 'unique_fraction', '_remove', 'plate_well', 'donor_id', 'donor_sex', 'tcr_stimulation']
#     elif m.shape[1] == 7:
#         m.columns = ['read', 'unique_umi', 'umi', 'gene', 'unique_fraction', '_remove', 'plate_well']
#     m = m.drop('_remove', axis=1).sort_values("umi", ascending=False, ignore_index=True)
#     m.index += 1
#     metrics[sample] = m

# samples = samples1 + samples2

# fig, axis = plt.subplots(2, len(samples), figsize=(6 * len(samples), 4 * 2))
# for i, sample in enumerate(samples):
#     m = metrics[sample]
#     axis[0, i].set_title(sample)
#     axis[0, i].plot(m.index, m['umi'])
#     axis[1, i].scatter(m['read'], m['unique_fraction'], s=1, alpha=0.05, rasterized=True)
#     axis[0, i].loglog()
#     axis[1, i].set_xscale("log")

# axis[0, 0].set_ylabel("Barcodes")
# axis[1, 0].set_ylabel("Unique fraction")
# axis[0, int(len(samples) / 2)].set_xlabel("Barcodes")
# axis[1, int(len(samples) / 2)].set_xlabel("Reads per barcode")
# fig.savefig("20200327.first_qc.svg", dpi=300, bbox_inches="tight")


# # Expression


# # for sample in samples1 + samples2:
# for sample in samples1:
#     if os.path.exists(sample + '.expression.h5ad'):
#         continue
#     else:
#         print(sample)
#     args.output_prefix = sample + ".expression."
#     expression_file = pjoin("~/Downloads", sample + ".expression.csv.gz")
#     metrics_file = pjoin("~/Downloads", sample + ".metrics.csv.gz")
#     try:
#         annotation = (
#             pd.read_csv(pjoin("metadata", "sciRNA-seq." + sample + '.csv'))
#             .drop(['sample_name', 'combinatorial_barcode'], axis=1))
#         attributes = annotation.drop(['plate_well'], axis=1).columns.tolist()
#     except FileNotFoundError:
#         annotation = None
#         attributes = []
#     output_prefix = pjoin(run, sample + ".")
#     expr = pd.read_csv(expression_file)
#     a = write_gene_expression_matrix(expr, annotation=annotation)


#

parameters = json.load(
    open(pjoin("metadata", "PD205.analysis_thresholds.json"), "r")
)

#

run = "results/202004_results"
os.makedirs(run, exist_ok=True)


for sample in samples1 + samples2:
    print(sample)
    params = Params()
    # set parameters as variables
    for k, v in parameters[sample].items():
        print(k, v)
        setattr(params, k, v)

    output_dir = pjoin(run, sample)  # + ".raw.")
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = pjoin(output_dir, sample)
    args.output_prefix = sample + ".expression."
    try:
        annotation = pd.read_csv(
            pjoin("metadata", "sciRNA-seq." + sample + ".csv")
        ).drop(["sample_name", "combinatorial_barcode"], axis=1)
        attributes = annotation.drop(["plate_well"], axis=1).columns.tolist()
    except FileNotFoundError:
        annotation = None
        attributes = []
    # expression_file = pjoin("~/Downloads", sample + ".expression.csv.gz")
    # metrics_file = pjoin("~/Downloads", sample + ".metrics.csv.gz")

    # r1_barcodes = pd.read_csv("metadata/737K-april-2014_rc.txt", header=None, squeeze=True)
    # r1_barcodes = pd.read_csv("metadata/737K-august-2016.txt", header=None, squeeze=True)
    r1_barcodes = pd.read_csv(
        "metadata/4M-with-alts-february-2016.txt", header=None, squeeze=True
    )
    r2_barcodes = pd.read_csv(
        "metadata/737K-cratac-v1.txt", header=None, squeeze=True
    )

    a = sc.read(pjoin(run, sample + ".expression.h5ad"))
    # This filtering here is just to get rid of the low end of barcodes
    sc.pp.filter_cells(a, min_counts=50)

    if "10xscRNA" not in sample:
        a.obs["r1"] = a.obs["plate_well"] = a.obs.index.str.slice(0, 3)
        a.obs["r2"] = a.obs.index.str.slice(4)
        a = a[a.obs["r2"].isin(r2_barcodes)]
    else:
        pass
        # a.obs['r1'] = a.obs.index
        # a = a[a.obs['r1'].isin(r1_barcodes)]
    a = a.copy()
    # a.obs = a.obs.reset_index().set_index("plate_well").join(
    #     annotation.set_index("plate_well")).reset_index().set_index("index")

    if "Mouse" in sample:
        try:
            biomart = pd.read_csv("mouse.biomart.csv", index_col=0)
        except FileNotFoundError:
            biomart = query(
                species="mmusculus", ensembl_version="grch38"
            ).set_index("ensembl_gene_id")
            biomart.to_csv("mouse.biomart.csv")
    else:
        try:
            biomart = pd.read_csv("human.biomart.csv", index_col=0)
        except FileNotFoundError:
            biomart = query(
                species="hsapiens", ensembl_version="grch38"
            ).set_index("ensembl_gene_id")
            biomart.to_csv("human.biomart.csv")
    for col in biomart.columns:
        biomart[col] = biomart[col].replace("nan", np.nan)
    a.var = a.var.join(biomart)
    # replace genes with no gene symbol with original Ensembl ID
    null = a.var.loc[a.var["external_gene_name"].isnull(), "external_gene_name"]
    null.update(null.index.to_series())
    a.var.loc[a.var["external_gene_name"].isnull(), "external_gene_name"] = null

    # QC
    # sc.pl.highest_expr_genes(a, n_top=20)
    a.var.loc[:, "n_counts"] = a.X.sum(axis=0).A1
    a.var["log_counts"] = np.log10(a.var["n_counts"])

    a.var.loc[:, "mito"] = a.var["external_gene_name"].str.contains(
        r"^MT-", case=False
    )
    a.obs.loc[:, "percent_mito"] = (
        np.sum(a[:, a.var["mito"]].X, axis=1).A1 / np.sum(a.X, axis=1).A1
    ) * 100
    a.var.loc[:, "ribo"] = a.var["external_gene_name"].str.contains(
        r"^RP", case=False
    )
    a.obs.loc[:, "percent_ribo"] = (
        np.sum(a[:, a.var["ribo"]].X, axis=1).A1 / np.sum(a.X, axis=1).A1
    ) * 100
    a.obs.loc[:, "n_counts"] = a.X.sum(axis=1).A1
    a.obs.loc[:, "log_counts"] = np.log10(a.obs.loc[:, "n_counts"])
    a.obs.loc[:, "n_genes"] = (a.X != 0).sum(1).A1
    a.obs.loc[:, "log_genes"] = np.log10(a.obs.loc[:, "n_genes"])

    a.obs.loc[:, "efficiency_ratio"] = a.obs["log_genes"] / a.obs["log_counts"]

    a.obs.loc[:, "malat1"] = (
        a.X[
            :,
            a.var["external_gene_name"]
            .str.contains("MALAT1", case=False)
            .values,
        ]
        .sum(axis=1)
        .A1
    )
    a.obs.loc[:, "percent_malat1"] = (a.obs["malat1"] / a.obs["n_counts"]) * 100

    tech_attributes = [
        "log_counts",
        "log_genes",
        "efficiency_ratio",
        "percent_mito",
        "percent_ribo",
        "percent_malat1",
    ]

    thresholds = [
        np.log10(params.min_umis_per_cell),
        np.log10(params.min_genes_per_cell),
        (params.min_efficiency_ratio, params.max_efficiency_ratio),
        (params.min_percent_mito, params.max_percent_mito),
        (params.min_percent_ribo, params.max_percent_ribo),
        (params.min_percent_malat1, params.max_percent_malat1),
    ]

    grid = sns.pairplot(
        a.obs[tech_attributes],
        plot_kws=dict(s=5, alpha=0.75, linewidths=0, rasterized=True),
    )
    grid.fig.suptitle(f"{sample}\nPre-filtering")
    for f, g in [("axhline", grid.axes), ("axvline", grid.axes.T)]:
        for ax, thresh in zip(g, thresholds):
            if thresh is None:
                continue
            for axx in ax:
                if isinstance(thresh, tuple):
                    for t in thresh:
                        getattr(axx, f)(t)
                else:
                    getattr(axx, f)(thresh)
    grid.fig.savefig(
        output_prefix + "quality_metrics.pre_filtering.svg",
        bbox_inches="tight",
        dpi=300,
    )

    a2 = a.copy()
    # # filter cells
    a = a[a.obs["n_counts"] >= params.min_umis_per_cell]
    a = a[a.obs["n_genes"] >= params.min_genes_per_cell]
    a = a[a.obs["efficiency_ratio"] >= params.min_efficiency_ratio]
    a = a[a.obs["efficiency_ratio"] <= params.max_efficiency_ratio]
    a = a[a.obs["percent_mito"] <= params.max_percent_mito]
    a = a[a.obs["percent_ribo"] >= params.min_percent_ribo]
    a = a[a.obs["percent_ribo"] <= params.max_percent_ribo]
    a = a[a.obs["percent_malat1"] >= params.min_percent_malat1]
    a = a[a.obs["percent_malat1"] <= params.max_percent_malat1]

    a.raw = a
    # # filter genes
    a = a[:, a.X.sum(0).A1 > 0]
    # plot_qc(a, suffix="pre_filtering")

    sc.pp.filter_genes(a, min_counts=a.var["n_counts"].quantile(0.1) + 1)
    a = a[:, a.X.sum(0).A1 > 0]

    # Get completely balanced design in terms of cell number
    # n = int(a.obs.groupby(attributes).count().min().min() / 3)
    # c = a.obs.groupby(attributes)['n_counts'].nlargest(n).index.get_level_values(-1).tolist()
    # a = a[c, :]

    # sns.distplot(np.log10(a.var['n_counts']), kde=False)

    sc.pp.filter_genes(a, min_cells=params.min_cells_per_gene)
    a = a[:, a.X.sum(0).A1 > 0]

    # plot_qc(a, suffix="post_filtering")

    grid = sns.pairplot(
        a.obs[tech_attributes],
        plot_kws=dict(s=5, alpha=0.75, linewidths=0, rasterized=True),
    )
    grid.fig.suptitle(f"{sample}\nPost-filtering")
    for f, g in [("axhline", grid.axes), ("axvline", grid.axes.T)]:
        for ax, thresh in zip(g, thresholds):
            if thresh is None:
                continue
            for axx in ax:
                if isinstance(thresh, tuple):
                    for t in thresh:
                        getattr(axx, f)(t)
                else:
                    getattr(axx, f)(thresh)
    grid.fig.savefig(
        output_prefix + "quality_metrics.post_filtering.svg",
        bbox_inches="tight",
        dpi=300,
    )

    sc.pp.normalize_per_cell(a, counts_per_cell_after=1e4)
    sc.pp.log1p(a)

    sc.pp.scale(a, max_value=10)
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)

    sc.tl.leiden(a, resolution=params.leiden["resolution"])
    if "leiden" not in attributes:
        attributes += ["leiden"]

    # Add sex
    if "Mouse" not in sample:
        attributes += ["sex_ratio"]
        a, fig = add_sex_ratio(a, plot=True)
        fig.savefig(
            output_prefix + ".single_cell.sex_ratio_estimate_over_mean.svg",
            dpi=300,
            bbox_inches="tight",
        )

        fig, ax = plt.subplots(1, 1, figsize=(3, 3))
        groupby = "donor_sex" if "donor_sex" in a.obs.columns else "leiden"
        sc.pl.violin(a, groupby=groupby, keys="sex_ratio", ax=ax, show=False)
        ax.axhline(0, linestyle="--", color="grey")
        fig.savefig(
            output_prefix + ".single_cell.sex_ratio.svg",
            dpi=300,
            bbox_inches="tight",
        )

    # Make sure plotting order is random
    a = a[
        np.random.choice(a.obs.index.tolist(), a.obs.shape[0], replace=False), :
    ]

    if not os.path.exists(output_prefix + ".filtered.h5ad"):
        sc.write(output_prefix + ".filtered.h5ad", a)

    a = sc.read(output_prefix + ".filtered.h5ad")

    sc.pl.pca_variance_ratio(a, log=True, show=False)
    plt.gca().figure.savefig(
        output_prefix + ".single_cell.pca_variance_ratio.svg",
        dpi=300,
        bbox_inches="tight",
    )

    fig = sc.pl.pca(
        a,
        color=tech_attributes + attributes,
        components=["1,2", "2,3", "3,4", "4,5"],
        return_fig=True,
    )
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix + ".single_cell.pca.svg", dpi=300, bbox_inches="tight"
    )

    fig = sc.pl.umap(a, color=tech_attributes + attributes, return_fig=True)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix + ".single_cell.umap.svg", dpi=300, bbox_inches="tight"
    )

    if "Mouse" not in sample:
        # Marker genes (tailored for a PBMC sample)
        mark1 = [
            "MALAT1",  # sc
            "CD34",  # HSC
            "CD3D",
            "CD3G",
            "CD247",  # T-cell
            "CD4",
            "FOXP3",
            "CCR7",
            "CTLA4",  # CD4
            "CD8A",
            "NKG7",
            "IL2RA",  # CD8
            "NCAM1",
            "GNLY",  # NK
            "CD14",
            "CST3",  # Monocytes
            "C1QA",
            "FOLR2",  # Myeloid
            "CD79A",
            "CD19",
            "IGHG1",  # B cells  (also MS4A1)
            "FCER1G",
            "CLEC10A",  # dendritic cells
            "GZMB",
            "CD68",  # monocyte/macrophage?
            "SELE",
            "CD93",
            "VWF",
            "CDH5",
            "PECAM1",
            "KDR",  # endothelial cells
            # "DCN", "GSN", "COL6A2", "PDGFRA",  # fibroblasts
            # "MYH11", "TAGLN", "ACTA2",  # smooth muscle
            # "GPAM", "LEP",  # adypocytes
            # "PLP1", "NRXN1", "NRXN3",  # Neuronal
            "CD27",
            "MS4A1",
            "CD24",
            "NCR1",
            "CD274",
        ]  # other immune
        mark2 = [
            "IL32",
            "IFNG",
            "IFNGR1",
            "IL4R",
            "IL4",
            "JUN",
            "JUNB",
            "JUND",
            "JAK1",
            "JAK2",
            "GATA1",
            "JARID2",
            "KRAS",
            "MYC",
        ]
        mark3 = ["BTK", "LCK", "E2F4", "CXCR4", "ITGA4", "HBA1", "PTPRC"]
        mark4 = [
            # T cell receptors
            "IL1B",
            "SELL",
            "CCR2",
            "GRN",
            "HIF1A",
            "HMOX1",
            "TIMP1",
            "AIF1",
            "CEBPB",
            "CXCL2",
            "CXCL8",
            "IL10",
            "FCER1A",  # anti-inflamatory
            "ADAM8",
            "NLRP3",
            "CCL2",
            "CCL3",
            "CCL4",
            "CCL5",
            "LCK",
            "ZAP70",
            "IGHA1",
            "IGHM",
            "IL2RB",
            "IL2RG",
            "IFNG",  # pro-inflamatory
        ]
        red_mark = ["CST3", "TCL1A", "GZMB", "NKG7", "CD3D"]
        marker_genes = mark1 + mark2 + mark3
    else:
        # Marker genes (tailored for mouse intestine)
        mark1 = [
            "Pecam1",
            "Cdh5",
            "Vcam1",
            "Madcam1",  # endothelium  Cdh5 == CD144
            "Epcam",
            "Thy1",
            "Ltbr",  # epithelium
            "Pdpn",  # fibroblast
        ]
        mark2 = [
            "Cd3e",
            "Cd3d",
            "Cd3g",
            "Cd4",
            "Cd8a",
            "Cd79a",
            "Cd19",
            "Cd14",
            "Cst3",  # Monocytes
        ]
        marker_genes = mark1 + mark2

    if params.plot_raw:
        g = [
            y
            for x in mark1
            for y in a.raw.var.loc[
                a.raw.var["external_gene_name"] == x
            ].index.tolist()
        ]
    else:
        g = [x for x in mark1 if x in a.var["external_gene_name"].tolist()]
    color = tech_attributes + attributes + g
    kwargs = dict(
        hspace=0.1,
        wspace=0,
        return_fig=True,
        use_raw=params.plot_raw,
        gene_symbols="external_gene_name" if not params.plot_raw else None,
    )

    fig = sc.pl.pca(a, color=color, **kwargs)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    if params.plot_raw:
        for ax in fig.axes:
            try:
                ax.set_title(
                    a.raw.var.loc[ax.get_title(), "external_gene_name"]
                )
            except KeyError:
                pass
    fig.savefig(
        output_prefix + ".single_cell.markers_red.pca.svg",
        dpi=300,
        bbox_inches="tight",
    )

    fig = sc.pl.umap(a, color=color, **kwargs)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    if params.plot_raw:
        for ax in fig.axes:
            try:
                ax.set_title(
                    a.raw.var.loc[ax.get_title(), "external_gene_name"]
                )
            except KeyError:
                pass
    fig.savefig(
        output_prefix + ".single_cell.markers_red.umap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    if "Mouse" not in sample:
        if params.plot_raw:
            g = [
                y
                for x in marker_genes
                for y in a.raw.var.loc[
                    a.raw.var["external_gene_name"] == x
                ].index.tolist()
            ]
        else:
            g = [
                x
                for x in marker_genes
                if x in a.var["external_gene_name"].tolist()
            ]
        color = tech_attributes + attributes + g
        kwargs = dict(
            hspace=0.1,
            wspace=0,
            return_fig=True,
            use_raw=params.plot_raw,
            gene_symbols="external_gene_name" if not params.plot_raw else None,
        )

        fig = sc.pl.umap(a, color=color, **kwargs)
        for ax in fig.axes:
            ax.get_children()[0].set_rasterized(True)
        if params.plot_raw:
            for ax in fig.axes:
                try:
                    ax.set_title(
                        a.raw.var.loc[ax.get_title(), "external_gene_name"]
                    )
                except KeyError:
                    pass
        fig.savefig(
            output_prefix + ".single_cell.markers_extended.umap.svg",
            dpi=300,
            bbox_inches="tight",
        )

    else:
        for label, plot_raw, markers in [
            ("structural_cells", True, mark1),
            ("immune_cells", True, mark2),
        ]:
            if plot_raw:
                g = [
                    y
                    for x in markers
                    for y in a.raw.var.loc[
                        a.raw.var["external_gene_name"] == x
                    ].index.tolist()
                ]
            else:
                g = [
                    x
                    for x in markers
                    if x in a.var["external_gene_name"].tolist()
                ]
            color = tech_attributes + attributes + g
            kwargs = dict(
                hspace=0.1,
                wspace=0.05,
                return_fig=True,
                use_raw=plot_raw,
                gene_symbols="external_gene_name" if not plot_raw else None,
            )

            fig = sc.pl.umap(a, color=color, **kwargs)
            for ax in fig.axes:
                ax.get_children()[0].set_rasterized(True)
            if plot_raw:
                for ax in fig.axes:
                    try:
                        ax.set_title(
                            a.raw.var.loc[ax.get_title(), "external_gene_name"]
                        )
                    except KeyError:
                        pass
            fig.savefig(
                output_prefix + f".single_cell.markers.{label}.umap.svg",
                dpi=300,
                bbox_inches="tight",
            )

    # Differential expression

    a2 = AnnData(a.raw.X)
    a2.obs = a.obs
    a2.var = a.raw.var
    sc.pp.filter_genes(a2, min_cells=10)
    # a2 = a2[:, a2.var.index.isin(a.var.index)]
    a2 = a2[:, a2.X.sum(0).A1 > 0]
    a2 = a2[~a2.obs["leiden"].isin(["-1", -1]), :]
    a2.raw = a2
    sc.pp.normalize_per_cell(a2, counts_per_cell_after=1e4)
    sc.pp.log1p(a2)

    for attribute in attributes:
        if not isinstance(a.obs[attribute].dtype, pd.CategoricalDtype):
            continue
        # # differential expression
        diff = differential_expression(
            a2, attribute, n_genes=params.diff_exp["max_genes"]
        )
        diff.to_csv(
            output_prefix + f".{attribute}.cluster_comparison.top_values.csv",
            index=True,
        )

        diff = pd.read_csv(
            output_prefix + f".{attribute}.cluster_comparison.top_values.csv",
            index_col=0,
        )

        fig = plot_differential_expression(diff)
        fig.savefig(
            output_prefix + f".{attribute}.differential_expression.ma_plot.svg",
            dpi=300,
            bbox_inches="tight",
        )

        # # differential enrichment
        groups = [x for x in diff["group"].unique() if x not in ["-1", -1]]

        enrichments = differential_enrichment(
            diff,
            groups,
            attribute,
            alpha=0.05,
            alpha_col="pvals_adj",
            max_n=params.diff_exp["max_genes"],
            sort_by="scores",
        )
        enrichments.to_csv(
            output_prefix + f"{attribute}.differential_enrichment.csv",
            index=False,
        )
        enrichments = pd.read_csv(
            output_prefix + f"{attribute}.differential_enrichment.csv",
            index_col=0,
        )

        g = (
            enrichments.set_index("description")
            .groupby(["group"])["combined_score"]
            .nlargest(5)
        )
        print("combined_score:\n", g)
        g = (
            enrichments.set_index("description")
            .groupby(["group"])["p_value"]
            .nsmallest(5)
        )
        print("p_value:\n", g)

        plot_differential_enrichment(enrichments, output_prefix)

        # # # display gene signatures heatmap with marker genes
        # sc.settings.figdir = os.path.dirname(output_prefix)
        # a2.obs['group'] = a.obs[attributes[::-1]].apply(lambda x: ", ".join(x), axis=1)
        # for n in [10, 20, 50, 75, 200, 350, 500, 1000]:
        #     # get genes
        #     g = diff.query("pvals < 1e-5 & abs(logfoldchanges) > 2").groupby('group')['scores'].nlargest(n).index.get_level_values(1)
        #     # replace with gene symbol
        #     g = a2.var.reindex(g)['external_gene_name']
        #     # g = a2.var.loc[a2.var['external_gene_name'].isin(g.values.tolist())].index.values.tolist()
        #     for kwargs, label in [
        #             ({}, "."),
        #             ({"standard_scale": "var", "vmax": 0.65}, ".zscore_gene.")
        #     ]:
        #         for cmap in [
        #                 None,
        #                 "inferno",
        #                 "ocean_r"
        #         ]:
        #             kwargs.update({"cmap": cmap})
        #             sc.pl.heatmap(
        #                 a2, groupby="group", var_names=g,
        #                 use_raw=False,
        #                 log=False, show=False,
        #                 gene_symbols="external_gene_name",
        #                 show_gene_labels=None,
        #                 save="." + os.path.basename(output_prefix) +
        #                      f"{attribute}.differential_test.top{str(n).zfill(2)}_genes.clustermap{label}{cmap}.svg",
        #                 **kwargs)


def stats():
    for sample in samples1 + samples2:
        output_dir = pjoin(run, sample)  # + ".raw.")
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = pjoin(output_dir, sample)
        a = sc.read(output_prefix + ".filtered.h5ad")
        print(
            sample,
            a.shape,
            a.obs["n_counts"].mean(),
            a.obs["n_genes"].mean(),
            a.obs["efficiency_ratio"].mean(),
        )


def TCR_activation_10X():
    ans = list()
    for sample, activation in zip(
        samples1[1:3], ["unstimulated", "stimulated"]
    ):
        print(sample)
        output_prefix = pjoin("results", run, sample, sample)
        b = sc.read(output_prefix + ".filtered.h5ad")
        b.obs = b.obs.assign(activation=activation)
        ans.append(b)

    output_prefix = (
        run + "/PD2XX1_10xscRNA_Human_Tcells_2S3Qmixed.samples_joined"
    )
    a = AnnData.concatenate(*ans)
    assert (
        np.asarray(a.var["external_gene_name-0"])
        == np.asarray(a.var["external_gene_name-1"])
    ).all()
    a.var["external_gene_name"] = a.var["external_gene_name-0"]
    a = a[
        np.random.choice(a.obs.index.tolist(), a.obs.shape[0], replace=False), :
    ]

    tech_attributes = [
        "log_counts",
        "log_genes",
        "efficiency_ratio",
        "percent_mito",
        "percent_ribo",
        "percent_malat1",
    ]
    attributes = ["activation"]

    params = Params()
    params.plot_raw = False

    sc.pp.scale(a)
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)

    # Plot

    # marker genes (tailored for a PBMC sample)
    mark1 = [
        "MALAT1",  # sc
        "CD34",  # HSC
        "CD3D",
        "CD3G",
        "CD247",  # T-cell
        "CD4",
        "FOXP3",
        "CCR7",
        "CTLA4",  # CD4
        "CD8A",
        "NKG7",
        "IL2RA",  # CD8
        "NCAM1",
        "GNLY",  # NK
        "CD14",
        "CST3",  # Monocytes
        "CD79A",
        "CD19",
        "IGHG1",  # B cells  (also MS4A1)
        "FCER1G",
        "CLEC10A",  # dendritic cells
        "SELE",
        "CD93",
        "PECAM1",
        "KDR",  # endothelial cells
        "DCN",
        "COL6A2",  # fibroblasts
        "GZMB",
        "CD68",
        "CD27",
        "MS4A1",
        "CD24",
        "NCR1",
        "CD274",
    ]  # Immune
    mark2 = [
        "IL32",
        "IFNG",
        "IFNGR1",
        "IL4R",
        "IL4",
        "JUN",
        "JUNB",
        "JUND",
        "JAK1",
        "JAK2",
        "GATA1",
        "JARID2",
        "KRAS",
        "MYC",
    ]
    mark3 = ["BTK", "LCK", "E2F4", "CXCR4", "ITGA4", "HBA1", "PTPRC"]
    red_mark = ["CST3", "TCL1A", "GZMB", "NKG7", "CD3D"]
    marker_genes = mark1 + mark2 + mark3

    sc.pl.pca_variance_ratio(a, log=True, show=False)
    plt.gca().figure.savefig(
        output_prefix + ".single_cell.pca_variance_ratio.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # fig = sc.pl.pca(a, color=tech_attributes + attributes, components=['1,2', '2,3', '3,4', '4,5'], return_fig=True)
    # for ax in fig.axes:
    #     ax.get_children()[0].set_rasterized(True)
    # fig.savefig(output_prefix + ".single_cell.pca.svg", dpi=300, bbox_inches="tight")

    fig = sc.pl.umap(a, color=tech_attributes + attributes, return_fig=True)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix + ".single_cell.umap.svg", dpi=300, bbox_inches="tight"
    )

    sc.pp.combat(a, "activation")

    sc.pp.scale(a)
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)
    fig = sc.pl.umap(a, color=tech_attributes + attributes, return_fig=True)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix + ".single_cell.post_combat.umap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    sc.pp.regress_out(a, keys=["log_counts"])

    sc.pp.scale(a)
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)
    fig = sc.pl.umap(a, color=tech_attributes + attributes, return_fig=True)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix
        + ".single_cell.post_combat.post_regressout_counts.umap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # sc.write(output_prefix + ".regress_out_counts.h5ad", a)

    # Now regress out also efficiency_ratio
    sc.pp.regress_out(a, keys=["efficiency_ratio"])

    sc.pp.scale(a)
    sc.pp.pca(a)
    sc.pp.neighbors(a)
    sc.tl.umap(a)
    fig = sc.pl.umap(a, color=tech_attributes + attributes, return_fig=True)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    fig.savefig(
        output_prefix
        + ".single_cell.post_combat.post_regressout_counts.post_regressout_efficiency.umap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # sc.write(output_prefix + ".regress_out_counts.regress_out_efficiency_ratio.h5ad", a)
    if params.plot_raw:
        g = [
            y
            for x in marker_genes
            for y in a.raw.var.loc[
                a.raw.var["external_gene_name"] == x
            ].index.tolist()
        ]
    else:
        g = [
            x for x in marker_genes if x in a.var["external_gene_name"].tolist()
        ]
    color = tech_attributes + attributes + g
    kwargs = dict(
        hspace=0.1,
        wspace=0,
        return_fig=True,
        use_raw=params.plot_raw,
        gene_symbols="external_gene_name" if not params.plot_raw else None,
    )

    fig = sc.pl.umap(a, color=color, **kwargs)
    for ax in fig.axes:
        ax.get_children()[0].set_rasterized(True)
    if params.plot_raw:
        for ax in fig.axes:
            try:
                ax.set_title(
                    a.raw.var.loc[ax.get_title(), "external_gene_name"]
                )
            except KeyError:
                pass
    fig.savefig(
        output_prefix
        + ".single_cell.post_combat.post_regressout_counts.post_regressout_efficiency.markers_extended.umap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # differntial expression

    max_genes = 500
    attribute = "activation"
    attributes = [attribute]
    a.X += abs(a.X.min())

    # # differential expression
    diff = differential_expression(a, attribute, n_genes=max_genes)
    diff.to_csv(
        output_prefix + f".{attribute}.cluster_comparison.top_values.csv",
        index=True,
    )

    diff = pd.read_csv(
        output_prefix + f".{attribute}.cluster_comparison.top_values.csv",
        index_col=0,
    )

    fig = plot_differential_expression(diff)
    fig.savefig(
        output_prefix + f".{attribute}.differential_expression.ma_plot.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # # differential enrichment
    groups = [x for x in diff["group"].unique() if x not in ["-1", -1]]

    enrichments = differential_enrichment(
        diff,
        groups,
        attribute,
        alpha=0.05,
        alpha_col="pvals_adj",
        max_n=max_genes,
        sort_by="scores",
    )
    enrichments.to_csv(
        output_prefix + f"{attribute}.differential_enrichment.csv", index=False
    )
    enrichments = pd.read_csv(
        output_prefix + f"{attribute}.differential_enrichment.csv", index_col=0
    )

    g = (
        enrichments.set_index("description")
        .groupby(["group"])["combined_score"]
        .nlargest(5)
    )
    print("combined_score:\n", g)
    g = (
        enrichments.set_index("description")
        .groupby(["group"])["p_value"]
        .nsmallest(5)
    )
    print("p_value:\n", g)

    plot_differential_enrichment(enrichments, output_prefix, ntop_terms=20)
