#!/usr/bin/env python

"""
Formerly "src/20200525.analysis.py"
"""

import os
import time
from os.path import join as pjoin
import json
from pathlib import Path

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


sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)


class Params:
    pass


@lru_cache(maxsize=None)
def query(species, ensembl_version):
    return query_biomart(
        attributes=["ensembl_gene_id", "external_gene_name"],
        species=species,
        ensembl_version=ensembl_version,
    )


def download_unsafe(url, output_file):
    from urllib.request import urlopen
    import ssl

    response = urlopen(url, context=ssl._create_unverified_context())
    with open(output_file, "wb") as handle:
        handle.write(response.read())


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
    elif file_format == "loom":
        a.write_loom(output_file)
    return a


def differential_expression(
    an,
    attribute,
    groups="all",
    method="t-test_overestim_var",
    n_genes=250,
    use_raw=False,
):
    sc.tl.rank_genes_groups(
        an,
        attribute,
        groups=groups,
        method=method,
        n_genes=n_genes,
        use_raw=use_raw,
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
    if m.min() < 0:
        mean = pd.Series(m, index=an.var.index, name="mean")
    else:
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
    attributes,
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
            "ARCHS4_Cell-lines",
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
    enrichments,
    output_prefix,
    gene_set_libraries=None,
    ntop_terms=5,
    robust=False,
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

                if not robust:
                    vmax = combs.max().max()
                    vmax += vmax * 0.15
                    kwargs = dict(vmax=vmax)
                else:
                    kwargs = dict(robust=True)
                grid = sns.clustermap(
                    combs2,
                    metric="correlation",
                    cbar_kws={"label": f"{metric}"},
                    square=True,
                    col_cluster=False,
                    row_cluster=False if n < 10 else True,
                    yticklabels=True,
                    **kwargs,
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
                    **kwargs,
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
        ax.set_ylabel("Sex ratio estimate (Y / X)")
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


def get_relevant_annotation(annotations, sample):
    try:
        annotation = annotations.loc[
            annotations["sample_name"].str.contains(sample)
        ].drop(["sample_name", "combinatorial_barcode"], axis=1)
        if annotation.empty:
            raise FileNotFoundError
        # drop empty columns
        annotation = annotation.drop(
            annotation.columns[annotation.isnull().all()], axis=1
        )
        # drop single columns
        annotation = annotation.drop(
            annotation.columns[
                annotation.apply(pd.Series.nunique, axis=0) == 1
            ],
            axis=1,
        )
        # get attribute names
        attributes = annotation.drop(["plate_well"], axis=1).columns.tolist()
        if (
            annotation.shape[1] == 1
            and (annotation.columns == ["plate_well"]).all()
        ):
            annotation = None
    except FileNotFoundError:
        annotation = None
        attributes = []
    return annotation, attributes


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


class Args:
    pass


args = Args()
args.droplet_column = "r2"


base_url = "http://biomedical-sequencing.at/projects/BSA_0383_PD_SCIFI_72203e85f90a4b0d8ad6cee221d59e6d"
samples1 = [
    "T46_10xGenomics_1_4lines_7650_nuclei",
    "T46_10xGenomics_2_4lines_7650_MeOH-cells",
    "T46_10xGenomics_3_4lines_7650_intact-cells",
]
samples2 = [
    "PD212_scifi_1N_4lines_7650_nuclei",
    "PD212_scifi_2C_4lines_7650_MeOH-cells",
    "PD213_scifi_1_CRISPR-TCR_15300_MeOH-cells",
    "PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells",
    "PD213_scifi_3_CRISPR-TCR_15300_nuclei",
]
samples = samples1 + samples2

annotations = pd.read_csv(
    Path(
        "metadata", "PD000212_PD000213_multiplexing_barcodes.scifi_specific.csv"
    )
).drop(["multiplexing_barcode", "gRNA_well", "gRNA_seq"], axis=1)
annotations["gRNA"] = (
    annotations["gRNA_ID"]
    .str.split("_")
    .apply(lambda x: x[0] if isinstance(x, list) else np.nan)
    .str.replace("CTRL.*", "CTRL")
)

# Analysis parameters
parameters = json.load(
    open(pjoin("metadata", "PD212-213.analysis_thresholds.json"), "r")
)


# Outputs
run = "results/202005_results"
os.makedirs(run, exist_ok=True)

#
# r1_barcodes = pd.read_csv("metadata/737K-april-2014_rc.txt", header=None, squeeze=True)
# r1_barcodes = pd.read_csv("metadata/737K-august-2016.txt", header=None, squeeze=True)
r1_barcodes = pd.read_csv(
    "metadata/4M-with-alts-february-2016.txt", header=None, squeeze=True
)
r2_barcodes = pd.read_csv(
    "metadata/737K-cratac-v1.txt", header=None, squeeze=True
)

REPL = {
    "REH": "NALM6",
    "JURKAT": "Jurkat-Cas9-TCR",
    "HEK293": "HEK293T",
    "K562": "K562",
}


def download_write():

    # # Metrics
    for sample in samples1 + samples2:
        for ftype in ["metrics", "expression"]:
            print(sample, ftype)
            metrics_url = (
                "/".join([base_url, "data", sample, sample])
                + f".{ftype}.csv.gz"
            )
            metrics_file = Path("data", sample, sample + f".{ftype}.csv.gz")
            metrics_file.parent.mkdir(exist_ok=True)
            if not os.path.exists(metrics_file):
                download_unsafe(metrics_url, metrics_file)

    # # Expression
    for sample in samples:
        print(sample)
        expression_file = Path("data", sample, sample + ".expression.csv.gz")
        h5ad_file = Path("data", sample, sample + ".expression.h5ad")
        args.output_prefix = Path(
            "data", sample, sample + ".expression."
        ).as_posix()
        if os.path.exists(h5ad_file):
            continue

        annotation, attributes = get_relevant_annotation(annotations, sample)
        output_prefix = pjoin("202005_results", sample + ".")
        print(f"# {time.asctime()} - Reading expression CSV file.")
        expr = pd.read_csv(expression_file)
        a = write_gene_expression_matrix(expr, annotation=annotation)


def quick_qc():
    # Quick QC plots
    if not os.path.exists("20200525.first_qc.svg"):
        metrics = dict()
        for sample in samples1:
            metrics_file = Path("data", sample, sample + ".metrics.csv.gz")
            m = (
                pd.read_csv(metrics_file)
                .set_index("r1")
                .sort_values("umi", ascending=False)
            )
            # m.index += 1
            metrics[sample] = m
        for sample in samples2:
            metrics_file = Path("data", sample, sample + ".metrics.csv.gz")
            m = pd.read_csv(metrics_file)
            nulls = m.columns[m.isnull().all()].tolist()
            if nulls:
                cols = m.columns.drop(["r2"] + nulls[1:])
                m = m.drop(nulls, axis=1)
                m.columns = cols
            metrics[sample] = m.sort_values("umi", ascending=False)

        fig, axis = plt.subplots(
            2, len(samples), figsize=(6 * len(samples), 4 * 2)
        )
        for i, sample in enumerate(samples):
            m = metrics[sample]
            axis[0, i].set_title(sample)
            axis[0, i].plot(m.index, m["umi"])
            axis[1, i].scatter(
                m["read"],
                m["unique_fraction"],
                s=1,
                alpha=0.05,
                rasterized=True,
            )
            axis[0, i].loglog()
            axis[1, i].set_xscale("log")

        axis[0, 0].set_ylabel("Barcodes")
        axis[1, 0].set_ylabel("Unique fraction")
        axis[0, len(samples) // 2].set_xlabel("Barcodes")
        axis[1, len(samples) // 2].set_xlabel("Reads per barcode")
        fig.savefig("20200525.first_qc.svg", dpi=300, bbox_inches="tight")

        # 20200720
        # replot with right colors
        from src.scifi_utils import inflection_point

        output_dir = Path(run) / "qc"

        fig, ax = plt.subplots(1, 1, figsize=(2 * 1.372, 4))
        for sample in [s for s in samples if "lines" in s]:
            print(sample)
            # a = sc.read(Path(run) / sample / (sample + ".filtered.h5ad"))
            m = metrics[sample]
            if "10x" not in sample:
                m["r1"] = m["plate_well"] + "-" + m.index

            m["rank"] = m["umi"].rank(ascending=False)
            ip = inflection_point(m["umi"])
            fore = m.iloc[:ip]
            back = m.iloc[ip:]

            c = "#3aa29a" if "10xGenomics" in sample else "#931e80"
            ls = (
                "-"
                if "intact" in sample
                else "--"
                if "MeOH" in sample
                else "-."
            )
            label = f"{sample}, n = {ip}"

            ax.plot(back["rank"], back["umi"], color="grey", linestyle=ls)
            ax.plot(
                fore["rank"], fore["umi"], color=c, label=label, linestyle=ls
            )

        ax.set(xscale="log", yscale="log", xlabel="Barcodes", ylabel="UMIs")
        fig.legend(loc="lower center")
        sns.despine(fig)
        fig.savefig(
            output_dir / "20200720.Fig2e.kneeplot.svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)


def process():

    # for sample in samples:
    for sample in [s for s in samples if "lines" in s]:
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
        annotation, attributes = get_relevant_annotation(annotations, sample)
        # expression_file = pjoin("~/Downloads", sample + ".expression.csv.gz")
        # metrics_file = pjoin("~/Downloads", sample + ".metrics.csv.gz")

        a = sc.read(Path("data", sample, sample + ".expression.h5ad"))
        # This filtering here is just to get rid of the low end of barcodes
        sc.pp.filter_cells(a, min_counts=50)

        if "10x" not in sample:
            a.obs["r1"] = a.obs["plate_well"] = a.obs.index.str.slice(0, 3)
            a.obs["r2"] = a.obs.index.str.slice(4)
            r2_in = a.obs["r2"].isin(r2_barcodes)
            print(
                f"{(r2_in.sum() / a.obs.shape[0]) * 100:.3f}% barcodes in reference"
            )
            a = a[r2_in]
        else:
            pass
            # a.obs['r1'] = a.obs.index
            # a = a[a.obs['r1'].isin(r1_barcodes)]
        # a = a.copy()
        # a.obs = a.obs.reset_index().set_index("plate_well").join(
        #     annotation.set_index("plate_well")).reset_index().set_index("index")

        if "Mouse" in sample:
            try:
                biomart = pd.read_csv(
                    pjoin("metadata", "mouse.biomart.csv"), index_col=0
                )
            except FileNotFoundError:
                biomart = query(
                    species="mmusculus", ensembl_version="grch38"
                ).set_index("ensembl_gene_id")
                biomart.to_csv(pjoin("metadata", "mouse.biomart.csv"))
        else:
            try:
                biomart = pd.read_csv(
                    pjoin("metadata", "human.biomart.csv"), index_col=0
                )
            except FileNotFoundError:
                biomart = query(
                    species="hsapiens", ensembl_version="grch38"
                ).set_index("ensembl_gene_id")
                biomart.to_csv(pjoin("metadata", "human.biomart.csv"))
        for col in biomart.columns:
            biomart[col] = biomart[col].replace("nan", np.nan)
        a.var = a.var.join(biomart)
        # replace genes with no gene symbol with original Ensembl ID
        null = a.var.loc[
            a.var["external_gene_name"].isnull(), "external_gene_name"
        ]
        null.update(null.index.to_series())
        a.var.loc[
            a.var["external_gene_name"].isnull(), "external_gene_name"
        ] = null

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

        a.obs.loc[:, "efficiency_ratio"] = (
            a.obs["log_genes"] / a.obs["log_counts"]
        )

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
        a.obs.loc[:, "percent_malat1"] = (
            a.obs["malat1"] / a.obs["n_counts"]
        ) * 100

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

        if not os.path.exists(
            output_prefix + ".quality_metrics.pre_filtering.svg"
        ):
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
                output_prefix + ".quality_metrics.pre_filtering.svg",
                bbox_inches="tight",
                dpi=300,
            )

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

        # a.raw = a  # this was the previous place to store raw
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
            output_prefix + ".quality_metrics.post_filtering.svg",
            bbox_inches="tight",
            dpi=300,
        )

        a.raw = a

        sc.pp.normalize_per_cell(a, counts_per_cell_after=1e4)
        sc.pp.log1p(a)

        # sc.pp.regress_out(a, keys=["n_counts"], n_jobs=4)

        # sc.pp.scale(a, max_value=10)
        sc.pp.highly_variable_genes(a)
        sc.pp.pca(a, zero_center=True, use_highly_variable=True)
        sc.pp.neighbors(a)
        sc.tl.umap(a)

        sc.tl.leiden(a, resolution=params.leiden["resolution"])
        if "leiden" not in attributes:
            attributes += ["leiden"]

        # # Add sex
        # if "Mouse" not in sample:
        #     attributes += ["sex_ratio"]
        #     a, fig = add_sex_ratio(a, plot=True)
        #     fig.savefig(
        #         output_prefix + ".single_cell.sex_ratio_estimate_over_mean.svg",
        #         dpi=300,
        #         bbox_inches="tight",
        #     )

        #     fig, ax = plt.subplots(1, 1, figsize=(3, 3))
        #     groupby = "donor_sex" if "donor_sex" in a.obs.columns else "leiden"
        #     sc.pl.violin(a, groupby=groupby, keys="sex_ratio", ax=ax, show=False)
        #     ax.axhline(0, linestyle="--", color="grey")
        #     fig.savefig(output_prefix + ".single_cell.sex_ratio.svg", dpi=300, bbox_inches="tight")

        # if not os.path.exists(output_prefix + ".filtered.h5ad"):
        sc.write(output_prefix + ".filtered.h5ad", a)
        # if not os.path.exists(output_prefix + ".filtered.obs.csv"):
        a.obs.to_csv(output_prefix + ".filtered.obs.csv")

        # Make sure plotting order is random
        a = a[
            np.random.choice(
                a.obs.index.tolist(), a.obs.shape[0], replace=False
            ),
            :,
        ]

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
            output_prefix + ".single_cell.umap.svg",
            dpi=300,
            bbox_inches="tight",
        )

        if "lines" not in sample:
            if params.plot_raw:
                g = [
                    y
                    for x in mark1
                    for y in a.raw.var.loc[
                        a.raw.var["external_gene_name"] == x
                    ].index.tolist()
                ]
            else:
                g = [
                    x
                    for x in mark1
                    if x in a.var["external_gene_name"].tolist()
                ]
            color = tech_attributes + attributes + g
            kwargs = dict(
                hspace=0.1,
                wspace=0,
                return_fig=True,
                use_raw=params.plot_raw,
                gene_symbols="external_gene_name"
                if not params.plot_raw
                else None,
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
                    gene_symbols="external_gene_name"
                    if not params.plot_raw
                    else None,
                )

                fig = sc.pl.umap(a, color=color, **kwargs)
                for ax in fig.axes:
                    ax.get_children()[0].set_rasterized(True)
                if params.plot_raw:
                    for ax in fig.axes:
                        try:
                            ax.set_title(
                                a.raw.var.loc[
                                    ax.get_title(), "external_gene_name"
                                ]
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
                        gene_symbols="external_gene_name"
                        if not plot_raw
                        else None,
                    )

                    fig = sc.pl.umap(a, color=color, **kwargs)
                    for ax in fig.axes:
                        ax.get_children()[0].set_rasterized(True)
                    if plot_raw:
                        for ax in fig.axes:
                            try:
                                ax.set_title(
                                    a.raw.var.loc[
                                        ax.get_title(), "external_gene_name"
                                    ]
                                )
                            except KeyError:
                                pass
                    fig.savefig(
                        output_prefix
                        + f".single_cell.markers.{label}.umap.svg",
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

        if "variables" in params.diff_exp:
            diff_attributes = params.diff_exp["variables"]
        else:
            diff_attributes = attributes

        for attribute in diff_attributes:
            if not isinstance(a.obs[attribute].dtype, pd.CategoricalDtype):
                continue
            if a.obs[attribute].nunique() == 1:
                continue
            # # differential expression
            diff = differential_expression(
                a2, attribute, n_genes=params.diff_exp["max_genes"]
            )
            diff.to_csv(
                output_prefix
                + f".{attribute}.cluster_comparison.top_values.csv",
                index=True,
            )

            diff = pd.read_csv(
                output_prefix
                + f".{attribute}.cluster_comparison.top_values.csv",
                index_col=0,
            )

            fig = plot_differential_expression(diff)
            fig.savefig(
                output_prefix
                + f".{attribute}.differential_expression.ma_plot.svg",
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
                alpha_col="pvals_adj"
                if "threshold_variable" not in params.diff_exp
                else params.diff_exp["threshold_variable"],
                max_n=params.diff_exp["max_genes"],
                sort_by="scores",
            )
            enrichments.to_csv(
                output_prefix + f".{attribute}.differential_enrichment.csv",
                index=False,
            )
            enrichments = pd.read_csv(
                output_prefix + f".{attribute}.differential_enrichment.csv",
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


def fourlines():
    side = 3.1

    output_dir = Path(run, "qc")
    output_dir.mkdir(exist_ok=True)

    metrics = dict()
    for sample in [s for s in samples1 if "lines" in s]:
        metrics_file = Path("data", sample, sample + ".metrics.csv.gz")
        m = (
            pd.read_csv(metrics_file)
            .set_index("r1")
            .sort_values("umi", ascending=False)
        )
        m["efficiency"] = m["umi"] / m["read"]
        metrics[sample] = m
    for sample in [s for s in samples2 if "lines" in s]:
        metrics_file = Path("data", sample, sample + ".metrics.csv.gz")
        m = pd.read_csv(metrics_file)
        nulls = m.columns[m.isnull().all()].tolist()
        if nulls:
            cols = m.columns.drop(["r2"] + nulls[1:])
            m = m.drop(nulls, axis=1)
            m.columns = cols
        m = m.loc[m.index.isin(r2_barcodes)]
        m["efficiency"] = m["umi"] / m["read"]
        m.index = m["plate_well"] + "-" + m.index
        metrics[sample] = m.sort_values("umi", ascending=False)

    ss = [s for s in samples if "lines" in s]

    metrics_dict = metrics
    fig, axes = plt.subplots(
        len(ss),
        4,
        figsize=(side * 4, side * len(ss)),
        sharex="col",
        sharey="col",
    )
    for i, sample in enumerate(ss):
        m = metrics_dict[sample]
        axes[i, 0].set_title(sample)
        # knee
        # # reads, umis, genes
        for variable in ["read", "umi", "gene"]:
            mt = m.sort_values(variable, ascending=False)
            axes[i, 0].plot(range(mt.shape[0]), mt[variable])
        axes[i, 0].loglog()
        axes[i, 0].set_xlabel("Barcodes")
        axes[i, 0].set_ylabel("Value")
        # reads vs UMI
        v = max(m["read"].max(), m["umi"].max())
        axes[i, 1].plot([1, v], [1, v], linestyle="--", color="grey")
        axes[i, 1].scatter(
            data=m,
            x="read",
            y="umi",
            c="efficiency",
            s=2,
            alpha=0.1,
            rasterized=True,
        )
        axes[i, 1].loglog()
        axes[i, 1].set_xlabel("Reads")
        axes[i, 1].set_ylabel("UMIs")
        # UMIs vs genes
        v = max(m["umi"].max(), m["gene"].max())
        axes[i, 2].plot([1, v], [1, v], linestyle="--", color="grey")
        axes[i, 2].scatter(
            data=m,
            x="umi",
            y="gene",
            c="efficiency",
            s=2,
            alpha=0.1,
            rasterized=True,
        )
        axes[i, 2].loglog()
        axes[i, 2].set_xlabel("Barcodes")
        axes[i, 2].set_ylabel("Unique fraction")
        # Reads vs unique fraction
        axes[i, 3].scatter(
            data=m,
            x="umi",
            y="unique_fraction",
            c="unique_fraction",
            s=1,
            alpha=0.05,
            rasterized=True,
        )
        axes[i, 3].set_xlabel("Reads")
        axes[i, 3].set_ylabel("Unique fraction")
        axes[i, 3].set_xscale("log")
    fig.savefig(
        output_dir / "4lines.metrics.all_cell_lines.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # get cell type labels for 10X samples based on clusting and its indentity in enrichments
    # it's done for both 10X and scifi just to have exactly the same procedure
    labels = dict()
    for sample in [s for s in samples if "lines" in s]:
        labels[sample] = label_clusters(sample, repl=REPL)

    for cell_line in REPL.values():
        print(cell_line)

        metrics_dict = dict()
        for sample, metric in metrics.items():
            print(sample)
            metrics_dict[sample] = metric.loc[labels[sample][cell_line]]

        fig, axes = plt.subplots(
            len(ss),
            4,
            figsize=(side * 4, side * len(ss)),
            sharex="col",
            sharey="col",
        )
        for i, sample in enumerate(ss):
            m = metrics_dict[sample]
            axes[i, 0].set_title(sample)
            # knee
            # # reads, umis, genes
            for variable in ["read", "umi", "gene"]:
                mt = m.sort_values(variable, ascending=False)
                axes[i, 0].plot(range(mt.shape[0]), mt[variable])
            axes[i, 0].loglog()
            axes[i, 0].set_xlabel("Barcodes")
            axes[i, 0].set_ylabel("Value")
            # reads vs UMI
            vmin = min(m["read"].min(), m["umi"].min())
            vmax = max(m["read"].max(), m["umi"].max())
            axes[i, 1].plot(
                [vmin, vmax], [vmin, vmax], linestyle="--", color="grey"
            )
            axes[i, 1].scatter(
                data=m,
                x="read",
                y="umi",
                c="efficiency",
                s=2,
                alpha=0.1,
                rasterized=True,
            )
            axes[i, 1].loglog()
            axes[i, 1].set_xlabel("Reads")
            axes[i, 1].set_ylabel("UMIs")
            # UMIs vs genes
            vmin = min(m["umi"].min(), m["gene"].min())
            vmax = max(m["umi"].max(), m["gene"].max())
            axes[i, 2].plot(
                [vmin, vmax], [vmin, vmax], linestyle="--", color="grey"
            )
            axes[i, 2].scatter(
                data=m,
                x="umi",
                y="gene",
                c="efficiency",
                s=2,
                alpha=0.1,
                rasterized=True,
            )
            axes[i, 2].loglog()
            axes[i, 2].set_xlabel("UMIs")
            axes[i, 2].set_ylabel("Unique fraction")
            # Reads vs unique fraction
            axes[i, 3].scatter(
                data=m,
                x="umi",
                y="unique_fraction",
                c="unique_fraction",
                s=1,
                alpha=0.05,
                rasterized=True,
            )
            axes[i, 3].set_xlabel("Reads")
            axes[i, 3].set_ylabel("Unique fraction")
            axes[i, 3].set_xscale("log")
            axes[i, 3].set_ylim((0, 1))
        fig.savefig(
            output_dir / f"4lines.metrics.{cell_line}.svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close("all")
        # gc.collect()


def fourlines_replot_umap_colors():
    from collections import defaultdict

    colors = {"1": "#2278b5", "2": "#ef7f20", "0": "#279e69", "3": "#d62728"}

    output_dir = Path(run) / "qc"
    output_dir.mkdir(exist_ok=True)

    for sample in [s for s in samples if "lines" in s]:
        a = sc.read(Path(run) / sample / (sample + ".filtered.h5ad"))

        fig, ax = plt.subplots(1, 1, figsize=(4 * 1.372, 4))
        for i, cl in enumerate(a.obs["leiden"].unique()):
            df = a.obs.join(
                pd.DataFrame(
                    a.obsm["X_umap"], index=a.obs.index, columns=["0", "1"]
                )
            ).query(f"leiden == '{cl}'")
            ax.scatter(
                data=df,
                x="0",
                y="1",
                c=colors.get(cl, "#878787"),
                rasterized=True,
                s=2,
                alpha=0.5,
            )
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.savefig(
            output_dir / (sample + ".umap.svg"), dpi=300, bbox_inches="tight"
        )
        plt.close(fig)


def fourlines_representation():
    from scipy.stats import mannwhitneyu

    output_dir = Path(run, "qc")
    output_dir.mkdir(exist_ok=True)

    ss = [s for s in samples if "lines" in s]
    cell_lines = dict()
    for sample in ss:
        cell_lines[sample] = {
            k: len(v) for k, v in label_clusters(sample).items()
        }

    df = pd.DataFrame.from_dict(cell_lines).rename_axis(
        index="cell_line", columns="sample"
    )
    df = (df / df.sum()) * 100
    df = df.reset_index().melt(id_vars="cell_line", value_name=r"% cells")

    df["error"] = (df["% cells"] - 25) ** 2
    df["10x"] = df["sample"].str.contains("10x")
    np.sqrt(df.groupby("10x")["error"].mean())

    stat, p = mannwhitneyu(
        df.loc[df["sample"].str.contains("10x"), "error"],
        df.loc[~df["sample"].str.contains("10x"), "error"],
    )

    for label, ddf in [
        ("", df),
        (
            ".exclude_intact_cells",
            df.query("~sample.str.contains('intact')", engine="python"),
        ),
    ]:
        fig, axis = plt.subplots(figsize=(6, 4))
        axis = sns.barplot(data=ddf, x="cell_line", y=r"% cells", hue="sample")
        axis.axhline(25, linestyle="--", color="grey")
        fig.savefig(
            output_dir / f"4lines.cell_line_representation{label}.svg",
            dpi=300,
            bbox_inches="tight",
        )

        fig, axis = plt.subplots(figsize=(6, 4))
        axis = sns.barplot(data=ddf, x="sample", y=r"% cells", hue="cell_line")
        axis.axhline(25, linestyle="--", color="grey")
        fig.savefig(
            output_dir
            / f"4lines.cell_line_representation.groupby_sample{label}.svg",
            dpi=300,
            bbox_inches="tight",
        )

    fig, axis = plt.subplots(figsize=(6, 4))
    axis = sns.barplot(data=df, x="sample", y="error")
    axis.axhline(25, linestyle="--", color="grey")
    fig.savefig(
        output_dir / f"4lines.cell_line_representation{label}.svg",
        dpi=300,
        bbox_inches="tight",
    )


# HEK293T: #2278b5, Jurkat: #ef7f20, K562: #279e69, NALM6: #d62728


def transcriptome_similarity():
    # import pyupset as pyu

    output_dir = Path(run, "qc")
    output_dir.mkdir(exist_ok=True)

    ss = [s for s in samples if "lines" in s]

    # correlate expression aggregated for each cluster
    sample_expr = list()
    for sample in ss:
        a = sc.read(pjoin(run, sample, sample) + ".filtered.h5ad")
        expr = a.to_df().join(a.obs["leiden"]).groupby("leiden").mean()
        expr = expr.loc[expr.index.astype(int) < 4]
        lab = {
            v: k for k, v in label_clusters(sample, return_labels=True).items()
        }
        expr.index = [lab[int(x)] for x in expr.index]
        expr.index = sample + " - " + expr.index.astype(str)
        sample_expr.append(expr)

    expr = pd.concat(sample_expr).T.dropna()
    # expr = expr.loc[:, expr.columns.str.split(" - ").map(lambda x: int(x[1])) < 4]

    grid = sns.clustermap(
        expr.corr(),
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws=dict(label="Pearson correlation"),
    )
    grid.fig.savefig(
        output_dir / "4lines.cell_line.transcriptome_correlation.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # overlap of top genes

    # # across all cell lines
    cell_lines = dict()
    for sample in ss:
        diff = pd.read_csv(
            pjoin(run, sample, sample)
            + ".leiden.cluster_comparison.top_values.csv",
            index_col=0,
        )
        cell_lines[sample] = diff.index

    over = pd.DataFrame()
    for sample1 in cell_lines:
        for sample2 in cell_lines:
            union = len(
                set(
                    cell_lines[sample1].unique().tolist()
                    + cell_lines[sample2].unique().tolist()
                )
            )
            over.loc[sample1, sample2] = (
                cell_lines[sample1]
                .unique()
                .isin(cell_lines[sample2].unique())
                .sum()
                / union
            )

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    sns.heatmap(
        over * 100,
        cbar_kws=dict(label="% top 100 differential genes overlap"),
        annot=True,
        square=True,
        fmt=".1f",
        ax=ax,
        vmin=0,
        vmax=100,
    )
    fig.savefig(
        output_dir / "4lines.cell_line.top_100_gene_overlap.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # # for each cell line
    top_overlap = {
        "K562": {},
        "NALM6": {},
        "Jurkat-Cas9-TCR": {},
        "HEK293T": {},
    }
    for sample in ss:
        labels = label_clusters(sample, return_labels=True)
        rev_labels = {v: k for k, v in labels.items()}
        diff = pd.read_csv(
            pjoin(run, sample, sample)
            + ".leiden.cluster_comparison.top_values.csv",
            index_col=0,
        )
        diff["group"] = diff["group"].replace(rev_labels)
        for cl in top_overlap:
            top_overlap[cl][sample] = diff.query(
                f"group == '{cl}'"
            ).index.drop_duplicates()

    _over = dict()
    for cl in top_overlap:
        over = pd.DataFrame()
        for sample1 in ss:
            for sample2 in ss:
                union = len(
                    set(top_overlap[cl][sample1] + top_overlap[cl][sample2])
                )
                over.loc[sample1, sample2] = (
                    top_overlap[cl][sample1]
                    .isin(top_overlap[cl][sample2])
                    .sum()
                    / union
                )
        _over[cl] = over.assign(cell_line=cl)

    fig, axes = plt.subplots(
        2, 2, figsize=(2 * 3, 2 * 3), sharex=True, sharey=True
    )
    for ax, cl in zip(axes.flat, _over):
        sns.heatmap(
            _over[cl].drop("cell_line", 1) * 100,
            cbar_kws=dict(label="% top 100 differential genes overlap"),
            annot=True,
            square=True,
            fmt=".1f",
            ax=ax,
            vmin=0,
            vmax=100,
        )
        ax.set(title=cl)
    fig.savefig(
        output_dir / "4lines.cell_line.top_100_gene_overlap.by_cell_line.svg",
        dpi=300,
        bbox_inches="tight",
    )

    top_overlap2 = {
        k + " - " + k2: v2
        for k, v in top_overlap.items()
        for k2, v2 in v.items()
    }

    over = pd.DataFrame()
    for g1 in top_overlap2:
        for g2 in top_overlap2:
            union = len(set(top_overlap2[g1] + top_overlap2[g2]))
            over.loc[g1, g2] = (
                top_overlap2[g1].isin(top_overlap2[g2]).sum() / union
            )
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    sns.heatmap(
        (over * 100).astype(int),
        vmin=0,
        vmax=100,
        annot=True,
        ax=ax,
        square=True,
    )
    fig.savefig(
        output_dir
        / "4lines.cell_line.top_100_gene_overlap.by_cell_line-joint.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # Across all genes
    # # for each cell line
    for metric in ["scores", "logfoldchanges", "pvals"]:
        sample_expr = list()
        _diffs = list()
        for sample in ss:
            a = sc.read(pjoin(run, sample, sample) + ".filtered.h5ad")
            sc.tl.rank_genes_groups(a, "leiden", method="t-test_overestim_var")
            labels = label_clusters(sample, return_labels=True)
            rev_labels = {str(v): k for k, v in labels.items()}

            nms = pd.DataFrame.from_records(a.uns["rank_genes_groups"]["names"])
            scr = pd.DataFrame.from_records(a.uns["rank_genes_groups"][metric])
            _diffs.append(
                pd.concat(
                    [
                        nms[cl]
                        .to_frame(name="gene")
                        .join(scr[cl].rename(metric))
                        .assign(cell_line=rev_labels[cl])
                        for cl in nms.columns[
                            nms.columns.isin(rev_labels.keys())
                        ]
                    ]
                ).assign(sample=sample)
            )

        diffs = pd.concat(_diffs)

        if metric == "pvals":
            # avoid infinite when pvals == 0
            v = diffs[metric].replace(0, 100).min()
            diffs[metric] = -np.log10(diffs[metric].replace(0, v))

        # correlate scores
        fig, axes = plt.subplots(
            2, 2, figsize=(3 * 2, 3 * 2), sharex=True, sharey=True
        )
        for ax, cl in zip(axes.flat, diffs["cell_line"].unique()):
            corrs = (
                diffs.query(f"cell_line == '{cl}'")
                .pivot_table(index="gene", columns="sample", values=metric)
                .corr()
            )  # .corr(method='spearman') -> too sloooow
            corrs = corrs.reindex(ss, axis=0).reindex(ss, axis=1)
            sns.heatmap(
                corrs,
                cbar_kws=dict(
                    label=f"Pearson correlation in differential expression {metric}"
                ),
                annot=True,
                square=True,
                fmt=".1f",
                ax=ax,
            )
            ax.set(title=cl, xlabel=None, ylabel=None)
        fig.savefig(
            output_dir
            / f"4lines.cell_line.correlation_in_diffexp_{metric}.by_cell_line.svg",
            dpi=300,
            bbox_inches="tight",
        )

        d = diffs.pivot_table(
            index="gene", columns=["cell_line", "sample"], values=metric
        )
        # grid = sns.clustermap(d.corr(), center=0, cmap="coolwarm")

        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        sns.heatmap(
            d.corr(),
            center=0,
            cmap="coolwarm",
            ax=ax,
            square=True,
            vmin=-1,
            vmax=1,
        )
        fig.savefig(
            output_dir
            / f"4lines.cell_line.correlation_in_diffexp_{metric}.joint.svg",
            dpi=300,
            bbox_inches="tight",
        )


def joint_embedding_across_experiments():
    ss = [s for s in samples if "lines" in s]
    output_dir = Path(run) / "qc"

    # Make one anndata with all experiments
    anndatas = list()
    for sample in ss:
        a = sc.read(pjoin(run, sample, sample) + ".filtered.h5ad")
        expr = a.to_df().join(a.obs["leiden"]).groupby("leiden").mean()
        expr = expr.loc[expr.index.astype(int) < 4]
        lab = {
            str(v): k
            for k, v in label_clusters(sample, return_labels=True).items()
        }
        a.obs["leiden_name"] = a.obs["leiden"].replace(lab)
        a.obs["sample"] = sample
        a = a[a.obs["leiden_name"].str.len() > 1, :].copy()
        anndatas.append(a)

    ann = sc.concat(anndatas)

    # # original embeddings (just an overlay, no meaning)
    # axes = sc.pl.pca(ann, color=['sample', 'leiden_name'], show=False)
    # axes = sc.pl.umap(ann, color=['sample', 'leiden_name'], show=False)

    # new joint embedding
    sc.pp.pca(ann)
    axes = sc.pl.pca(ann, color=["sample", "leiden_name"], show=False)
    for ax in axes:
        ax.get_children()[0].set(rasterized=True)
    axes[0].figure.savefig(
        output_dir / f"4lines.cell_line.joint_embedding.pca.svg",
        dpi=300,
        bbox_inches="tight",
    )

    sc.pp.neighbors(ann)
    sc.tl.umap(ann)

    # let's try other embeddings just for fun
    sc.tl.diffmap(ann)
    sc.tl.tsne(ann, n_jobs=12)
    sc.tl.draw_graph(ann, n_jobs=12)

    sc.write(run + "/all_experiments.h5ad", ann)

    # let's compute silhouette scores for all embeddings between samples and between cell lines
    variables = ["sample", "leiden_name"]
    spaces = ["pca", "umap", "diffmap", "tsne", "draw_graph_fa"]
    scores = pd.DataFrame(index=variables, columns=spaces, dtype=float)
    for var in variables:
        for space in spaces:
            print(f"{var} - {space}")
            scores.loc[var, space] = sklearn.metrics.silhouette_score(
                ann.obsm[f"X_{space}"][:, :2], ann.obs[var]
            )

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.heatmap(
        scores.T,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        annot=True,
        cbar_kws=dict(label="Silhouette score"),
        ax=ax,
    )
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    fig.savefig(
        output_dir / f"4lines.cell_line.joint_embedding.silhouette_scores.svg",
        dpi=300,
        bbox_inches="tight",
    )

    # Plot the embeddings
    for method in [
        sc.pl.pca,
        sc.pl.umap,
        sc.pl.diffmap,
        sc.pl.tsne,
        sc.pl.draw_graph,
    ]:
        name = method.__name__
        name = name if "draw" not in name else name + "_fa"
        score = scores[name]

        axes = method(ann, color=["sample", "leiden_name"], show=False)

        for ax in axes:
            ax.get_children()[0].set(rasterized=True)
        axes[0].set(title=f"Silhouette:\nSamples: {score['sample']:.3f}")
        axes[1].set(title=f"Silhouette:\nCell line: {score['leiden_name']:.3f}")

        axes[0].figure.savefig(
            output_dir / f"4lines.cell_line.joint_embedding.{name}.svg",
            dpi=300,
            bbox_inches="tight",
        )

        sc.tl.embedding_density(ann, basis=name, groupby="leiden_name")
        axes = sc.pl.embedding_density(
            ann, basis=name, key="umap_density_leiden_name", show=False
        )
        for ax in axes:
            ax.get_children()[0].set(rasterized=True)
        axes[0].figure.savefig(
            output_dir
            / f"4lines.cell_line.joint_embedding.{name}.density.by_cell_line.svg",
            dpi=300,
            bbox_inches="tight",
        )

    # add gene names to anndata
    try:
        biomart = pd.read_csv(
            pjoin("metadata", "human.biomart.csv"), index_col=0
        )
    except FileNotFoundError:
        biomart = query(species="hsapiens", ensembl_version="grch38").set_index(
            "ensembl_gene_id"
        )
        biomart.to_csv(pjoin("metadata", "human.biomart.csv"))
    for col in biomart.columns:
        biomart[col] = biomart[col].replace("nan", np.nan)
    ann.var = ann.var.join(biomart)

    # Differential expression
    for variable in ["sample", "leiden_name"]:
        diff_dir = output_dir / f"differential_{variable}"
        diff_dir.mkdir(exist_ok=True)
        output_prefix = (diff_dir / f"differential_{variable}.").as_posix()

        diff = differential_expression(
            ann, variable, n_genes=5000, use_raw=True
        )
        diff.to_csv(output_prefix + "diff_results.csv")

        fig = plot_differential_expression(diff)
        fig.savefig(output_prefix + "diff_results.ma_plots.svg")

        # # differential enrichment
        enrichments = differential_enrichment(
            diff,
            diff["group"].unique(),
            ["group"],
            alpha=0.05,
            alpha_col="pvals",
            max_n=100,
            sort_by="scores",
        )
        enrichments.to_csv(
            output_prefix + "enrichment_results.csv", index=False
        )
        enrichments = pd.read_csv(
            output_prefix + "enrichment_results.csv", index_col=0
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

        plot_differential_enrichment(enrichments, output_prefix, robust=True)


def label_clusters(
    sample, n_clusters=4, attribute="leiden", repl=REPL, return_labels=False
):
    obs = pd.read_csv(pjoin(run, sample, sample) + ".filtered.obs.csv")
    enrichments = pd.read_csv(
        pjoin(run, sample, sample)
        + f".{attribute}.differential_enrichment.csv",
        index_col=0,
    ).reset_index()
    enrichments = enrichments.loc[enrichments["group"] < n_clusters]
    enrichments = enrichments.query("gene_set_library == 'ARCHS4_Cell-lines'")
    labels = (
        enrichments.loc[enrichments.groupby("group")["combined_score"].idxmax()]
        .set_index("description")["group"]
        .to_dict()
    )
    repl = {k: k for k, v in labels.items()} if repl is None else repl
    labels = {repl[k]: v for k, v in labels.items()}
    if return_labels:
        return labels
    groups = {
        cell: obs.query(f"leiden == '{group}'").index.tolist()
        for cell, group in labels.items()
    }
    return groups


def crop_seq_analysis():
    crop_samples = [s for s in samples if "CRISPR" in s]

    crsig = pd.read_csv(
        "../archive/crop-seq/results/differential_expression.CROP-seq_Jurkat_TCR_stimulation.allcells.group_means.signature_cells_annotated.csv"
    )
    crsig = crsig.columns[2:-2]

    # At single cell
    for sample in crop_samples:
        # sample = crop_samples[-1]

        output_dir = pjoin(run, sample)  # + ".raw.")
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = pjoin(output_dir, sample)
        a = sc.read(output_prefix + ".filtered.h5ad")

        # # Look into activation
        # sc.tl.rank_genes_groups(a, groupby='TCR_status', n_genes=a.shape[1])

        # fc = pd.DataFrame.from_records(aa.uns['rank_genes_groups']['logfoldchanges'], index=a.var.index)
        # fc.columns = "logfc_" + fc.columns
        # p = pd.DataFrame.from_records(aa.uns['rank_genes_groups']['pvals'], index=a.var.index)
        # p.columns = "pval_" + p.columns
        # res = fc.join(p)
        # res['logfc_change'] = res['logfc_stimulated'] - res['logfc_unstimulated']
        # res = res.sort_values("logfc_change")

        #

        sc.pp.regress_out(a, keys=["n_counts"], n_jobs=4)
        sc.pp.scale(a, max_value=10)
        sc.pp.pca(a)
        sc.pp.neighbors(a)
        sc.tl.umap(a)
        sc.tl.leiden(a, resolution=0.35)
        a.write(output_prefix + ".filtered.regressout.h5ad")
        a = sc.read(output_prefix + ".filtered.regressout.h5ad")

        a = a[a.obs.sample(frac=1).index, :]

        a.obs["log_counts"] = np.log10(a.obs["n_counts"])
        fig = sc.pl.umap(
            a,
            color=["log_counts", "TCR_status", "efficiency_ratio", "gRNA"],
            return_fig=True,
            palette="tab20c",
        )
        a.obs["pos_control"] = (
            a.obs["gRNA"].isin(["LAT", "LCK", "ZAP70"]).astype(str)
        )
        fig = sc.pl.umap(
            a,
            color=[
                "log_counts",
                "TCR_status",
                "efficiency_ratio",
                "gRNA",
                "pos_control",
            ],
            return_fig=True,
            palette="Paired",
        )
        fig = sc.pl.umap(
            a,
            color=["leiden"],
            return_fig=True,
            palette="tab20b",
        )
        [
            c.set_rasterized(True)
            for c in fig.axes[0].get_children()
            if isinstance(c, matplotlib.collections.PathCollection)
        ]
        fig.savefig(
            output_prefix + ".filtered.regressout.umap.leiden_clusters.svg",
            dpi=300,
            bbox_inches="tight",
        )

        # gRNA overwrepresentation
        import scipy
        from statsmodels.stats.multitest import multipletests

        t = -np.log10(0.01)
        for group, state in [("TCR_status", "stimulated"), ("leiden", "0")]:
            for variable in ["gRNA_ID", "gRNA"]:
                grnas = a.obs[variable].unique()
                table = pd.DataFrame(index=grnas, columns=["a", "b", "c", "d"])
                for i, grna in enumerate(grnas):
                    # KOs in stimulated
                    table.loc[grna, "a"] = a.obs.query(
                        f"{group}  == '{state}' & {variable} == '{grna}'"
                    ).shape[0]
                    # KO without stimulation
                    table.loc[grna, "b"] = a.obs.query(
                        f"{group}  != '{state}' & {variable} == '{grna}'"
                    ).shape[0]
                    # unstimulated with KO
                    table.loc[grna, "c"] = a.obs.query(
                        f"{group}  == '{state}' & {variable} != '{grna}'"
                    ).shape[0]
                    table.loc[grna, "d"] = a.obs.query(
                        f"{group}  != '{state}' & {variable} != '{grna}'"
                    ).shape[0]
                    odds, p = scipy.stats.fisher_exact(
                        table.loc[grna].dropna().values.reshape((2, 2)),
                        alternative="two-sided",
                    )
                    table.loc[grna, "odds"] = odds
                    table.loc[grna, "p"] = p

                table["fdr"] = multipletests(table["p"], method="fdr_bh")[1]
                table["log_odds"] = np.log(table["odds"])
                table["logp"] = -np.log10(table["p"])
                table["log_fdr"] = -np.log10(table["fdr"])
                table.to_csv(
                    output_prefix
                    + f".CRISPR.single_cell.{group}.{variable}.fisher_results.csv"
                )
                table = table.replace(np.inf, np.nan)

                fig, axes = plt.subplots(1, figsize=(4, 4))
                axes.scatter(data=table, x="log_odds", y="logp")
                axes.axhline(
                    table.query(f"log_fdr >= {t}")["log_fdr"].min(),
                    linestyle="--",
                    color="grey",
                )
                axes.axvline(0, linestyle="--", color="grey")
                v = table["log_odds"].abs().max()
                v += v * 0.1
                axes.set_xlim(-v, v)
                for grna, row in table.query(f"log_fdr >= {t}").iterrows():
                    axes.text(row["log_odds"], row["logp"], grna)
                axes.set_xlabel("Log(Odds-ratio)")
                axes.set_ylabel("-log10(p-value)")
                axes.set_title(
                    "Enrichment of genes in Stimulated/Unstimulated cells"
                )
                fig.savefig(
                    output_prefix
                    + f".CRISPR.single_cell.{group}.{variable}.fisher_results.volcano.svg",
                    bbox_inches="tight",
                    dpi=300,
                )

                # 20200831 version: only target cells (no background)
                # illustrate
                # group = "leiden"  # TCR_status
                # state = "0"  # "stimulated"
                # variable = "gRNA_ID"
                table = pd.read_csv(
                    output_prefix
                    + f".CRISPR.single_cell.{group}.{variable}.fisher_results.csv",
                    index_col=0,
                )
                # sigs = table.query(f"log_fdr >= {t}")
                sigs = table

                fig, axes = plt.subplots(6, 8, figsize=(8 * 4, 6 * 4))
                axes = axes.flatten()
                colors = sns.color_palette("colorblind")
                stim = a.obs[group] == state
                for i in range(sigs.shape[0]):
                    axes[i].scatter(
                        *a.obsm["X_umap"][~stim, :].T,
                        color=colors[0],
                        s=2,
                        alpha=0.05,
                        rasterized=True,
                    )
                    axes[i].scatter(
                        *a.obsm["X_umap"][stim, :].T,
                        color=colors[1],
                        s=2,
                        alpha=0.05,
                        rasterized=True,
                    )

                for i, (grna, row) in enumerate(
                    sigs.sort_values("log_odds").iterrows()
                ):
                    gmask = a.obs[variable] == grna
                    axes[i].scatter(
                        *a.obsm["X_umap"][gmask, :].T,
                        color="black",
                        s=5,
                        alpha=0.5,
                        rasterized=True,
                    )
                    axes[i].set_title(
                        f"{grna}\n"
                        f"Log-odds: {table.loc[grna, 'log_odds']:.2f}\n"
                        f"FDR: {table.loc[grna, 'fdr']:.2e}"
                    )
                for ax in axes.flat:
                    ax.axis("off")

                fig.savefig(
                    output_prefix
                    + f".CRISPR.single_cell.{group}.{variable}.fisher_results.umap_demonstration.single_gRNAs.svg",
                    bbox_inches="tight",
                    dpi=300,
                )

        # Compare over-representation of gRNAs in stimulated state vs UMAP space
        fig, axes = plt.subplots(1, 2, figsize=(4 * 2, 4), tight_layout=True)
        for i, variable in enumerate(["gRNA_ID", "gRNA"]):
            tcr = pd.read_csv(
                output_prefix
                + f".CRISPR.single_cell.TCR_status.{variable}.fisher_results.csv",
                index_col=0,
            )

            leiden = pd.read_csv(
                output_prefix
                + f".CRISPR.single_cell.leiden.{variable}.fisher_results.csv",
                index_col=0,
            )
            leiden = leiden.reindex(tcr.index)

            _a = tcr["log_odds"]
            _b = leiden["log_odds"]
            c = _b - _a
            v = pd.DataFrame([_a, _b]).T.abs().max().max()
            v2 = v + v * 0.1
            axes[i].plot((-v2, v2), (-v2, v2), linestyle="--", color="grey")
            axes[i].scatter(_a, _b, c=c, cmap="coolwarm", vmin=-v, vmax=v)
            axes[i].set(
                xlabel="Effect on survival\n(Cell depletion)",
                ylabel="Effect on TCR activation\n(Transcriptome change)",
                xlim=(-v2, v2),
                ylim=(-v2, v2),
                title="Gene" if i == 1 else "gRNA",
            )
            done = list()
            n = 8
            for t in c.sort_values().head(n).index:
                axes[i].text(_a.loc[t], _b.loc[t], s=t)
                done.append(t)
            for t in c.sort_values().tail(n).index:
                if t not in done:
                    axes[i].text(_a.loc[t], _b.loc[t], s=t)
                    done.append(t)
            for t in _a.sort_values().head(n).index:
                if t not in done:
                    axes[i].text(_a.loc[t], _b.loc[t], s=t)
                    done.append(t)
            for t in _a.sort_values().tail(n).index:
                if t not in done:
                    axes[i].text(_a.loc[t], _b.loc[t], s=t)
                    done.append(t)
            for t in _b.sort_values().head(n).index:
                if t not in done:
                    axes[i].text(_a.loc[t], _b.loc[t], s=t)
                    done.append(t)
            for t in _b.sort_values().tail(n).index:
                if t not in done:
                    axes[i].text(_a.loc[t], _b.loc[t], s=t)
                    done.append(t)
        fig.savefig(
            output_prefix
            + f".CRISPR.single_cell.survival_effect_vs_activation_effect.svg",
            bbox_inches="tight",
            dpi=300,
        )

        # # Use coefficients
        # sc.tl.rank_genes_groups(a, groupby="gRNA_ID", n_genes=a.shape[1])
        # c = pd.DataFrame.from_records(a.uns["rank_genes_groups"]["pvals_adj"], index=a.var.index)

        # m = c.abs().max(1)
        # m = c.loc[m >= m.quantile(0.9)]

        # sns.clustermap(m, center=0, metric="correlation", cmap="RdBu_r", robust=True)

        # sig = a2.var.loc[a2.var["external_gene_name"].isin(crsig)].index

        # a[:, sig]
        # import diffxpy.api as de

        # wald = de.test.wald(
        #     data=a.raw.to_adata()[:, a.var.index],
        #     formula_loc="~ 1 + TCR_status",
        #     factor_loc_totest="TCR_status",
        # )
        # walds = test.summary()
        # walds.to_csv(output_prefix + f".CRISPR.single_cell.{variable}.diffxpy.TCR_status.wald.csv",)
        # ttest = de.test.t_test(data=a.raw.to_adata()[:, a.var.index], grouping="TCR_status")
        # ttests = test.summary()
        # ttests = ttests.query("abs(log2fc) < 300")
        # ttests.to_csv(
        #     output_prefix + f".CRISPR.single_cell.{variable}.diffxpy.TCR_status.ttest.csv",
        # )

        # plt.scatter(np.log1p(ttests["mean"]), ttests["log2fc"], c=-np.log10(ttests["pval"]))

        # araw = a.raw.to_adata()[:, a.var.index]
        # araw = a
        # ttests = list()
        # for state in ["stimulated", "unstimulated"]:
        #     braw = araw[(araw.obs["TCR_status"] == state), :]
        #     for grna in a.obs["gRNA_ID"].unique():
        #         print(state, grna)
        #         braw.obs["test"] = braw.obs["gRNA_ID"] == grna
        #         ttest = de.test.t_test(data=braw, grouping="test")
        #         ttests.append(ttest.summary().assign(TCR_status=state, gRNA=grna))
        # ttests = pd.concat(ttests)
        # ttests = ttests.query("abs(log2fc) < 300")
        # ttests.to_csv(
        #     output_prefix + f".CRISPR.single_cell.{variable}.diffxpy.TCR_status+gRNA.ttest.csv",
        # )

        # p = ttests.pivot_table(
        #     index="gene", columns=["TCR_status", "gRNA"], values="log2fc", fill_value=0
        # )

        # sns.clustermap(
        #     p.corr(),
        #     mask=np.eye(p.shape[1], dtype=bool),
        #     center=0,
        #     cmap="RdBu_r",
        #     metric="correlation",
        #     xticklabels=True,
        #     yticklabels=True,
        # )

        # q = p.T
        # sns.clustermap(
        #     q.corr(),
        #     mask=np.eye(q.shape[1], dtype=bool),
        #     center=0,
        #     cmap="RdBu_r",
        #     metric="correlation",
        # )

        # # a[:, sig]
        # import diffxpy.api as de

        # wald = de.test.wald(
        #     data=a.raw.to_adata()[:, a.var.index],
        #     formula_loc="~ 1 + gRNA + TCR_status",
        #     factor_loc_totest="gRNA",
        # )

    # Aggregated
    min_cells = 25
    cd69 = "ENSG00000110848"
    for sample in crop_samples:
        # sample = crop_samples[-1]

        output_dir = pjoin(run, sample)  # + ".raw.")
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = pjoin(output_dir, sample)
        a = sc.read(output_prefix + ".filtered.h5ad")

        for label, agg_vars in [
            ("TCR_status+gene", ["TCR_status", "gRNA"]),
            ("TCR_status+gRNA", ["TCR_status", "gRNA", "gRNA_ID"]),
        ]:
            # a2 = AnnData(a.raw.to_adata().to_df().join(a.obs[agg_vars]).groupby(agg_vars).mean())
            # get number of cells per group
            n_cells = (
                a.obs.groupby(agg_vars)
                .count()
                .iloc[:, 0]
                .rename("Cells per group")
            )
            n_cells = n_cells[n_cells >= min_cells]
            matrix = a.to_df().join(a.obs[agg_vars]).groupby(agg_vars).mean()

            matrix = matrix.loc[n_cells.index]
            matrix += abs(matrix.min().min())

            a2 = AnnData(matrix)
            a2.var = a2.var.join(a.var)
            a2.obs = a2.obs.index.to_frame()
            a2.obs.index = pd.MultiIndex.from_frame(
                a2.obs.index.to_series().apply(pd.Series).set_axis(agg_vars, 1)
            )
            a2._obs = a2.obs.join(n_cells)

            sc.pp.normalize_per_cell(
                a2, counts_per_cell_after=1e4, min_counts=0
            )
            sc.pp.regress_out(a2, keys=["n_counts"], n_jobs=4)

            # sc.pp.log1p(a2)
            sc.pp.scale(a2, max_value=10)
            sc.pp.pca(a2, zero_center=0)
            sc.pp.neighbors(a2)
            sc.tl.diffmap(a2)
            sc.tl.umap(a2)
            sc.tl.leiden(a2, resolution=0.8)

            fig = sc.pl.pca(a2, color=a2.obs.columns, return_fig=True)
            for variable, ax in zip(a2.obs.columns, fig.axes):
                pca = a2.obsm["X_pca"]
                for i, lab in enumerate(a2.obs.index):
                    ax.text(
                        pca[i, 0],
                        pca[i, 1],
                        s=lab[2 if variable == "gRNA_ID" else 1],
                    )
            fig.savefig(output_prefix + f".CRISPR.groupby_{label}.PCA.svg")

            fig = sc.pl.diffmap(a2, color=a2.obs.columns, return_fig=True)
            for variable, ax in zip(a2.obs.columns, fig.axes):
                dm = a2.obsm["X_diffmap"]
                for i, lab in enumerate(a2.obs.index):
                    ax.text(dm[i, 1], dm[i, 2], s=lab[1])
            fig.savefig(output_prefix + f".CRISPR.groupby_{label}.diffmap.svg")

            fig = sc.pl.umap(a2, color=a2.obs.columns, return_fig=True)
            for variable, ax in zip(a2.obs.columns, fig.axes):
                umap = a2.obsm["X_umap"]
                for i, lab in enumerate(a2.obs.index):
                    ax.text(umap[i, 0], umap[i, 1], s=lab[1])
            fig.savefig(output_prefix + f".CRISPR.groupby_{label}.umap.svg")

            # pairwise correlation
            c = a2.to_df().T.corr()
            grid = sns.clustermap(
                c,
                xticklabels=True,
                yticklabels=True,
                robust=True,
                metric="correlation",
                rasterized=True,
                center=0,
                cmap="RdBu_r",
                row_colors=n_cells,
                col_colors=n_cells,
            )
            grid.savefig(
                output_prefix
                + f".CRISPR.groupby_{label}.pairwise_correlation.clustermap.svg"
            )

            # Gene subsets

            # # all
            # # highly variable
            hvg = sc.pp.highly_variable_genes(a2, inplace=False)
            hvg.index = a2.var.index
            hvgs = hvg.query("highly_variable").index
            # # differential
            sc.tl.rank_genes_groups(a2, groupby="TCR_status", n_genes=1500)
            degs = np.asarray(
                a2.uns["rank_genes_groups"]["names"].tolist()
            ).flatten()

            # # crop-seq signature
            sig = a2.var.loc[a2.var["external_gene_name"].isin(crsig)].index

            for label2, genes in [
                ("all_genes", a2.var.index),
                ("highly_variable_genes", hvgs),
                ("differential_genes", degs),
                ("original_cropseq_signature", sig),
            ]:
                diff = pd.DataFrame(a2.uns["rank_genes_groups"]["pvals"])[
                    "stimulated"
                ]
                diff.index = pd.DataFrame(a2.uns["rank_genes_groups"]["names"])[
                    "stimulated"
                ]
                diff = diff.sort_values()
                for g in diff.head(50).index:
                    gname = a.var.loc[g, "external_gene_name"]
                    gs = pd.Series(
                        a2[:, g].X.squeeze(), index=n_cells.index, name=gname
                    )
                    p = a2.to_df()[genes].T
                    grid = sns.clustermap(
                        p,
                        metric="correlation",
                        center=0,
                        cmap="RdBu_r",
                        robust=True,
                        z_score=1,
                        xticklabels=True,
                        yticklabels=False,
                        rasterized=True,
                        row_colors=p.mean(1).rename("Mean expression"),
                        col_colors=n_cells.to_frame().join(gs),
                    )
                    grid.savefig(
                        output_prefix
                        + f".CRISPR.groupby_{label}.{label2}.clustermap.{gname}_expression.svg"
                    )

                # # # sort by difference to respective CTRL

                _order = dict()
                for grna_state in p.columns.levels[0][::-1]:
                    for ctrl_state in p.columns.levels[0][::-1]:
                        for y1, y2 in zip(
                            p.columns.get_level_values(1),
                            p.columns.get_level_values(2),
                        ):
                            try:
                                _order[
                                    (grna_state, ctrl_state, y1, y2)
                                ] = np.sqrt(
                                    sum(
                                        abs(
                                            p[(grna_state, y1, y2)]
                                            - p[ctrl_state, "CTRL"].mean(1)
                                        )
                                    )
                                )
                            except KeyError:
                                pass
                # Make into DataFrame
                order = (
                    (
                        pd.Series(_order).rename_axis(
                            index=[
                                "TCR_status",
                                "ctrl_state",
                                "gRNA",
                                "gRNA_ID",
                            ]
                        )
                    )
                    .to_frame()
                    .pivot_table(
                        index=["TCR_status", "gRNA", "gRNA_ID"],
                        columns="ctrl_state",
                        values=0,
                    )
                )
                # # Replace Zero of controls with minimum distance
                # order = order.replace(0, order.replace(0, np.nan).min().min())
                # Calculate difference distances to controls
                # order = order - order.min()
                order = (order.T - order.mean(1)).T
                order = order - order.mean()
                order["stimulation"] = (
                    order["unstimulated"] - order["stimulated"]
                )
                order = order.sort_values("stimulation")

                order.to_csv(
                    output_prefix
                    + f".CRISPR.groupby_{label}.{label2}.deviation_from_ctrl.csv"
                )

                for g in diff.head(50).index:
                    gname = a.var.loc[g, "external_gene_name"]
                    gs = pd.Series(
                        a2[:, g].X.squeeze(), index=n_cells.index, name=gname
                    )

                    grid = sns.clustermap(
                        p.loc[:, order.index],
                        metric="correlation",
                        center=0,
                        cmap="RdBu_r",
                        robust=True,
                        z_score=1,
                        col_cluster=False,
                        xticklabels=True,
                        yticklabels=False,
                        rasterized=True,
                        col_colors=order.join(gs),
                        cbar_kws=dict(label="Z-score"),
                    )
                    grid.savefig(
                        output_prefix
                        + f".CRISPR.groupby_{label}.{label2}.clustermap.ordered.{gname}_expression.svg"
                    )

                fig, ax = plt.subplots(1, 1, figsize=(10, 2))
                ax.set(xlabel="TCR activation (rank)", ylabel="TCR activation")
                ax.axhline(0, linestyle="--", color="grey")
                rank = order["stimulation"].rank()
                ax.scatter(
                    rank,
                    order["stimulation"],
                    c=order["stimulation"],
                    cmap="coolwarm",
                    vmin=-10,
                    vmax=10,
                )
                for t in order.loc[
                    order.index.get_level_values(0) == "unstimulated"
                ].index:
                    ax.text(
                        rank.loc[t],
                        order.loc[t, "stimulation"],
                        s=t[-1],
                        rotation=90,
                        va="bottom",
                    )
                for t in order.loc[
                    order.index.get_level_values(0) == "stimulated"
                ].index:
                    ax.text(
                        rank.loc[t],
                        order.loc[t, "stimulation"],
                        s=t[-1],
                        rotation=90,
                        va="top",
                    )
                fig.savefig(
                    output_prefix
                    + f".CRISPR.groupby_{label}.{label2}.activation.rank_vs_score.scatter.svg",
                    dpi=300,
                    bbox_inches="tight",
                )

                # Groupby gene
                pagg = p.T.groupby(level=[0, 1]).mean().T
                orderagg = (
                    order.groupby(level=[0, 1])
                    .mean()
                    .sort_values("stimulation")
                )

                orderagg.to_csv(
                    output_prefix
                    + f".CRISPR.groupby_{label}.{label2}.deviation_from_ctrl.aggregated_by_gene.csv"
                )

                grid = sns.clustermap(
                    pagg.loc[:, orderagg.index],
                    metric="correlation",
                    center=0,
                    cmap="RdBu_r",
                    robust=True,
                    z_score=1,
                    col_cluster=False,
                    xticklabels=True,
                    yticklabels=False,
                    rasterized=True,
                    col_colors=orderagg,
                    cbar_kws=dict(label="Z-score"),
                )
                grid.savefig(
                    output_prefix
                    + f".CRISPR.groupby_{label}.{label2}.clustermap.ordered.aggregated_by_gene.svg"
                )

                fig, ax = plt.subplots(1, 1, figsize=(10, 2))
                ax.set(xlabel="TCR activation (rank)", ylabel="TCR activation")
                ax.axhline(0, linestyle="--", color="grey")
                rank = orderagg["stimulation"].rank()
                ax.scatter(
                    rank,
                    orderagg["stimulation"],
                    c=orderagg["stimulation"],
                    cmap="coolwarm",
                    vmin=-10,
                    vmax=10,
                )
                for t in orderagg.loc[
                    orderagg.index.get_level_values(0) == "unstimulated"
                ].index:
                    ax.text(
                        rank.loc[t],
                        orderagg.loc[t, "stimulation"],
                        s=t[-1],
                        rotation=90,
                        va="bottom",
                    )
                for t in orderagg.loc[
                    orderagg.index.get_level_values(0) == "stimulated"
                ].index:
                    ax.text(
                        rank.loc[t],
                        orderagg.loc[t, "stimulation"],
                        s=t[-1],
                        rotation=90,
                        va="top",
                    )
                fig.savefig(
                    output_prefix
                    + f".CRISPR.groupby_{label}.{label2}.activation.rank_vs_score.aggregated_by_gene.scatter.svg",
                    dpi=300,
                    bbox_inches="tight",
                )


# The stuff below is to patch seaborn clustermap
from functools import wraps
from typing import Optional, Union, List, Literal, Callable

Series = Union[pd.Series]
DataFrame = Union[pd.DataFrame]

SEQUENCIAL_CMAPS = [
    "Purples",
    "Greens",
    "Oranges",
    "Greys",
    "Reds",
    "Blues",
    "YlOrBr",
    "YlOrRd",
    "OrRd",
    "PuRd",
    "RdPu",
    "BuPu",
    "GnBu",
    "PuBu",
    "YlGnBu",
    "PuBuGn",
    "BuGn",
    "YlGn",
    "binary",
    "gist_yarg",
    "gist_gray",
    "gray",
    "bone",
    "pink",
    "spring",
    "summer",
    "autumn",
    "winter",
    "cool",
    "Wistia",
    "hot",
    "afmhot",
    "gist_heat",
    "copper",
]


def minmax_scale(x):
    return (x - x.min()) / (x.max() - x.min())


def to_color_series(x: Series, cmap: Optional[str] = "Greens") -> Series:
    """Map a numeric pandas series to a series of RBG values."""
    return Series(
        plt.get_cmap(cmap)(minmax_scale(x)).tolist(), index=x.index, name=x.name
    )


def to_color_dataframe(
    x: Union[Series, DataFrame],
    cmaps: Optional[Union[str, List[str]]] = None,
    offset: int = 0,
) -> DataFrame:
    """Map a numeric pandas DataFrame to RGB values."""
    if isinstance(x, pd.Series):
        x = x.to_frame()
    if cmaps is None:
        # the offset is in order to get different colors for rows and columns by default
        cmaps = [plt.get_cmap(cmap) for cmap in SEQUENCIAL_CMAPS[offset:]]
    if isinstance(cmaps, str):
        cmaps = [cmaps]
    return pd.concat(
        [to_color_series(x[col], cmap) for col, cmap in zip(x, cmaps)], axis=1
    )


def _add_extra_colorbars_to_clustermap(
    grid: sns.matrix.ClusterGrid,
    datas: Union[Series, DataFrame],
    cmaps: Optional[Union[str, List[str]]] = None,
    location: str = Union[Literal["col"], Literal["row"]],
) -> None:
    """Add either a row or column colorbar to a seaborn Grid."""

    def add(
        data: Series, cmap: str, bbox: List[List[int]], orientation: str
    ) -> None:
        ax = grid.fig.add_axes(matplotlib.transforms.Bbox(bbox))
        norm = matplotlib.colors.Normalize(vmin=data.min(), vmax=data.max())
        cb1 = matplotlib.colorbar.ColorbarBase(
            ax,
            cmap=plt.get_cmap(cmap),
            norm=norm,
            orientation=orientation,
            label=data.name,
        )

    offset = 1 if location == "row" else 0

    if isinstance(datas, pd.Series):
        datas = datas.to_frame()
    if cmaps is None:
        cmaps = SEQUENCIAL_CMAPS[offset:]
    if isinstance(cmaps, str):
        cmaps = [cmaps]

    # get position to add new axis in existing figure
    # # get_position() returns ((x0, y0), (x1, y1))
    heat = grid.ax_heatmap.get_position()
    cbar_spacing = 0.05
    cbar_size = 0.025
    if location == "col":
        orientation = "vertical"
        dend = grid.ax_col_dendrogram.get_position()
        y0 = dend.y0
        y1 = dend.y1
        for i, (data, cmap) in enumerate(zip(datas, cmaps)):
            if i == 0:
                x0 = heat.x1
                x1 = heat.x1 + cbar_size
            else:
                x0 += cbar_size + cbar_spacing
                x1 += cbar_size + cbar_spacing
            add(datas[data], cmap, [[x0, y0], [x1, y1]], orientation)
    else:
        orientation = "horizontal"
        dend = grid.ax_row_dendrogram.get_position()
        x0 = dend.x0
        x1 = dend.x1
        for i, (data, cmap) in enumerate(zip(datas, cmaps)):
            if i == 0:
                y0 = dend.y0 - cbar_size
                y1 = dend.y0
            else:
                y0 -= cbar_size + cbar_spacing
                y1 -= cbar_size + cbar_spacing
            add(datas[data], cmap, [[x0, y0], [x1, y1]], orientation)


def _add_colorbars(
    grid: sns.matrix.ClusterGrid,
    rows: DataFrame = None,
    cols: DataFrame = None,
    row_cmaps: Optional[List[str]] = None,
    col_cmaps: Optional[List[str]] = None,
) -> None:
    """Add row and column colorbars to a seaborn Grid."""
    if rows is not None:
        _add_extra_colorbars_to_clustermap(
            grid, rows, location="row", cmaps=row_cmaps
        )
    if cols is not None:
        _add_extra_colorbars_to_clustermap(
            grid, cols, location="col", cmaps=col_cmaps
        )


def colorbar_decorator(f: Callable) -> Callable:
    """
    Decorate seaborn.clustermap in order to have numeric values passed to the
    ``row_colors`` and ``col_colors`` arguments translated into row and column
    annotations and in addition colorbars for the restpective values.
    """
    # TODO: edit original seaborn.clustermap docstring to document {row,col}_colors_cmaps arguments.
    @wraps(f)
    def clustermap(*args, **kwargs):
        cmaps = {"row": None, "col": None}
        # capture "row_cmaps" and "col_cmaps" out of the kwargs
        for arg in ["row", "col"]:
            if arg + "_colors_cmaps" in kwargs:
                cmaps[arg] = kwargs[arg + "_colors_cmaps"]
                del kwargs[arg + "_colors_cmaps"]
        # get dataframe with colors and respective colormaps for rows and cols
        # instead of the original numerical values
        _kwargs = dict(rows=None, cols=None)
        for arg in ["row", "col"]:
            if arg + "_colors" in kwargs:
                if isinstance(
                    kwargs[arg + "_colors"], (pd.DataFrame, pd.Series)
                ):
                    _kwargs[arg + "s"] = kwargs[arg + "_colors"]
                    kwargs[arg + "_colors"] = to_color_dataframe(
                        x=kwargs[arg + "_colors"],
                        cmaps=cmaps[arg],
                        offset=1 if arg == "row" else 0,
                    )
        grid = f(*args, **kwargs)
        _add_colorbars(
            grid, **_kwargs, row_cmaps=cmaps["row"], col_cmaps=cmaps["col"]
        )
        return grid

    return clustermap


sns.clustermap = colorbar_decorator(sns.clustermap)
