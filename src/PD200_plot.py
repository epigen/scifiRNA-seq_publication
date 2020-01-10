#!/usr/bin/env python


import os
from functools import partial

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np

sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
matplotlib.rcParams["svg.fonttype"] = "none"
# Don't use LaTeX for rendering
matplotlib.rcParams["text.usetex"] = False


def inflection_point(curve):
    """Return the index of the inflection point of a curve"""
    # knee
    from numpy.matlib import repmat

    n_points = len(curve)
    all_coord = np.vstack((range(n_points), curve)).T
    line_vec = (all_coord[-1] - all_coord[0])
    line_vec_norm = line_vec / np.sqrt(np.sum(line_vec ** 2))
    vec_from_first = all_coord - all_coord[0]
    scalar_product = np.sum(
        vec_from_first * repmat(line_vec_norm, n_points, 1), axis=1)
    vec_to_line = vec_from_first - np.outer(scalar_product, line_vec_norm)
    return np.argmax(np.sqrt(np.sum(vec_to_line ** 2, axis=1)))


annot = pd.read_csv(os.path.join("metadata", "annotation.csv"))
experiments = annot.query(
    "sample_name.str.contains('PD200').values", parser="python"
)

# Collect metrics
metrics_a = dict()
metrics_f = dict()
inflex = dict()
for experiment in experiments["sample_name"]:
    print(experiment)
    # read metrics
    m = pd.read_csv(
        os.path.join("data", experiment, experiment + ".metrics.csv.gz")
    ).sort_values("umi")
    m["umi_read_ratio"] = m["umi"] / m["read"]
    m["gene_umi_ratio"] = m["gene"] / m["umi"]
    if "material_type" not in m.columns:
        m.loc[:, 'material_type'] = experiments.loc[
            experiments['sample_name'] == experiment, 'material'].squeeze()
    inflex[experiment] = inflection_point(m['umi'])
    metrics_a[experiment] = m.query("umi > 10")
    metrics_f[experiment] = m.iloc[inflex[experiment]:]

metrics_a = (
    pd.concat(metrics_a)
    .reset_index(level=1, drop=True)
    .reset_index()
    .rename(columns={"index": "sample"})
    .sort_values(["sample", "umi"])
)
metrics_f = (
    pd.concat(metrics_f)
    .reset_index(level=1, drop=True)
    .reset_index()
    .rename(columns={"index": "sample"})
    .sort_values(["sample", "umi"])
)
# metrics.to_csv("results/PD200.metrics.join_experiments.csv.gz")

fig_kws = dict(dpi=300, bbox_inches="tight")

# quick summary
vars_ = [
    "read",
    "umi",
    "gene",
    "unique_fraction",
    "umi_read_ratio",
    "gene_umi_ratio",
    "doublet",
]

p75 = partial(np.percentile, q=75)
p75.__name__ = "75_percentile"
p25 = partial(np.percentile, q=25)
p25.__name__ = "25_percentile"

summary = (
    metrics_f.query("umi > 100")
    # .nlargest(20000, columns="umi")
    .groupby(["sample", "material_type"])[vars_]
    .agg([p25, np.mean, np.median, p75])
    .T
)
summary.to_csv("results/PD200/PD200_comparison.summary.csv")

# summary
for label in ["per_sample", "per_material"]:
    fig, axis = plt.subplots(1, 7, figsize=(3 * 7, 3 * 1))
    ax = iter(axis.flatten())
    b_kws = dict(whis=1e9)
    if label == "per_sample":
        b_kws.update(
            dict(
                x="sample",
                hue="material_type",
                order=metrics_f["sample"].unique(),
                hue_order=metrics_f["material_type"].unique(),
            )
        )
    else:
        b_kws.update(
            dict(
                x="material_type",
                hue="sample",
                order=metrics_f["material_type"].unique(),
                hue_order=metrics_f["sample"].unique(),
            )
        )
    for var_ in vars_[:-1]:
        sns.boxplot(data=metrics_f, y=var_, ax=next(ax), **b_kws)
    del b_kws["whis"]
    sns.violinplot(data=metrics_f, y=vars_[-1], ax=next(ax), **b_kws)
    for ax in axis.flatten()[:3]:
        ax.set_yscale("log")
    for t, ax in zip(vars_, axis.flatten()):
        ax.set_title(t)
        ax.set_ylabel(ax.get_ylabel(), visible=False)
    for ax in axis.flatten()[1:]:
        ax.get_legend().set_visible(False)
    if label == "per_sample":
        for ax in axis.flatten():
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.savefig(
        f"results/PD200/PD200_comparison.summary_boxplot.{label}.svg",
        **fig_kws,
    )

# Illustrations
s = metrics_f[['sample', 'material_type']].drop_duplicates()
n = s.shape[0]

fig, axis = plt.subplots(2, 5, figsize=(3 * 5, 3 * 2))
for i, (sample, material) in s.iterrows():
    m_a = metrics_a.query(
        f"(sample == '{sample}') & (umi > 100) & (material_type == '{material}')"
    )
    m_f = metrics_f.query(
        f"(sample == '{sample}') & (material_type == '{material}')"
    )
    for ax, m in zip(axis, (m_a, m_f)):
        kwds = dict(label=sample + " - " + material)
        s_kwds = dict(alpha=0.1, s=1, rasterized=True)
        rr = m["read"].rank(ascending=False)
        ur = m["umi"].rank(ascending=False)
        ax[0].plot(rr, m["read"], **kwds)
        ax[0].loglog()
        ax[0].set_xlabel("Barcodes")
        ax[0].set_ylabel("Reads")
        ax[1].plot(ur, m["umi"], **kwds)
        ax[1].loglog()
        ax[1].set_xlabel("Barcodes")
        ax[1].set_ylabel("UMIs")
        ax[2].scatter(m["read"], m["umi"], **s_kwds, **kwds)
        ax[2].loglog()
        ax[2].set_xlabel("Reads")
        ax[2].set_ylabel("UMIs")
        ax[3].scatter(m["umi"], m["gene"], **s_kwds, **kwds)
        ax[3].loglog()
        ax[3].set_xlabel("UMIs")
        ax[3].set_ylabel("Genes")
        ax[4].scatter(m["umi"], m["unique_fraction"], **s_kwds, **kwds)
        ax[4].set_xscale("log")
        ax[4].set_xlabel("UMIs")
        ax[4].set_ylabel("Unique fraction")

axis[0, 0].legend()
axis[0, 0].set_ylabel(
    "All barcodes\n" + axis[0, 0].get_ylabel(), ha="center", va="bottom"
)
axis[1, 0].set_ylabel(
    "Real cells\n" + axis[1, 0].get_ylabel(), ha="center", va="bottom"
)
fig.savefig(
    "results/PD200/PD200_comparison.joint_comparison.svg",
    **fig_kws,
)

# Species mixing
for label in ["log", "linear"]:
    for fixed_axis, l2 in [(False, "free_scale"), (True, "fixed")]:
        fig, axis = plt.subplots(
            1, n, figsize=(3 * n, 3 * 1), sharex=fixed_axis, sharey=fixed_axis
        )
        axis = iter(axis.flatten())
        for i, (sample, material) in s.iterrows():
            m_a = metrics_f.query(
                f"(sample == '{sample}') & (material_type == '{material}')"
            )
            m = m_a.nlargest(10000, columns="umi")
            kwds = dict(label=sample + " - " + material)
            s_kwds = dict(
                alpha=0.1 if label == "log" else 0.8, s=1, rasterized=True
            )
            ax = next(axis)
            ax.set_title(kwds["label"])
            ax.scatter(m["mouse"], m["human"], **s_kwds, **kwds)
            if label == "log":
                ax.loglog()
            ax.set_xlabel("Mouse UMIs")
            ax.set_ylabel("Human UMIs")
        fig.savefig(
            f"results/PD200/PD200_comparison.species_mixing.{label}.{l2}.separate.svg",
            **fig_kws,
        )
