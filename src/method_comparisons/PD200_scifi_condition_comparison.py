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

annot = pd.read_csv(os.path.join("metadata", "annotation.csv"))
experiments = annot.query(
    "sample_name.str.contains('PD200').values", parser="python"
).query(
    "~sample_name.str.contains('PD200_10xscRNA').values", parser="python"
)
# Collect metrics
metrics_d = dict()
for experiment in experiments["sample_name"]:
    # read metrics
    metrics_d[experiment] = pd.read_csv(
        os.path.join("data", experiment, experiment + ".metrics.csv.gz")
    )
metrics = (
    pd.concat(metrics_d)
    .reset_index(level=1, drop=True)
    .reset_index()
    .rename(columns={"index": "sample"})
    .sort_values(["sample", "umi"])
)
metrics["umi_read_ratio"] = metrics["umi"] / metrics["read"]
metrics["gene_umi_ratio"] = metrics["gene"] / metrics["umi"]


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
    metrics.query("umi > 100")
    .nlargest(20000, columns="umi")
    .groupby(["sample", "material_type"])[vars_]
    .agg([p25, np.mean, np.median, p75])
    .T
)
summary.to_csv("results/PD200/PD200_scifi_condition_comparison.summary.csv")

# summary
for label in ["per_sample", "per_material"]:
    fig, axis = plt.subplots(1, 7, figsize=(3 * 7, 3 * 1))
    ax = iter(axis.flatten())
    m = metrics.query("umi > 100").nlargest(20000, columns="umi")
    b_kws = dict(whis=1e9)
    if label == "per_sample":
        b_kws.update(
            dict(
                x="sample",
                hue="material_type",
                order=metrics["sample"].unique(),
                hue_order=metrics["material_type"].unique(),
            )
        )
    else:
        b_kws.update(
            dict(
                x="material_type",
                hue="sample",
                order=metrics["material_type"].unique(),
                hue_order=metrics["sample"].unique(),
            )
        )
    for var_ in vars_[:-1]:
        sns.boxplot(data=m, y=var_, ax=next(ax), **b_kws)
    del b_kws["whis"]
    sns.violinplot(data=m, y=vars_[-1], ax=next(ax), **b_kws)
    for ax in axis.flatten()[:2]:
        ax.set_yscale("log")
    for t, ax in zip(vars_, axis.flatten()):
        ax.set_title(t)
        ax.set_ylabel(ax.get_ylabel(), visible=False)
    fig.savefig(
        f"results/PD200/PD200_scifi_condition_comparison.summary_boxplot.{label}.svg",
        **fig_kws,
    )

# Illustrations
fig, axis = plt.subplots(2, 5, figsize=(3 * 5, 3 * 2))
for e in experiments["sample_name"]:
    for material in ["cells", "nuclei"]:
        m_a = metrics.query(
            f"(sample == '{e}') & (umi > 100) & (material_type == '{material}')"
        )
        m_f = m_a.nlargest(10000, columns="umi")
        for ax, m in zip(axis, (m_a, m_f)):
            kwds = dict(label=e + " - " + material)
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
    "results/PD200/PD200_scifi_condition_comparison.joint_comparison.svg",
    **fig_kws,
)

# Species mixing
for label in ["log", "linear"]:
    fig, axis = plt.subplots(
        2, 2, figsize=(3 * 2, 3 * 2), sharex=True, sharey=True
    )
    axis = iter(axis.flatten())
    for e in experiments["sample_name"]:
        for material in ["cells", "nuclei"]:
            m_a = metrics.query(
                f"(sample == '{e}') & (umi > 100) & (material_type == '{material}')"
            )
            m = m_a.nlargest(10000, columns="umi")
            kwds = dict(label=e + " - " + material)
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
        f"results/PD200/PD200_scifi_condition_comparison.species_mixing.{label}.separate.svg",
        **fig_kws,
    )
