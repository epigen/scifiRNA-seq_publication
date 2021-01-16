a = sc.read(
    "results/PD195-1_Tcells_765k-sample1_P7/PD195-1_Tcells_765k-sample1_P7.single_cell.expression.min_200_umis.384_wells.filtered.processed.h5ad"
)

a = a.raw.to_adata()

a.var = a.var.drop(["external_gene_name"], axis=1)
a.var = a.var.join(biomart)

genes = ["CD3G", "CD247", "CD4", "CD8A", "FOXP3"]
genes = [g for g in genes if g in a.var["external_gene_name"].tolist()]

fig, axes = plt.subplots(
    2, 4, figsize=(4 * 4, 2 * 4), gridspec_kw=dict(hspace=0.5, wspace=0.5)
)
axes = axes.flatten()
for ax, f in zip(axes[:3], ["leiden", "donor_id", "activation"]):
    sc.pl.umap(a, color=f, show=False, ax=ax)

for ax, gene in zip(axes[3:], genes):
    x = (
        a[:, a.var.query(f"external_gene_name == '{gene}'").index]
        .X.todense()
        .A1
    )
    vmin, vmax = np.percentile(x, 5), np.percentile(x, 95)
    sc.pl.umap(
        a,
        color=gene,
        gene_symbols="external_gene_name",
        vmin=-vmax,
        vmax=vmax,
        color_map="coolwarm",  #  s=50, alpha=0.75,
        show=False,
        ax=ax,
    )
    ax.set(title=gene)

for ax in axes:
    ax.get_children()[0].set_rasterized(True)
fig.savefig(
    "results/Reviewer_figure.T-cell_CD4_CD8_expression.umap.svg",
    dpi=300,
    bbox_inches="tight",
)

s = pd.Series(a.X.mean(0).A1, index=a.var.index).sort_values(ascending=False)
s = np.log1p(s)
r = s.rank(ascending=False)

fig, axes = plt.subplots(1, 2, figsize=(2 * 4, 3))
for ax in axes:
    ax.plot(r, s)
    ax.set(xlabel="Genes (rank mean expression)", ylabel="Mean log expression")
    for g in genes:
        e = a.var.query(f"external_gene_name == '{g}'").index[0]
        ax.axvline(r.loc[e], s.loc[e], linestyle="--", color="grey", alpha=0.5)
        ax.text(r.loc[e], s.loc[e], s=g)
axes[1].set_xscale("log")
fig.savefig(
    "results/Reviewer_figure.T-cell_CD4_CD8_expression.rank.svg",
    dpi=300,
    bbox_inches="tight",
)
