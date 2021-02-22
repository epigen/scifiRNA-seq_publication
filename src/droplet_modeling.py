#! /usr/bin/env python

"""
Modeling droplet fill rates and collision rates of the 10X Chromium device.
"""

import sys
from typing import List, Union, Tuple, Dict, Callable, Optional, Collection
from pathlib import Path
import itertools

try:
    from functools import cache  # Python 3.9
except ImportError:
    from functools import lru_cache as cache  # Python <3.9
import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import scipy
from sklearn.linear_model import LinearRegression

import pymc3 as pm
from pymc3.variational.callbacks import CheckParametersConvergence

# Custom types
Array = Union[np.ndarray]
Fig = Union[matplotlib.figure.Figure]
DataFrame = Union[pd.DataFrame]

distributions = [
    "Poisson",
    "ZeroInflatedPoisson",
    "NegativeBinomial",
    "ZeroInflatedNegativeBinomial",
]


def main():

    for version in droplet_files:
        print(version)

        # Load observed counts of nuclei per droplet
        droplet_counts = get_droplet_counts(droplet_files[version])
        # counts[version] = droplet_counts

        # # some preparations for plotting
        exprs = droplet_counts["loaded_nuclei"].unique()

        droplet_counts_p = droplet_counts.pivot_table(
            index="cells_per_droplet",
            columns="loaded_nuclei",
            values="count",
            fill_value=0,
        )
        print(droplet_counts_p)

        # Expand distributions into the actual observed data
        counts = counts_to_observations(droplet_counts_p, tailor_and_trim=True)
        counts = np.asarray(counts).T
        print(counts.shape)

        for distribution in distributions:
            print(version, distribution)
            output_prefix = (
                output_dir / f"droplets_{version}.{distribution}."
            ).as_posix()
            # if os.path.exists(output_prefix + "estimated_parameters.svg"):
            #     continue

            # Get model for distribution and shape compatible with dataset
            model = get_model(distribution, counts)

            # Sample with MCMC
            mcmc_trace, ppc_trace = sample_from_model(model)

            summary = pm.summary(mcmc_trace).round(2)
            print(np.allclose(summary["r_hat"], 1))
            diverg = mcmc_trace["diverging"].sum()
            print(f"Number of divergent: {diverg}")
            divperc = diverg / len(mcmc_trace) * 100
            print(f"Percent divergent: {divperc}")

            # Just for fun, let's also do Variational Inference too
            advi, vi_trace, mean_field, tracker = inference_with_model(model)

            # Plot
            # # MCMC
            fig = plot_mcmc(mcmc_trace)
            fig.savefig(output_prefix + "mcmc_trace.svg", **figkws)

            fig = plot_divergences(mcmc_trace)
            fig.savefig(output_prefix + "mcmc_divergences.svg", **figkws)

            # # VI
            fig = plot_vi(advi, mean_field, tracker)
            fig.savefig(output_prefix + "vi_trace.svg", **figkws)

            # Get parameters
            params = get_parameters(model, mcmc_trace, vi_trace, labels=exprs)
            params.to_csv(output_prefix + "estimated_parameters.csv")
            params = pd.read_csv(
                output_prefix + "estimated_parameters.csv", index_col=0
            )

            fig = plot_parameter_estimates(params, distribution)
            fig.savefig(output_prefix + "estimated_parameters.svg", **figkws)

            funcs, fig = interpolate_over_parameters(params, kind="linear")
            fig.savefig(output_prefix + "parameter_interpolation.svg", **figkws)

            collision_estimates = predict_collisions(funcs, dist=distribution)
            collision_estimates.to_csv(
                output_prefix + "collision_estimates.csv"
            )

            fig = plot_collision_estimates(
                collision_estimates, dist=distribution
            )
            fig.savefig(output_prefix + "predicted_collisions.svg", **figkws)

    # Plot mean estimate across datasets and models
    fig = joint_mean_plots()
    fig.savefig(output_dir / "predicted_means.all_models.svg", **figkws)

    # Plot collision estimates
    fig, axes = plt.subplots(2, 4, figsize=(4 * 5, 2 * 5))
    for i, version in enumerate(droplet_files):
        for j, distribution in enumerate(distributions):
            print(version, distribution)
            output_prefix = (
                output_dir / f"droplets_{version}.{distribution}."
            ).as_posix()
            colisions = pd.read_csv(
                output_prefix + "collision_estimates.csv", index_col=[0, 1]
            )
            plot_collision_estimates(
                colisions, dist=distribution, ax=axes[i][j]
            )
    fig.savefig(output_dir / "predicted_collisions.all_models.svg", **figkws)


def get_droplet_counts(droplet_file):
    droplet_counts = pd.read_csv(droplet_file)

    # # split apply combine for fraction normalization
    droplet_counts = droplet_counts.join(
        droplet_counts.groupby("loaded_nuclei")
        .apply(
            lambda x: pd.Series(
                (x["count"] / x["count"].sum()).tolist(),
                index=x.index,
                name="count (%)",
            )
        )
        .reset_index(drop=True)
    )
    # droplet_counts['norm_count'] += sys.float_info.epsilon

    # # split apply combine for % max normalization
    droplet_counts = droplet_counts.join(
        droplet_counts.groupby("loaded_nuclei")
        .apply(
            lambda x: pd.Series(
                ((x["count"] / x["count"].max()) * 100).tolist(),
                index=x.index,
                name="count (% max)",
            )
        )
        .reset_index(drop=True)
    )

    # # split apply combine for number of droplets counted
    droplet_counts = (
        droplet_counts.reset_index()
        .set_index("loaded_nuclei")
        .join(
            droplet_counts.groupby("loaded_nuclei")["count"]
            .sum()
            .rename("droplets_analyzed")
        )
        .reset_index()
        .drop("index", axis=1)
    )
    return droplet_counts


def counts_to_observations(
    droplet_counts_p, tailor_and_trim: bool = False
) -> List[Array]:
    counts = [
        np.array(
            [
                i
                for i in droplet_counts_p.index
                for _ in range(droplet_counts_p.loc[i, j])
            ]
        )
        for j in droplet_counts_p.columns
    ]

    # # # # cap all to min(observations) ~= 609 for first dataset
    if tailor_and_trim:
        n_exp = len(counts)
        m = min(map(len, counts))
        counts = [np.random.choice(x, m) for x in counts]
    return counts


# @cache
def get_model(dist, data) -> pm.Model:
    means = data.mean(0)
    n_exp = data.shape[1]
    if dist == "Poisson":
        with pm.Model() as poi_model:
            lam = pm.Exponential("lam", lam=means, shape=(1, n_exp))
            poi = pm.Poisson(
                "poi",
                mu=lam,
                observed=data,
            )
        return poi_model
    if dist == "ZeroInflatedPoisson":
        with pm.Model() as zip_model:
            psi = pm.Uniform("psi", shape=(1, n_exp))
            lam = pm.Exponential("lam", lam=means, shape=(1, n_exp))
            zip = pm.ZeroInflatedPoisson(
                "zip",
                psi=psi,
                theta=lam,
                observed=data,
            )
        return zip_model
    if dist == "NegativeBinomial":
        with pm.Model() as nb_model:
            gamma = pm.Gamma("gm", 0.01, 0.01, shape=(1, n_exp))
            lam = pm.Exponential("lam", lam=means, shape=(1, n_exp))
            nb = pm.NegativeBinomial(
                "nb",
                alpha=gamma,
                mu=lam,
                observed=data,
            )
        return nb_model
    if dist == "ZeroInflatedNegativeBinomial":
        with pm.Model() as zinb_model:
            gamma = pm.Gamma("gm", 0.01, 0.01, shape=(1, n_exp))
            lam = pm.Exponential("lam", lam=means, shape=(1, n_exp))
            psi = pm.Uniform("psi", shape=(1, n_exp))
            zinb = pm.ZeroInflatedNegativeBinomial(
                "zinb",
                psi=psi,
                alpha=gamma,
                mu=lam,
                observed=data,
            )
        return zinb_model


def sample_from_model(model):
    with model:
        trace = pm.sample(**sampler_params)[TAKE_AFTER:]

    # We can now sample from the model posterior
    with model:
        ppc_trace = pm.sample_posterior_predictive(
            trace=trace,
            samples=sampler_params["draws"] * 2,
            random_seed=sampler_params["random_seed"],
        )
    return trace, ppc_trace


def inference_with_model(model):
    with model:
        advi = pm.ADVI()
        tracker = pm.callbacks.Tracker(
            mean=advi.approx.mean.eval, std=advi.approx.std.eval
        )
        mean_field = advi.fit(
            n=vi_params["n"],
            callbacks=[CheckParametersConvergence(), tracker],
        )
    vi_trace = mean_field.sample(draws=sampler_params["draws"])
    return advi, vi_trace, mean_field, tracker


def plot_mcmc(mcmc_trace) -> Fig:

    plt.close("all")
    axes = pm.traceplot(mcmc_trace)
    return axes[0, 0].figure


def plot_divergences(trace) -> Fig:
    params = trace.varnames
    combs = list(itertools.combinations(params, 2))

    n, m = get_grid_dims(len(combs))
    fig, axes = plt.subplots(n, m, figsize=(m * 3, n * 3), squeeze=False)
    for i, (a, b) in enumerate(combs):
        ax = axes.flatten()[i]
        x = trace.get_values(varname=a, combine=True)[:, 0].mean(1)
        y = trace.get_values(varname=b, combine=True)[:, 0].mean(1)
        ax.scatter(x, y, color="grey", alpha=0.1, rasterized=True)
        divergent = trace["diverging"]
        ax.scatter(x[divergent], y[divergent], color="red", alpha=0.5)
        ax.set(xlabel=a, ylabel=b)
        p = divergent.sum() / divergent.shape[0] * 100
    fig.suptitle(f"{divergent.sum()} divergences ({p:.3f}%)")
    return fig


def get_grid_dims(
    dims: Union[int, Collection], _nstart: Optional[int] = None
) -> Tuple[int, int]:
    """
    Given a number of `dims` subplots, choose optimal x/y dimentions of plotting
    grid maximizing in order to be as square as posible and if not with more
    columns than rows.
    """
    if not isinstance(dims, int):
        dims = len(dims)
    if _nstart is None:
        n = min(dims, 1 + int(np.ceil(np.sqrt(dims))))
    else:
        n = _nstart
    if (n * n) == dims:
        m = n
    else:
        a = pd.Series(n * np.arange(1, n + 1)) / dims
        m = a[a >= 1].index[0] + 1
    assert n * m >= dims

    if n * m % dims > 1:
        try:
            n, m = get_grid_dims(dims=dims, _nstart=n - 1)
        except IndexError:
            pass
    return n, m


def plot_vi(advi, mean_field, tracker) -> Fig:
    fig = plt.figure(figsize=(8, 6))
    mu_ax = fig.add_subplot(221)
    std_ax = fig.add_subplot(222)
    hist_ax = fig.add_subplot(212)
    mu_ax.plot(tracker["mean"])
    mu_ax.set_title("Mean")
    std_ax.plot(tracker["std"])
    std_ax.set_title("SD")
    hist_ax.plot(advi.hist)
    hist_ax.set_title("Negative ELBO")
    return fig


def get_parameters(model, mcmc_trace, vi_trace, labels=None) -> DataFrame:
    # # Let's gather the model parameters in a dataframe
    _res = dict()
    vars_ = [t.name for t in model.unobserved_RVs[1:]]
    for param in vars_:
        for red_func, metric in [(np.mean, "_mean"), (np.std, "_std")]:
            for sampler, s_label in [(mcmc_trace, "_mcmc"), (vi_trace, "_vi")]:
                _res[param + metric + s_label] = red_func(
                    sampler[param], axis=0
                ).flatten()
    res = pd.DataFrame(
        _res, index=pd.Series(labels) if labels is not None else None
    ).rename_axis("loaded_nuclei", axis=0)
    return res


def plot_parameter_estimates(res: DataFrame, dist: str) -> Fig:
    # # # Plot estimated lamba
    params = res.columns.str.extract(r"(.*?)_")[0].unique().tolist()
    rows = len(params)
    fig, axis = plt.subplots(
        rows,
        2,
        figsize=(2 * 4, rows * 4),
        sharey="row",
        sharex="row",
        squeeze=False,
    )
    fig.suptitle(
        f"Nuclei loading counts modeled as output of a {dist} function",
        fontsize=16,
    )
    for i, param in enumerate(params):
        for j, method in enumerate(["mcmc", "vi"]):
            param_name = param.replace("lam", "lambda").replace("gm", "gamma")
            axis[i, j].set_title("\n" + method.upper(), fontsize=14)
            axis[i, j].plot(
                res.index,
                res[f"{param}_mean_{method}"],
                color=colors[0],
                marker="o",
            )
            axis[i, j].fill_between(
                res.index,
                res[f"{param}_mean_{method}"]
                - res[f"{param}_std_{method}"] * 3,
                res[f"{param}_mean_{method}"]
                + res[f"{param}_std_{method}"] * 3,
                color=colors[i],
                alpha=0.2,
            )
            axis[i, j].tick_params(axis="y", labelcolor=colors[i])
            axis[i, j].set_xlabel("Nuclei loaded")
            axis[i, j].set_ylabel(f"$\\{param_name}$", color=colors[i])
    fig.tight_layout()
    return fig


def interpolate_over_parameters(
    params: DataFrame, kind: str = "quadratic"
) -> Tuple[Dict[str, Callable], Fig]:
    # # Okay, so we interpolate over loading concentrations and get the values of each parameter
    ps = params.columns.str.extract(r"(.*?)_")[0].unique().tolist()

    x = np.linspace(1000, params.index.max())
    funcs = dict()
    # # # Plot estimated lamba
    rows = len(ps)
    fig, axes = plt.subplots(rows, 1, figsize=(5, rows * 5), squeeze=False)
    for i, param in enumerate(ps):
        param_name = param.replace("lam", "lambda").replace("gm", "gamma")
        f = scipy.interpolate.interp1d(
            params.index,
            params[f"{param}_mean_mcmc"],
            kind=kind,
            fill_value="extrapolate",
        )
        y = f(x)
        f.__name__ = param
        funcs[param] = f
        ax = axes.flatten()[i]
        ax.plot(
            params.index,
            params[f"{param}_mean_mcmc"],
            color=colors[i],
            marker="o",
            linestyle="--",
        )
        ax.plot(x, y, color=colors[0], linestyle="-")
        ax.tick_params(axis="y", labelcolor=colors[i])
        ax.set_xlabel("Nuclei loaded")
        ax.set_ylabel(f"$\\{param_name}$", color=colors[i])
    return funcs, fig


def predict_collisions(funcs: Dict[str, Callable], dist: str) -> DataFrame:
    # # Now we sample across the X
    # x = np.linspace(1000, res.index.max(), 10)
    x = np.logspace(3, 6.5, num=10, base=10)
    n = int(1e5)
    _collision_estimates = dict()
    for barcodes in BARCODE_COMBOS:
        for input_nuclei in x:
            dist_inputs = dict()
            for param in funcs:
                val = float(funcs[param](input_nuclei / barcodes))
                if param == "psi":
                    # bounding PSI to [0, 1] since the extrapolation is unbounded
                    val = min(val, 1.0)
                    val = max(val, 0.0)
                if param in ["theta", "gm"]:
                    val = max(val, 0.01)

                # adapt parameter names to PyMC3 convention
                if dist == "ZeroInflatedPoisson":
                    key = param.replace("lam", "theta")
                else:
                    key = param.replace("lam", "mu").replace("gm", "alpha")

                dist_inputs[key] = val
            distf = getattr(pm.distributions, dist)

            # sample from distribution with interpolated parameters
            s = distf.dist(**dist_inputs).random(size=n)
            # Now we simply count fraction of collisions (droplets with more than one nuclei)
            _collision_estimates[(barcodes, input_nuclei)] = list(
                dist_inputs.values()
            ) + [(s > 1).sum() / n]

    collision_estimates = pd.DataFrame(
        _collision_estimates, index=list(funcs.keys()) + ["collision_rate"]
    ).T
    collision_estimates.index.names = ["barcodes", "loaded_nuclei"]
    return collision_estimates


def plot_collision_estimates(
    collision_estimates: DataFrame, dist: str, ax=None
) -> Optional[Fig]:
    if ax is None:
        fig, axis = plt.subplots(1, 1, figsize=(1 * 5, 5))
    else:
        axis = ax
    ax.set_title(f"Prediction with {dist}", fontsize=16)
    for barcodes in BARCODE_COMBOS:
        axis.plot(
            collision_estimates.loc[barcodes].index,
            collision_estimates.loc[barcodes, "collision_rate"] * 100,
            linestyle="-",
            marker="o",
            label=f"{barcodes} round1 barcodes",
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Nuclei loaded")
    axis.set_ylabel("% collisions")
    axis.legend()
    for x_ in [100, 10, 5, 1]:
        axis.axhline(x_, color="grey", linestyle="--")
    for y_ in [10500, 191000, 383000, 765000, 1530000]:
        axis.axvline(y_, color="grey", linestyle="--")
    return fig if ax is None else None


def joint_mean_plots() -> Fig:
    _means = list()
    _stds = list()
    for version in droplet_files:
        for distribution in distributions:
            output_prefix = (
                output_dir / f"droplets_{version}.{distribution}."
            ).as_posix()
            est = pd.read_csv(
                output_prefix + "estimated_parameters.csv", index_col=0
            )
            est.index = (
                est.index.to_series()
                .replace(191250, 191000)
                .replace(382500, 383000)
            )
            m = est[["lam_mean_mcmc", "lam_mean_vi"]].mean(1).to_frame("lambda")
            _means.append(m.assign(version=version, dist=distribution))
            s = est[["lam_std_mcmc", "lam_std_vi"]].mean(1).to_frame("lambda")
            _stds.append(s.assign(version=version, dist=distribution))

    means = pd.concat(_means).sort_index(1)
    means = means.pivot_table(
        index="loaded_nuclei", columns=["version", "dist"], values="lambda"
    )
    stds = pd.concat(_stds).sort_index(1)
    stds = stds.pivot_table(
        index="loaded_nuclei", columns=["version", "dist"], values="lambda"
    )

    fig, axes = plt.subplots(1, 4, figsize=(4 * 5.3, 1 * 5))
    for ax in axes[0:3]:
        for i, dist in enumerate(means.columns.levels[1]):
            for version in means.columns.levels[0]:
                ax.plot(
                    means.index,
                    means[version][dist],
                    "-.",
                    color=colors[i],
                    marker="v" if version == "v1" else "p",
                    label=version,
                )
                s = stds[version][dist] * 3
                ax.fill_between(
                    means.index,
                    means[version][dist] - s,
                    means[version][dist] + s,
                    label=dist,
                    alpha=0.2,
                    color=colors[i],
                )
        ax.set(xlabel="Nuclei loaded", ylabel="Mean nuclei per droplet")
    axes[1].set_xscale("log")
    axes[1].legend()
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")

    sns.heatmap(
        means.T, cbar_kws=dict(label="Mean nuclei per droplet"), ax=axes[3]
    )
    axes[3].set(xlabel="Nuclei loaded", ylabel="Model")
    return fig


if __name__ == "__main__":

    sns.set(
        context="paper", style="ticks", palette="colorblind", color_codes=True
    )
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["text.usetex"] = False
    figkws = dict(bbox_inches="tight", dpi=300)

    metadata_dir = Path("metadata").absolute()
    metadata_dir.mkdir(exist_ok=True)
    output_dir = Path("results").absolute() / "droplet_modelling"
    output_dir.mkdir(exist_ok=True, parents=True)

    droplet_files = {
        "v1": metadata_dir / "droplet_counts.v1.csv",
        "NextGEM": metadata_dir / "droplet_counts.NextGEM.csv",
    }

    colors = sns.color_palette("colorblind")

    # MCMC parameters
    sampler_params = dict(
        random_seed=0,  # random seed for reproducibility
        n_init=200_000,  # number of iterations of initializer (this is actually the default)
        tune=10_000,  # number of tuning iterations (this is probably the most critical)
        draws=5_000,  # number of iterations to sample (these will be used for our parameter estimates)
        # target_accept=0.99,
    )
    TAKE_AFTER = 1_000  # number of initial iterations to discard (just as precaution we'll exclude these)

    # VI parameters
    vi_params = dict(n=100_000)  # iterations

    BARCODE_COMBOS = [1, 96, 96 * 4, 96 * 16]

    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(1)
