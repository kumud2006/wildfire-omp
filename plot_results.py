#!/usr/bin/env python3
"""
plot_results.py -- turns results/*.csv into figures (figs/*.png) and summary
tables (results/*_summary.csv + printed markdown tables for the report).

Usage:   python3 plot_results.py            (needs pandas, matplotlib)

Every figure uses the MEDIAN of the repeated runs.
  speedup    = T_serial / T_parallel          (same N, steps, terrain)
  efficiency = speedup / threads
  weak eff.  = T(1 thread, base size) / T(p threads, scaled size)
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figs")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.spines.top": False, "axes.spines.right": False})


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        print(f"[skip] {name} not found")
        return None
    df = pd.read_csv(p)
    print(f"[ok]   {name}: {len(df)} rows")
    return df


def md(df, floatfmt="{:.4f}"):
    """tiny markdown-table printer (no 'tabulate' dependency)"""
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = [floatfmt.format(v) if isinstance(v, (float, np.floating)) else str(v) for v in r]
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"       -> figs/{name}")


# ---------------------------------------------------------------- strong ---
def strong():
    df = load("strong_combined.csv")
    if df is None:
        return
    keys = ["kernel", "terrain", "n", "threads"]
    med = df.groupby(keys, as_index=False).agg(time_s=("time_s", "median"),
                                               mcups=("mcups", "median"))
    serial = med[med.kernel == "serial"][["terrain", "n", "time_s"]].rename(columns={"time_s": "t_serial"})
    par = med[med.kernel != "serial"].merge(serial, on=["terrain", "n"], how="left")
    par["speedup_vs_serial"] = par.t_serial / par.time_s
    par["efficiency_vs_serial"] = par.speedup_vs_serial / par.threads

    # Parallel scaling: compare each kernel with its own one-thread run.
    own_base = par[par.threads == 1][["kernel", "terrain", "n", "time_s"]].rename(
        columns={"time_s": "t_kernel_1"})
    par = par.merge(own_base, on=["kernel", "terrain", "n"], how="left")
    par["speedup"] = par.t_kernel_1 / par.time_s
    par["efficiency"] = par.speedup / par.threads
    par = par.sort_values(["kernel", "n", "threads"])
    par.to_csv(os.path.join(RES, "strong_summary.csv"), index=False)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    styles = {"base": "--o", "opt": "-s"}
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, max(1, par.n.nunique())))
    for ci, (n, g_n) in enumerate(par.groupby("n")):
        for k, g in g_n.groupby("kernel"):
            lab = f"{k}, N={n}"
            ax[0].plot(g.threads, g.speedup, styles.get(k, "-o"), color=colors[ci], label=lab)
            ax[1].plot(g.threads, g.efficiency, styles.get(k, "-o"), color=colors[ci], label=lab)
            ax[2].plot(g.threads, g.mcups, styles.get(k, "-o"), color=colors[ci], label=lab)
    tmax = par.threads.max()
    ax[0].plot([1, tmax], [1, tmax], ":", color="gray", label="ideal")
    ax[1].axhline(1.0, ls=":", color="gray")
    ax[0].set(title="Strong scaling: speedup", xlabel="threads", ylabel="speedup vs same kernel, 1 thread")
    ax[1].set(title="Parallel efficiency", xlabel="threads", ylabel="parallel efficiency (vs 1 thread)")
    ax[2].set(title="Throughput", xlabel="threads", ylabel="MCUPS")
    ax[0].legend(fontsize=7)
    save(fig, "strong_scaling.png")

    print("\n### Strong scaling (median of repeats)\n")
    show = par[["kernel", "n", "threads", "time_s", "speedup", "efficiency", "mcups"]]
    print(md(show))


# ------------------------------------------------------------------ weak ---
def weak():
    df = load("weak.csv")
    if df is None:
        return
    med = df.groupby(["kernel", "threads", "n"], as_index=False).agg(time_s=("time_s", "median"),
                                                                      mcups=("mcups", "median"))
    rows = []
    for k, g in med.groupby("kernel"):
        g = g.sort_values("threads")
        t1 = g.iloc[0].time_s
        g = g.assign(efficiency=t1 / g.time_s)
        rows.append(g)
    res = pd.concat(rows)
    res.to_csv(os.path.join(RES, "weak_summary.csv"), index=False)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for k, g in res.groupby("kernel"):
        ax[0].plot(g.threads, g.time_s, "-o", label=k)
        ax[1].plot(g.threads, g.efficiency, "-o", label=k)
        for _, r in g.iterrows():
            ax[0].annotate(f"N={int(r.n)}", (r.threads, r.time_s), fontsize=7,
                           textcoords="offset points", xytext=(0, 6), ha="center")
    ax[1].axhline(1.0, ls=":", color="gray")
    ax[0].set(title="Weak scaling: run time (ideal = flat)", xlabel="threads", ylabel="time (s)")
    ax[1].set(title="Weak-scaling efficiency", xlabel="threads", ylabel="T(1) / T(p)", ylim=(0, 1.4))
    ax[0].legend()
    save(fig, "weak_scaling.png")
    print("\n### Weak scaling\n")
    print(md(res[["kernel", "threads", "n", "time_s", "efficiency", "mcups"]]))


# -------------------------------------------------------------- schedule ---
def schedule():
    df = load("schedule.csv")
    if df is None:
        return
    df["config"] = df.schedule + "," + df.chunk.astype(str)
    med = df.groupby(["terrain", "kernel", "config"], as_index=False, sort=False).agg(
        time_s=("time_s", "median"), imbalance=("imbalance", "median"), mcups=("mcups", "median"))
    med.to_csv(os.path.join(RES, "schedule_summary.csv"), index=False)

    terrains = list(med.terrain.unique())
    kernels = list(med.kernel.unique())
    for metric, ylabel, fname in [("time_s", "time (s)", "schedule_time.png"),
                                  ("imbalance", "imbalance (max/mean busy time)", "schedule_imbalance.png")]:
        fig, axes = plt.subplots(len(kernels), len(terrains),
                                 figsize=(6.2 * len(terrains), 3.6 * len(kernels)), squeeze=False)
        for r, k in enumerate(kernels):
            for c, t in enumerate(terrains):
                g = med[(med.kernel == k) & (med.terrain == t)]
                ax = axes[r][c]
                colors = ["#4c72b0" if s.startswith("static") else "#dd8452" if s.startswith("dynamic")
                          else "#55a868" for s in g.config]
                ax.bar(range(len(g)), g[metric], color=colors)
                ax.set_xticks(range(len(g)))
                ax.set_xticklabels(g.config, rotation=60, ha="right", fontsize=7)
                ax.set(title=f"{k} kernel, {t} terrain", ylabel=ylabel)
                if metric == "imbalance":
                    ax.axhline(1.0, ls=":", color="gray")
                    ax.set_ylim(bottom=0.9)
        save(fig, fname)

    print("\n### Scheduling (best 3 per terrain/kernel by time)\n")
    best = med.sort_values("time_s").groupby(["terrain", "kernel"]).head(3).sort_values(["terrain", "kernel", "time_s"])
    print(md(best[["terrain", "kernel", "config", "time_s", "imbalance", "mcups"]]))


# ------------------------------------------------------------------ simd ---
def simd():
    df = load("simd.csv")
    if df is None:
        return
    med = df.groupby(["n", "threads", "tag"], as_index=False).agg(time_s=("time_s", "median"),
                                                                    mcups=("mcups", "median"))
    order = ["base", "opt-novec", "opt-simd"]
    base = med[med.tag == "base"][["n", "threads", "time_s"]].rename(columns={"time_s": "t_base"})
    med = med.merge(base, on=["n", "threads"], how="left")
    med["speedup_vs_base"] = med.t_base / med.time_s

    # Direct comparison isolates SIMD's effect within the optimized kernel.
    novec = med[med.tag == "opt-novec"][["n", "threads", "time_s"]].rename(
        columns={"time_s": "t_opt_novec"})
    med = med.merge(novec, on=["n", "threads"], how="left")
    med["speedup_simd_vs_novec"] = med.t_opt_novec / med.time_s
    med.to_csv(os.path.join(RES, "simd_summary.csv"), index=False)

    groups = med[["n", "threads"]].drop_duplicates().values.tolist()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    w = 0.25
    for ti, tag in enumerate(order):
        g = med[med.tag == tag]
        vals_t = [g[(g.n == n) & (g.threads == t)].time_s.mean() for n, t in groups]
        vals_m = [g[(g.n == n) & (g.threads == t)].mcups.mean() for n, t in groups]
        ax[0].bar(np.arange(len(groups)) + ti * w, vals_t, w, label=tag)
        ax[1].bar(np.arange(len(groups)) + ti * w, vals_m, w, label=tag)
    for a, yl in zip(ax, ["time (s)", "MCUPS"]):
        a.set_xticks(np.arange(len(groups)) + w)
        a.set_xticklabels([f"N={n}\n{t} thr" for n, t in groups], fontsize=8)
        a.set_ylabel(yl)
    ax[0].set_title("Kernel comparison: time")
    ax[1].set_title("Kernel comparison: throughput")
    ax[0].legend()
    save(fig, "simd_comparison.png")
    print("\n### SIMD / memory optimisation\n")
    print(md(med[["n", "threads", "tag", "time_s", "mcups", "speedup_vs_base", "speedup_simd_vs_novec"]]))


# ------------------------------------------------------------------ bind ---
def bind():
    df = load("bind.csv")
    if df is None:
        return
    med = df.groupby("tag", as_index=False).agg(time_s=("time_s", "median"), mcups=("mcups", "median"),
                                                imbalance=("imbalance", "median"))
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.bar(med.tag, med.time_s, color="#8172b2")
    ax.set(title="Thread binding (opt kernel)", ylabel="time (s)")
    save(fig, "binding.png")
    print("\n### Thread binding\n")
    print(md(med))


if __name__ == "__main__":
    for fn in (strong, weak, schedule, simd, bind):
        fn()
    print("\nDone. Figures are in figs/, summary CSVs in results/.")
