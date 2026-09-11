#!/usr/bin/env python3
"""Saddle-node analysis of the stability monitor output.

Near a saddle-node (fold) instability at load alpha_c the lowest eigenvalue of the
Hessian at the equilibrium vanishes like

    lambda_min ~ sqrt(s * (alpha_c - alpha)),   i.e.   lambda_min^2 = s * (alpha_c - alpha),

so lambda_min^2 is linear in alpha just before the instability. For every instability of a
run this script fits a straight line to lambda_min^2 on the last points before it, which
gives alpha_c and s, and checks that alpha_c falls between the last stable load step and
the step where the instability happened.

Instabilities are taken from the stress: a decrease while loading larger than --stress-tol
typical elastic increments, so events that were not saved as avalanches are included; saved
avalanches are marked.

Input (run directory): eigen_log.csv (--eig-every / --eig-at-avalanche), fractional_drops.csv
Output (run directory): saddle_node_fits.csv, saddle_node.png

Usage:
    python3 plot_saddle_node.py [run_dir] [--max-points 8] [--window 30] [--r2 0.995]
                                [--stress-tol 2]
Best with --eig-every=5 --eig-refine-ahead=2 (default), which switches to every load step
when lambda_min^2 extrapolates to zero within two grid intervals.
"""

import argparse
import csv
import math
import os
import sys

import numpy as np


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def find_events(steps, alphas, stress, saved, rel_tol):
    """Clusters of consecutive steps where the stress decreases while loading."""
    direction = math.copysign(1.0, alphas[1] - alphas[0])
    dsig = direction * np.diff(stress)  # > 0 on elastic branches
    elastic = dsig[dsig > 0]
    tol = rel_tol * (np.median(elastic) if elastic.size else 0.0)
    drop = np.concatenate([[False], dsig < -tol])
    events = []
    i = 1
    while i < len(steps):
        if drop[i]:
            j = i
            while j + 1 < len(steps) and drop[j + 1]:
                j += 1
            events.append({"first": steps[i], "last": steps[j],
                           "alpha": alphas[i],
                           "saved": bool(saved[i:j + 1].any())})
            i = j + 1
        else:
            i += 1
    return events, direction


def fit_branch(st, al, lam, event, direction, dalpha, max_points, window, r2_min):
    """Largest tail (>= 3 points) of the branch with lambda^2 linear in alpha."""
    best = None
    for k in range(3, min(max_points, len(st)) + 1):
        s_k, a_k, l_k = st[-k:], al[-k:], lam[-k:]
        if event["first"] - s_k[0] > window:
            break
        y = l_k ** 2
        B, A = np.polyfit(a_k, y, 1)
        if direction * B >= 0:  # lambda^2 must decrease towards the instability
            continue
        pred = A + B * a_k
        ss = np.sum((y - y.mean()) ** 2)
        r2 = 1.0 - np.sum((y - pred) ** 2) / ss if ss > 0 else 0.0
        if r2 < r2_min:
            continue
        alpha_c = -A / B
        best = {"n": k, "alpha_c": alpha_c, "slope": abs(B), "r2": r2,
                "first_step": s_k[0],
                "alpha_first": a_k[0], "alpha_last": a_k[-1],
                "lambda_last": l_k[-1]}
    if best is not None:
        # position of alpha_c in load steps after the last stable step (expected in (0, 1])
        last_stable = event["alpha"] - direction * dalpha
        best["offset_steps"] = direction * (best["alpha_c"] - last_stable) / dalpha
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".")
    ap.add_argument("--max-points", type=int, default=8,
                    help="max points used in a fit (default 8)")
    ap.add_argument("--window", type=int, default=30,
                    help="fit points must lie within this many steps of the event (30)")
    ap.add_argument("--r2", type=float, default=0.995,
                    help="minimum R^2 of the linear fit of lambda^2 (0.995)")
    ap.add_argument("--stress-tol", type=float, default=2.0,
                    help="stress decrease counted as an instability, in units of the median "
                         "elastic stress increment per step (2; smaller values also catch the "
                         "slight softening one step before a jump and misplace the event)")
    args = ap.parse_args()

    d = args.run_dir
    eig_path = os.path.join(d, "eigen_log.csv")
    fd_path = os.path.join(d, "fractional_drops.csv")
    for p in (eig_path, fd_path):
        if not os.path.exists(p):
            sys.exit(f"missing {p} (run with --eig-every=N)")

    eig = {}
    for r in read_csv(eig_path):
        if r["lambda_1"]:
            eig[int(r["Iteration"])] = (float(r["Alpha"]), float(r["lambda_1"]))
    e_steps = np.array(sorted(eig))
    e_alpha = np.array([eig[s][0] for s in e_steps])
    e_lam = np.array([eig[s][1] for s in e_steps])

    fd = read_csv(fd_path)
    steps = np.array([int(r["Iteration"]) for r in fd])
    alphas = np.array([float(r["Alpha"]) for r in fd])
    stress = np.array([float(r["PostStress"]) for r in fd])
    saved = np.array([r["StressDropDetected"] == "1" for r in fd])
    dalpha = abs(alphas[1] - alphas[0])

    events, direction = find_events(steps, alphas, stress, saved, args.stress_tol)

    fits = []
    prev_last = -1
    for ev in events:
        sel = (e_steps > prev_last) & (e_steps < ev["first"])
        ev["prev_last"] = prev_last
        prev_last = ev["last"]
        f = None
        if sel.sum() >= 3:
            f = fit_branch(e_steps[sel], e_alpha[sel], e_lam[sel], ev, direction, dalpha,
                           args.max_points, args.window, args.r2)
        fits.append((ev, f))

    out_csv = os.path.join(d, "saddle_node_fits.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["EventStep", "EventAlpha", "SavedAvalanche", "FitPoints", "AlphaC",
                    "Slope_s", "R2", "AlphaCOffsetSteps", "LambdaLast"])
        for ev, fit in fits:
            if fit:
                w.writerow([ev["first"], ev["alpha"], int(ev["saved"]), fit["n"],
                            fit["alpha_c"], fit["slope"], fit["r2"], fit["offset_steps"],
                            fit["lambda_last"]])
            else:
                w.writerow([ev["first"], ev["alpha"], int(ev["saved"]), 0, "", "", "", "", ""])

    good = [(ev, fit) for ev, fit in fits if fit]
    consistent = [(ev, fit) for ev, fit in good if -0.5 <= fit["offset_steps"] <= 1.5]
    print(f"{len(events)} instabilities ({sum(e['saved'] for e in events)} saved avalanches); "
          f"lambda^2 linear fit found for {len(good)}, alpha_c within the last step "
          f"(+-0.5) for {len(consistent)}")
    print(f"fits written to {out_csv}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 8.5))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 2, 3)
    ax3 = fig.add_subplot(2, 2, 4)

    # (a) lambda_min along the loading with the fitted square roots
    ax1.semilogy(e_alpha, e_lam, ".", ms=2, color="0.6", label=r"$\lambda_{\min}$")
    for ev in events:
        ax1.axvline(ev["alpha"], color="tab:red" if ev["saved"] else "0.8",
                    lw=0.6 if ev["saved"] else 0.3, zorder=0)
    for ev, fit in good:
        a = np.linspace(fit["alpha_first"], fit["alpha_c"], 60)
        lam = np.sqrt(np.clip(fit["slope"] * direction * (fit["alpha_c"] - a), 0, None))
        ok = -0.5 <= fit["offset_steps"] <= 1.5
        ax1.plot(a, lam, "-", lw=0.8, color="tab:blue" if ok else "tab:orange")
    ax1.set_xlabel(r"$\alpha$")
    ax1.set_ylabel(r"$\lambda_{\min}$")
    ax1.set_title(r"$\lambda_{\min}$ at relaxed states; fits $\lambda_{\min}^2=s(\alpha_c-\alpha)$ "
                  "(blue: $\\alpha_c$ within the last load step); red: saved avalanches, "
                  "grey: other instabilities", fontsize=9)

    # (b) collapse: lambda / sqrt(s) against distance to alpha_c, slope 1/2; all points of
    #     the branch within --window steps before the instability (fitted or not)
    xs, ys = [], []
    for ev, fit in consistent:
        sel = (e_steps >= ev["first"] - args.window) & (e_steps < ev["first"]) & \
              (e_steps > ev.get("prev_last", -1))
        x = direction * (fit["alpha_c"] - e_alpha[sel]) / dalpha
        y = e_lam[sel] / math.sqrt(fit["slope"] * dalpha)
        keep = x > 0
        xs.append(x[keep]); ys.append(y[keep])
        ax2.loglog(x[keep], y[keep], ".-", lw=0.5, ms=3, alpha=0.7)
    if xs:
        xx = np.logspace(math.log10(max(min(np.concatenate(xs)), 1e-3)),
                         math.log10(max(np.concatenate(xs))), 50)
        ax2.loglog(xx, np.sqrt(xx), "k--", lw=1.5, label="slope 1/2")
        ax2.legend(fontsize=8)
    ax2.set_xlabel(r"$(\alpha_c-\alpha)/\Delta\alpha$  (load steps to the instability)")
    ax2.set_ylabel(r"$\lambda_{\min}/\sqrt{s\,\Delta\alpha}$")
    ax2.set_title(f"collapse of {len(consistent)} instabilities", fontsize=9)

    # (c) where alpha_c falls relative to the last stable step
    offs = [fit["offset_steps"] for ev, fit in good]
    if offs:
        ax3.hist(np.clip(offs, -3, 4), bins=np.arange(-3, 4.25, 0.25), color="tab:blue")
        ax3.axvspan(0, 1, color="tab:green", alpha=0.15, label="between last stable\n"
                                                                 "and unstable step")
        ax3.legend(fontsize=8)
    ax3.set_xlabel(r"$(\alpha_c-\alpha_{\mathrm{last\ stable}})/\Delta\alpha$")
    ax3.set_ylabel("instabilities")
    ax3.set_title(f"predicted $\\alpha_c$ ({len(good)} fits of {len(events)} instabilities)",
                  fontsize=9)

    fig.tight_layout()
    out_png = os.path.join(d, "saddle_node.png")
    fig.savefig(out_png, dpi=130)
    print(f"plot written to {out_png}")


if __name__ == "__main__":
    main()
