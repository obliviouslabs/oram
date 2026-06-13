import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress


def parse_counts(log_file):
    counts = {}
    with open(log_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                counts[int(parts[0])] = int(parts[1])
    return counts


def tail_points(counts):
    total_samples = sum(counts.values())
    if total_samples == 0:
        return [], []

    k_vals = []
    log_probs = []
    cumulative_count = total_samples
    for k in range(max(counts.keys()) + 1):
        freq = counts.get(k, 0)
        prob_exceed = cumulative_count / total_samples
        if prob_exceed > 0 and cumulative_count < total_samples:
            k_vals.append(k)
            log_probs.append(np.log2(prob_exceed))

        cumulative_count -= freq
        if cumulative_count <= 0:
            break

    return k_vals, log_probs


def distribution_points(counts):
    total_samples = sum(counts.values())
    if total_samples == 0:
        return [], []

    k_vals = list(range(max(counts.keys()) + 1))
    log_probs = []
    for k in k_vals:
        freq = counts.get(k, 0)
        log_probs.append(np.log2(freq / total_samples) if freq > 0 else np.nan)
    return k_vals, log_probs


def label_for(log_file):
    return Path(log_file).stem.replace("_", " ")


def default_output(log_files):
    if len(log_files) == 1:
        path = Path(log_files[0])
        return str(path.with_name(path.stem + "_tail_bound.pdf"))
    return "tail_bound_comparison.pdf"


def print_regression(log_file, slope, intercept, r_value):
    print(f"Linear Regression Results for {log_file}:")
    print(f"Slope (m): {slope:.4f}")
    print(f"Intercept (c): {intercept:.4f}")
    print(f"R-squared: {r_value**2:.4f}")
    print(f"Equation: log2(P(X >= k)) = {slope:.4f} * k + {intercept:.4f}")
    print("-" * 40)

    failure_probs_log2 = [-40, -60, -64, -80]
    print("Required Stash Size (k) to achieve specific failure probabilities:")
    for target_log2 in failure_probs_log2:
        required_k = (target_log2 - intercept) / slope
        print(
            f"Probability 2^{target_log2}: {required_k:.2f} "
            f"(ceiling: {int(np.ceil(required_k))})"
        )
    print()


def analyze_tail_bound(log_files, output_pdf):
    datasets = []
    for log_file in log_files:
        counts = parse_counts(log_file)
        total_samples = sum(counts.values())
        if total_samples == 0:
            print(f"No valid data found in {log_file}.")
            continue

        k_tail, log_tail = tail_points(counts)
        if len(k_tail) < 2:
            print(f"Not enough tail data points in {log_file}.")
            continue

        slope, intercept, r_value, _, _ = linregress(k_tail, log_tail)
        print_regression(log_file, slope, intercept, r_value)

        k_dist, log_dist = distribution_points(counts)
        datasets.append(
            {
                "label": label_for(log_file),
                "k_dist": np.array(k_dist),
                "log_dist": np.array(log_dist),
                "k_tail": np.array(k_tail),
                "log_tail": np.array(log_tail),
                "slope": slope,
                "intercept": intercept,
            }
        )

    if not datasets:
        raise RuntimeError("No plottable data found.")

    fig, (dist_ax, tail_ax) = plt.subplots(2, 1, figsize=(9, 10), sharex=True)

    for dataset in datasets:
        label = dataset["label"]
        dist_ax.plot(
            dataset["k_dist"],
            dataset["log_dist"],
            marker="o",
            linestyle="-",
            alpha=0.8,
            label=label,
        )

        tail_ax.plot(
            dataset["k_tail"],
            dataset["log_tail"],
            marker="o",
            linestyle="",
            alpha=0.7,
            label=f"{label} data",
        )

        fitted_line = dataset["slope"] * dataset["k_tail"] + dataset["intercept"]
        tail_ax.plot(
            dataset["k_tail"],
            fitted_line,
            linestyle="-",
            label=f"{label} fit ({dataset['slope']:.2f}k {dataset['intercept']:+.2f})",
        )

    dist_ax.set_ylabel("$\\log_2(P(X = k))$")
    dist_ax.set_title("Stash Load Distribution")
    dist_ax.legend()
    dist_ax.grid(True)

    tail_ax.set_xlabel("Stash Size (k)")
    tail_ax.set_ylabel("$\\log_2(P(X \\geq k))$")
    tail_ax.set_title("Tail Bound Analysis")
    tail_ax.legend()
    tail_ax.grid(True)

    fig.tight_layout()
    fig.savefig(output_pdf)
    print(f"Graph saved to {output_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot stash-load distributions and tail bounds from logs."
    )
    parser.add_argument("log_files", nargs="+", help="log files to analyze")
    parser.add_argument(
        "-o",
        "--output",
        help="output PDF path; defaults to a per-log or comparison filename",
    )
    args = parser.parse_args()

    analyze_tail_bound(args.log_files, args.output or default_output(args.log_files))


if __name__ == "__main__":
    main()
