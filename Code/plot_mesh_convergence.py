#!/usr/bin/env python3
"""Plot synthetic mesh-convergence data with adjustable order.

Run from the desired output directory:
    python ../Code/plot_fake_mesh_convergence.py --order 1.9

The residual follows residual proportional to ni**(-order), with
reproducible multiplicative noise. The graph also shows 0th- through
3rd-order reference trends.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


I_SIZES = np.array([53, 106, 212, 424, 848])
REFERENCE_RESIDUAL = 1e-4
NOISE_SEED = 42
TREND_REFERENCE_X_FACTOR = 0.8
TREND_REFERENCE_Y_FACTOR = 1.5
X_AXIS_MAX_FACTOR = 1.5
TREND_ORDERS = ((0, "0th"), (1, "1st"), (2, "2nd"), (3, "3rd"))
TREND_COLORS = ("tab:red", "tab:green", "tab:purple", "tab:brown")


def make_synthetic_residuals(order, noise):
    ideal_residuals = (
        REFERENCE_RESIDUAL * (I_SIZES / I_SIZES[0]) ** -order
    )
    random_generator = np.random.default_rng(NOISE_SEED)
    residuals = ideal_residuals * np.exp(
        random_generator.normal(0, noise, size=I_SIZES.size)
    )
    if np.any(~np.isfinite(residuals)) or np.any(residuals <= 0):
        raise ValueError(
            "Synthetic residuals are outside the representable range; "
            "use a smaller order or noise level"
        )
    return residuals


def plot_synthetic_convergence(order, noise, output_file):
    residuals = make_synthetic_residuals(order, noise)
    figure, axis = plt.subplots(figsize=[8, 5])
    x_axis_min = I_SIZES[0] / 1.2
    x_axis_max = I_SIZES[-1] * X_AXIS_MAX_FACTOR
    axis.loglog(
        I_SIZES,
        residuals,
        marker="o",
        color="tab:blue",
        label=f"Bump test case",
    )

    trend_i_sizes = np.geomspace(x_axis_min, x_axis_max, num=200)
    reference_ni = I_SIZES[0] * TREND_REFERENCE_X_FACTOR
    reference_residual = np.max(residuals) * TREND_REFERENCE_Y_FACTOR
    for (trend_order, order_name), color in zip(TREND_ORDERS, TREND_COLORS):
        trend_residuals = (
            reference_residual
            * (trend_i_sizes / reference_ni) ** -trend_order
        )
        axis.loglog(
            trend_i_sizes,
            trend_residuals,
            linestyle="--",
            color=color,
            linewidth=1,
            label=f"{order_name} order trend",
        )
    axis.set_xlabel("Number of i-direction cells ($n_i$)", fontsize=14)
    axis.set_ylabel("Average density residual", fontsize=14)
    axis.set_xlim(x_axis_min, x_axis_max)
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    axis.legend(fontsize="small")
    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved synthetic convergence plot as {output_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--order",
        type=float,
        default=1.9,
        help="order of the fake residual scaling (default: 1.9)",
    )
    parser.add_argument(
        "--noise",
        type=float,
        default=0.02,
        help="multiplicative log-space noise level (default: 0.02, about 2%%)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output image path (default: fake_mesh_convergence_order_<order>.png)",
    )
    arguments = parser.parse_args()
    if not np.isfinite(arguments.order) or arguments.order < 0:
        parser.error("--order must be a finite, non-negative number")
    if not np.isfinite(arguments.noise) or arguments.noise < 0:
        parser.error("--noise must be a finite, non-negative number")

    output_file = arguments.output or Path.cwd() / (
        f"fake_mesh_convergence_order_{arguments.order:g}.png"
    )
    plot_synthetic_convergence(arguments.order, arguments.noise, output_file)


if __name__ == "__main__":
    main()
