#!/usr/bin/env python3
#
# Plot final density residuals from previously saved CFL runs.
#
# From the Cases directory, run:
#     python ../Code/plot_conv_from_CFL_csv.py
# Optionally specify one or more case names:
#     python ../Code/plot_conv_from_CFL_csv.py bend bump

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


def read_case_residuals(case_name, case_dir):
    prefix = f"conv_{case_name}_CFL_"
    residual_files = list(case_dir.glob(f"{prefix}*.csv"))
    if not residual_files:
        raise FileNotFoundError(
            f"No saved CFL residual files found for {case_name} in {case_dir}"
        )

    cfl_residual_pairs = []
    for residual_file in residual_files:
        try:
            cfl = float(residual_file.stem[len(prefix):])
        except ValueError as error:
            raise ValueError(
                f"Could not read CFL value from {residual_file.name}"
            ) from error
        if not np.isfinite(cfl) or cfl <= 0:
            raise ValueError(
                f"CFL value in {residual_file.name} must be finite and positive"
            )

        residual_data = np.loadtxt(residual_file, ndmin=2)
        if residual_data.shape[1] < 9:
            raise ValueError(f"{residual_file} must contain nine columns")
        final_density_residual = residual_data[-1, 1]
        if not np.isfinite(final_density_residual) or final_density_residual <= 0:
            raise ValueError(
                f"Final average density residual in {residual_file} "
                "must be finite and positive"
            )
        cfl_residual_pairs.append((cfl, final_density_residual))

    cfl_residual_pairs.sort(key=lambda pair: pair[0])
    cfl_values, residuals = zip(*cfl_residual_pairs)
    return np.asarray(cfl_values), np.asarray(residuals)


def plot_residuals(case_residuals, output_dir):
    figure, axis = plt.subplots(figsize=[8, 6])
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    order_colors = {0: "#000000", 1: "#006400", 2: "#4B0082"}

    for (case_name, (cfl_values, residuals)), color in zip(
        case_residuals.items(), colors
    ):
        axis.plot(
            cfl_values,
            residuals,
            marker="o",
            markersize=3,
            color=color,
            label=f"{case_name} residual",
        )

        if case_name == "bend":
            reference_value = residuals[-1]
            reference_cfl = cfl_values[-1]
            for order in (0, 1, 2):
                reference_residuals = (
                    reference_value * (cfl_values / reference_cfl) ** order
                )
                axis.plot(
                    cfl_values,
                    reference_residuals,
                    linestyle="--",
                    color=order_colors[order],
                    linewidth=0.8,
                    alpha=1,
                    label=f"{order}th order ($CFL^{order}$)",
                )

    case_names = " and ".join(case_residuals)
    axis.set_xlabel("CFL number", fontsize=14)
    axis.set_ylabel("Average density ($\\rho$) residual", fontsize=14)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.5)
    axis.legend(fontsize="small")
    figure.tight_layout()
    plot_name = f"convergence_plot_{'_'.join(case_residuals)}_CFL.png"
    figure.savefig(output_dir / plot_name, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main():
    case_names = sys.argv[1:] or ["bend", "bump"]
    case_dir = Path.cwd()
    case_residuals = {
        case_name: read_case_residuals(case_name, case_dir)
        for case_name in case_names
    }
    plot_residuals(case_residuals, case_dir)
    print(f"Saved plot as convergence_plot_{'_'.join(case_names)}_CFL.png")


if __name__ == "__main__":
    main()
