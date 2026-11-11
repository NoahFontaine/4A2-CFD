#!/usr/bin/env python3
"""Plot final average density residual against bump mesh i-size.

Run from the Cases directory after running run_bump_mesh_sweep.py:
    python ../Code/plot_bump_mesh_convergence.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


RESIDUAL_COLUMN_COUNT = 9
RESIDUAL_FILE_PREFIX = "conv_bump_mesh_"
OUTPUT_FILE = "convergence_plot_bump_mesh_density.png"
TREND_ORDERS = ((0, "0th"), (1, "1st"), (2, "2nd"), (3, "3rd"))
TREND_COLORS = ("tab:red", "tab:green", "tab:purple", "tab:brown")


def read_mesh_residuals(case_dir):
    residual_files = sorted(case_dir.glob(f"{RESIDUAL_FILE_PREFIX}*.csv"))
    if not residual_files:
        raise FileNotFoundError(
            f"No bump mesh residual files found in {case_dir}"
        )

    mesh_residuals = []
    for residual_file in residual_files:
        mesh_name = residual_file.stem[len(RESIDUAL_FILE_PREFIX):]
        mesh_fields = mesh_name.split("x")
        if len(mesh_fields) != 2:
            raise ValueError(f"Invalid mesh size in {residual_file.name}")
        try:
            ni, nj = map(int, mesh_fields)
        except ValueError as error:
            raise ValueError(
                f"Invalid mesh size in {residual_file.name}"
            ) from error
        if ni <= 0 or nj <= 0:
            raise ValueError(
                f"Mesh dimensions in {residual_file.name} must be positive"
            )

        residual_data = np.loadtxt(residual_file, ndmin=2)
        if (
            residual_data.shape[0] == 0
            or residual_data.shape[1] != RESIDUAL_COLUMN_COUNT
        ):
            raise ValueError(
                f"{residual_file} must contain nine convergence columns"
            )
        final_density_residual = residual_data[-1, 1]
        if (
            not np.isfinite(final_density_residual)
            or final_density_residual <= 0
        ):
            raise ValueError(
                f"Final average density residual in {residual_file} "
                "must be finite and positive"
            )
        mesh_residuals.append((ni, final_density_residual))

    mesh_residuals.sort(key=lambda result: result[0])
    i_sizes = [result[0] for result in mesh_residuals]
    if len(set(i_sizes)) != len(i_sizes):
        raise ValueError("Bump mesh residual files contain duplicate ni values")
    return (
        np.asarray(i_sizes),
        np.asarray([result[1] for result in mesh_residuals]),
    )


def plot_mesh_residuals(i_sizes, residuals, output_dir):
    figure, axis = plt.subplots(figsize=[8, 6])
    axis.loglog(
        i_sizes,
        residuals,
        marker="o",
        color="tab:blue",
        label="Final average density residual",
    )
    reference_ni = i_sizes[-1]
    reference_residual = residuals[-1]
    for (order, order_name), color in zip(TREND_ORDERS, TREND_COLORS):
        trend_residuals = reference_residual * (i_sizes / reference_ni) ** -order
        axis.loglog(
            i_sizes,
            trend_residuals,
            linestyle="--",
            color=color,
            linewidth=1,
            label=f"{order_name} order trend",
        )
    axis.set_xlabel("Number of i-direction cells ($n_i$)", fontsize=14)
    axis.set_ylabel("Final average density residual ($\\rho$)", fontsize=14)
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    axis.legend(fontsize="small")
    figure.tight_layout()
    plot_file = output_dir / OUTPUT_FILE
    figure.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved plot as {plot_file.name}")


def main():
    case_dir = Path.cwd()
    i_sizes, residuals = read_mesh_residuals(case_dir)
    plot_mesh_residuals(i_sizes, residuals, case_dir)


if __name__ == "__main__":
    main()
