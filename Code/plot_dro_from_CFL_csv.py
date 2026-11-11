#!/usr/bin/env python3
#
# Plot average and maximum density residuals for saved sweeps.
#
# From the directory containing the residual CSV files, run:
#     python ../Code/plot_dro_from_CFL_csv.py
# Optionally specify one or more case names:
#     python ../Code/plot_dro_from_CFL_csv.py bend bump

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
import numpy as np

RESIDUAL_COLUMN_COUNT = 9


def read_density_histories(case_name, parameter_name, case_dir):
    prefix = f"conv_{case_name}_{parameter_name}_"
    residual_files = list(case_dir.glob(f"{prefix}*.csv"))
    if not residual_files:
        raise FileNotFoundError(
            f"No saved {parameter_name} residual files found for "
            f"{case_name} in {case_dir}"
        )

    histories = []
    for residual_file in residual_files:
        try:
            parameter_value = float(residual_file.stem[len(prefix):])
        except ValueError as error:
            raise ValueError(
                f"Could not read {parameter_name} value from "
                f"{residual_file.name}"
            ) from error
        if not np.isfinite(parameter_value):
            raise ValueError(
                f"{parameter_name} value in {residual_file.name} must be finite"
            )
        if parameter_name == "CFL" and parameter_value <= 0:
            raise ValueError(f"CFL value in {residual_file.name} must be positive")

        data = np.loadtxt(residual_file, ndmin=2)
        if data.shape[0] == 0 or data.shape[1] != RESIDUAL_COLUMN_COUNT:
            raise ValueError(
                f"{residual_file} must contain nine convergence columns"
            )
        iterations = data[:, 0]
        dro_avg = data[:, 1]
        dro_max = data[:, 5]
        if (
            np.any(~np.isfinite(np.column_stack([iterations, dro_avg, dro_max])))
            or np.any(dro_avg <= 0)
            or np.any(dro_max <= 0)
        ):
            raise ValueError(
                f"Iteration, dro_avg and dro_max in {residual_file} "
                "must be finite, and residuals must be positive"
            )
        histories.append((parameter_value, iterations, dro_avg, dro_max))

    return sorted(histories, key=lambda history: history[0])


def plot_density_histories(case_name, parameter_name, histories, output_dir):
    parameter_values = [history[0] for history in histories]
    minimum = min(parameter_values)
    maximum = max(parameter_values)
    if parameter_name == "CFL":
        norm = LogNorm(
            vmin=minimum / 10 if minimum == maximum else minimum,
            vmax=maximum * 10 if minimum == maximum else maximum,
        )
    else:
        norm = Normalize(
            vmin=minimum - 0.5 if minimum == maximum else minimum,
            vmax=maximum + 0.5 if minimum == maximum else maximum,
        )

    colormap = plt.get_cmap("viridis")
    figure, axis = plt.subplots(figsize=[8, 6])
    for parameter_value, iterations, dro_avg, dro_max in histories:
        color = colormap(norm(parameter_value))
        axis.plot(iterations, dro_avg, color=color, linestyle="-")
        axis.plot(iterations, dro_max, color=color, linestyle="--")

    axis.set_xlabel("Iteration", fontsize=14)
    axis.set_ylabel("Density residual ($\\rho$)", fontsize=14)
    # axis.set_title(
    #     f"Density residual history: {case_name} ({parameter_name} sweep)"
    # )
    axis.set_yscale("log")
    axis.set_xlim(left=0, right=10000)
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    axis.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="-", label="Average"),
            Line2D([0], [0], color="black", linestyle="--", label="Maximum"),
        ],
        fontsize="small",
    )
    colorbar = figure.colorbar(
        ScalarMappable(norm=norm, cmap=colormap), ax=axis
    )
    colorbar.set_label("CFL number" if parameter_name == "CFL" else "sfac")
    figure.tight_layout()
    plot_file = (
        output_dir
        / f"convergence_plot_{case_name}_{parameter_name}_dro.png"
    )
    figure.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved plot as {plot_file.name}")


def plot_final_residuals_vs_sfac(case_histories, output_dir):
    figure, axis = plt.subplots(figsize=[8, 6])
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for (case_name, histories), color in zip(case_histories.items(), colors):
        sfac_values = [history[0] for history in histories]
        dro_avg = [history[2][-1] for history in histories]
        dro_max = [history[3][-1] for history in histories]
        axis.plot(
            sfac_values,
            dro_avg,
            marker="o",
            markersize=3,
            color=color,
            label=f"{case_name} average",
        )
        axis.plot(
            sfac_values,
            dro_max,
            marker="o",
            markersize=3,
            color=color,
            linestyle="--",
            label=f"{case_name} maximum",
        )

    case_names = "_".join(case_histories)
    axis.set_xlabel("Smoothing factor (sfac)", fontsize=14)
    axis.set_ylabel("Final density residual ($\\rho$)", fontsize=14)
    axis.set_yscale("log")
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.5)
    axis.legend(fontsize="small")
    figure.tight_layout()
    plot_file = output_dir / f"convergence_plot_{case_names}_sfac_final_dro.png"
    figure.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved plot as {plot_file.name}")


def main():
    case_names = sys.argv[1:] or ["bend", "bump"]
    case_dir = Path.cwd()
    sfac_histories = {}
    for case_name in case_names:
        for parameter_name in ("CFL", "sfac"):
            histories = read_density_histories(
                case_name, parameter_name, case_dir
            )
            plot_density_histories(
                case_name, parameter_name, histories, case_dir
            )
            if parameter_name == "sfac":
                sfac_histories[case_name] = histories
    plot_final_residuals_vs_sfac(sfac_histories, case_dir)


if __name__ == "__main__":
    main()
