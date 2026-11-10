#!/usr/bin/env python3
#
# Run a case at ten CFL values and plot the convergence histories.
#
# From the Cases directory, run:
#     python ../Code/plot_conv_from_CFL.py bend

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np


def plot_residuals(case_name, cfl_values, residual_files, output_dir):
    fieldnames = [
        ["dro_avg", "droe_avg", "drovx_avg", "drovy_avg"],
        ["dro_max", "droe_max", "drovx_max", "drovy_max"],
    ]
    residual_columns = {
        "dro_avg": 1,
        "droe_avg": 2,
        "drovx_avg": 3,
        "drovy_avg": 4,
        "dro_max": 5,
        "droe_max": 6,
        "drovx_max": 7,
        "drovy_max": 8,
    }
    labels = ["Density", "Energy", "X momentum", "Y momentum"]
    figure, axes = plt.subplots(2, 4, figsize=(18, 9), sharex=True)
    colors = plt.cm.viridis(np.linspace(0, 1, len(cfl_values)))
    residual_data = [np.loadtxt(filename, ndmin=2) for filename in residual_files]
    if any(data.shape[1] < 9 for data in residual_data):
        raise ValueError("A convergence CSV must contain nine columns")

    for row, fields in enumerate(fieldnames):
        for column, fieldname in enumerate(fields):
            axis = axes[row, column]
            for cfl, data, color in zip(cfl_values, residual_data, colors):
                axis.plot(
                    data[:, 0],
                    data[:, residual_columns[fieldname]],
                    color=color,
                    label=f"CFL={cfl:.5f}",
                )

            statistic = "average" if row == 0 else "maximum"
            axis.set_title(f"{labels[column]} ({statistic})")
            axis.set_yscale("log")
            axis.tick_params(direction="in", which="both")
            axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.5)
            if column == 0:
                axis.set_ylabel("Residual")
            if row == 1:
                axis.set_xlabel("Iteration")

    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, legend_labels, loc="lower center", ncol=5)
    figure.tight_layout(rect=(0, 0.1, 1, 1))
    figure.savefig(
        output_dir / f"convergence_plot_{case_name}_CFL.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            f"Usage: python {Path(sys.argv[0]).name} <casename>"
        )

    case_name = sys.argv[1]
    case_dir = Path.cwd()
    input_file = case_dir / f"input_{case_name}.txt"
    solver = Path(__file__).resolve().parent / "solver.x"

    if not input_file.is_file():
        raise FileNotFoundError(f"Case input file not found: {input_file}")
    if not solver.is_file():
        raise FileNotFoundError(f"Solver executable not found: {solver}")

    input_lines = input_file.read_text().splitlines()
    if len(input_lines) < 5:
        raise ValueError(f"{input_file} must contain at least five rows")

    cfl_fields = input_lines[2].split()
    if len(cfl_fields) < 3:
        raise ValueError(
            f"The third row of {input_file} must contain CFL, smoothing, "
            "and convergence values"
        )

    mesh_fields = input_lines[4].split()
    if not mesh_fields:
        raise ValueError(f"The fifth row of {input_file} has no mesh size")
    geometry_input = (
        case_dir / f"geom_{case_name}.txt"
        if int(mesh_fields[0]) != -1
        else case_dir / f"mesh_{case_name}.bin"
    )
    if not geometry_input.is_file():
        raise FileNotFoundError(f"Solver input file not found: {geometry_input}")

    solver_inputs = input_lines.copy()
    cfl_values = np.linspace(0.1, 1.1, 11)
    residual_files = []
    residual_source_name = f"conv_{case_name}.csv"

    for cfl in cfl_values:
        cfl_text = f"{cfl:.5f}"
        row_fields = solver_inputs[2].split()
        row_fields[0] = cfl_text
        solver_inputs[2] = " ".join(row_fields)
        input_text = "\n".join(solver_inputs) + "\n"

        with tempfile.TemporaryDirectory(
            prefix=f".{case_name}_cfl_", dir=case_dir
        ) as working_directory:
            working_dir = Path(working_directory)
            shutil.copy2(geometry_input, working_dir / geometry_input.name)
            with (case_dir / "output.log").open("w") as log_file:
                print(f"Running {case_name} with CFL={cfl_text}")
                subprocess.run(
                    [str(solver)],
                    input=input_text,
                    text=True,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=working_dir,
                    check=True,
                )

            residual_source = working_dir / residual_source_name
            if not residual_source.is_file():
                raise FileNotFoundError(
                    f"Solver did not create residual file {residual_source}"
                )
            residual_destination = (
                case_dir / f"conv_{case_name}_CFL_{cfl_text}.csv"
            )
            shutil.copy2(residual_source, residual_destination)
            residual_files.append(residual_destination)

    plot_residuals(case_name, cfl_values, residual_files, case_dir)
    print(f"Saved residual histories as conv_{case_name}_CFL_*.csv")
    print(f"Saved plot as convergence_plot_{case_name}_CFL.png")


if __name__ == "__main__":
    main()
