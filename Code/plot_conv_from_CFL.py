#!/usr/bin/env python3
#
# Run cases at twenty CFL values and compare final density residuals.
#
# From the Cases directory, run:
#     python ../Code/plot_conv_from_CFL.py
# Optionally specify one or more case names:
#     python ../Code/plot_conv_from_CFL.py bend bump

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np


def plot_residuals(case_residuals, cfl_values, output_dir):
    figure, axis = plt.subplots(figsize=[8, 6])
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    order_colors = {0: "#000000", 1: "#006400", 2: "#4B0082"}

    for (case_name, residuals), color in zip(case_residuals.items(), colors):
        if np.any(~np.isfinite(residuals)) or np.any(residuals <= 0):
            raise ValueError(
                f"Final average density residuals for {case_name} "
                "must be finite and positive"
            )
        axis.plot(
            cfl_values,
            residuals,
            marker="o",
            color=color,
            label=f"{case_name} residual",
        )

        reference_value = residuals[-1]
        reference_cfl = cfl_values[-1]
        if case_name == "bend":
            for order in (0, 1, 2):
                reference_residuals = (
                    reference_value * (cfl_values / reference_cfl) ** order
                )
                axis.plot(
                    cfl_values,
                    reference_residuals,
                    linestyle="--",
                    color=order_colors[order],
                    alpha=0.55,
                    label=f"{order}th order (CFL^{order})",
                )

    case_names = " and ".join(case_residuals)
    axis.set_title(f"Final average density residual versus CFL: {case_names}")
    axis.set_xlabel("CFL number")
    axis.set_ylabel("Residual")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.5)
    axis.legend(fontsize="small")
    figure.tight_layout()
    plot_name = f"convergence_plot_{'_'.join(case_residuals)}_CFL.png"
    figure.savefig(output_dir / plot_name,
                   dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_case(case_name, cfl_values, case_dir, solver):
    input_file = case_dir / f"input_{case_name}.txt"
    if not input_file.is_file():
        raise FileNotFoundError(f"Case input file not found: {input_file}")

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
    solver_inputs = input_lines.copy()
    final_density_residuals = []
    residual_source_name = f"conv_{case_name}.csv"

    for cfl in cfl_values:
        cfl_text = f"{cfl:.8g}"
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
            residual_data = np.loadtxt(residual_destination, ndmin=2)
            if residual_data.shape[1] < 9:
                raise ValueError(
                    f"{residual_destination} must contain nine columns"
                )
            final_density_residuals.append(residual_data[-1, 1])

    print(f"Saved residual histories as conv_{case_name}_CFL_*.csv")
    return np.asarray(final_density_residuals)


def main():
    case_names = sys.argv[1:] or ["bend", "bump"]
    case_dir = Path.cwd()
    solver = Path(__file__).resolve().parent / "solver.x"
    if not solver.is_file():
        raise FileNotFoundError(f"Solver executable not found: {solver}")

    cfl_values = np.geomspace(1e-5, 0.9, 10)
    case_residuals = {
        case_name: run_case(case_name, cfl_values, case_dir, solver)
        for case_name in case_names
    }
    plot_residuals(case_residuals, cfl_values, case_dir)
    print(f"Saved plot as convergence_plot_{'_'.join(case_names)}_CFL.png")


if __name__ == "__main__":
    main()
