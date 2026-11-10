#!/usr/bin/env python3
#
# Run a case at a range of CFL values, time each solver run, and plot the times.
#
# From the Cases directory, run:
#     python ../Code/plot_time_from_CFL.py bend

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import matplotlib.pyplot as plt
import numpy as np


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
    elapsed_seconds = []
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
                start_time = time.perf_counter()
                subprocess.run(
                    [str(solver)],
                    input=input_text,
                    text=True,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=working_dir,
                    check=True,
                )
                elapsed_seconds.append(time.perf_counter() - start_time)

            residual_source = working_dir / residual_source_name
            if not residual_source.is_file():
                raise FileNotFoundError(
                    f"Solver did not create residual file {residual_source}"
                )
            residual_destination = (
                case_dir / f"conv_{case_name}_CFL_{cfl_text}.csv"
            )
            shutil.copy2(residual_source, residual_destination)

    timing_file = case_dir / f"simulation_times_{case_name}_CFL.csv"
    np.savetxt(
        timing_file,
        np.column_stack((cfl_values, elapsed_seconds)),
        delimiter=",",
        header="CFL,elapsed_seconds",
        comments="",
        fmt=("%.5f", "%.6f"),
    )

    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot(cfl_values, elapsed_seconds, marker="o")
    axis.set_xlabel("CFL number")
    axis.set_ylabel("Simulation time (seconds)")
    axis.set_title(f"Solver time versus CFL ({case_name})")
    axis.grid(True, linestyle="-", linewidth=0.5, alpha=0.6)
    figure.tight_layout()
    plot_file = case_dir / f"simulation_time_{case_name}_CFL.png"
    figure.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Saved residual histories as conv_{case_name}_CFL_*.csv")
    print(f"Saved timing data as {timing_file.name}")
    print(f"Saved timing plot as {plot_file.name}")


if __name__ == "__main__":
    main()
