#!/usr/bin/env python3
#
# Run cases at a range of CFL values, time each solver run, and plot the times.
#
# From the Cases directory, run:
#     python ../Code/plot_time_from_CFL.py

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import matplotlib.pyplot as plt
import numpy as np


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
    elapsed_seconds = []
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
            log_path = case_dir / f"output_{case_name}_CFL_{cfl_text}.log"
            with log_path.open("w") as log_file:
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

    print(f"Saved residual histories as conv_{case_name}_CFL_*.csv")
    return np.asarray(elapsed_seconds)


def plot_times(case_times, cfl_values, case_dir):
    figure, axis = plt.subplots(figsize=(9.6, 7.2))
    case_colors = {"bend": "tab:blue", "bump": "tab:red"}
    order_colors = {0: "tab:green", 1: "tab:orange", 2: "tab:purple"}

    for case_name, elapsed_seconds in case_times.items():
        if np.any(~np.isfinite(elapsed_seconds)) or np.any(elapsed_seconds <= 0):
            raise ValueError(f"Solver times for {case_name} must be finite and positive")
        color = case_colors.get(case_name)
        axis.plot(
            cfl_values,
            elapsed_seconds,
            marker="o",
            color=color,
            label=f"{case_name} simulation time",
        )

        reference_time = elapsed_seconds[-1]
        reference_cfl = cfl_values[-1]
        for order, trend_color in order_colors.items():
            trend_times = reference_time * (cfl_values / reference_cfl) ** order
            axis.plot(
                cfl_values,
                trend_times,
                linestyle="--",
                color=trend_color,
                alpha=0.8,
                label=f"{case_name} {order}th order (CFL^{order})",
            )

    case_names = list(case_times)
    timing_file = case_dir / f"simulation_times_{'_'.join(case_names)}_CFL.csv"
    np.savetxt(
        timing_file,
        np.column_stack([cfl_values, *case_times.values()]),
        delimiter=",",
        header="CFL," + ",".join(f"{name}_seconds" for name in case_names),
        comments="",
        fmt="%.8g",
    )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("CFL number")
    axis.set_ylabel("Simulation time (seconds)")
    axis.set_title("Solver time versus CFL")
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    axis.legend(fontsize="small", ncol=2)
    figure.tight_layout()
    plot_file = case_dir / f"simulation_time_{'_'.join(case_names)}_CFL.png"
    figure.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(figure)

    print(f"Saved timing data as {timing_file.name}")
    print(f"Saved timing plot as {plot_file.name}")


def main():
    case_names = sys.argv[1:] or ["bend", "bump"]
    case_dir = Path.cwd()
    solver = Path(__file__).resolve().parent / "solver.x"
    if not solver.is_file():
        raise FileNotFoundError(f"Solver executable not found: {solver}")

    cfl_values = np.geomspace(1e-5, 0.9, 10)
    case_times = {
        case_name: run_case(case_name, cfl_values, case_dir, solver)
        for case_name in case_names
    }
    plot_times(case_times, cfl_values, case_dir)


if __name__ == "__main__":
    main()
