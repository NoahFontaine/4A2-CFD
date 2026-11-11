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

import matplotlib.pyplot as plt
import numpy as np


def read_cpu_time(time_output):
    time_lines = time_output.read_text().splitlines()
    if len(time_lines) < 2:
        raise ValueError(f"Incomplete time output in {time_output}")
    cpu_fields = time_lines[-2].split()
    if len(cpu_fields) != 2 or cpu_fields[0] != "user":
        raise ValueError(
            f"Expected the second-to-last time output line to be user CPU "
            f"time in {time_output}"
        )
    try:
        cpu_seconds = float(cpu_fields[1])
    except ValueError as error:
        raise ValueError(f"Invalid CPU time in {time_output}") from error
    if not np.isfinite(cpu_seconds) or cpu_seconds <= 0:
        raise ValueError(f"CPU time in {time_output} must be finite and positive")
    return cpu_seconds


def run_case(case_name, parameter_name, parameter_values, case_dir, solver,
             time_command):
    input_file = case_dir / f"input_{case_name}.txt"
    if not input_file.is_file():
        raise FileNotFoundError(f"Case input file not found: {input_file}")

    input_lines = input_file.read_text().splitlines()
    if len(input_lines) < 5:
        raise ValueError(f"{input_file} must contain at least five rows")

    parameter_fields = input_lines[2].split()
    if len(parameter_fields) < 3:
        raise ValueError(
            f"The third row of {input_file} must contain CFL, smoothing, "
            "and convergence values"
        )
    parameter_index = {"CFL": 0, "sfac": 1}[parameter_name]

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
    cpu_times = []
    residual_source_name = f"conv_{case_name}.csv"

    for parameter_value in parameter_values:
        value_text = f"{parameter_value:.8g}"
        solver_inputs = input_lines.copy()
        row_fields = solver_inputs[2].split()
        row_fields[parameter_index] = value_text
        solver_inputs[2] = " ".join(row_fields)
        input_text = "\n".join(solver_inputs) + "\n"

        with tempfile.TemporaryDirectory(
            prefix=f".{case_name}_{parameter_name}_", dir=case_dir
        ) as working_directory:
            working_dir = Path(working_directory)
            shutil.copy2(geometry_input, working_dir / geometry_input.name)
            log_path = (
                case_dir
                / f"output_{case_name}_{parameter_name}_{value_text}.log"
            )
            time_output = working_dir / "time.txt"
            with log_path.open("w") as log_file:
                print(f"Running {case_name} with {parameter_name}={value_text}")
                subprocess.run(
                    [
                        time_command,
                        "-p",
                        "-o",
                        str(time_output),
                        str(solver),
                    ],
                    input=input_text,
                    text=True,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=working_dir,
                    check=True,
                )
                cpu_times.append(read_cpu_time(time_output))

            residual_source = working_dir / residual_source_name
            if not residual_source.is_file():
                raise FileNotFoundError(
                    f"Solver did not create residual file {residual_source}"
                )
            residual_destination = (
                case_dir
                / f"conv_{case_name}_{parameter_name}_{value_text}.csv"
            )
            shutil.copy2(residual_source, residual_destination)

    print(
        f"Saved residual histories as "
        f"conv_{case_name}_{parameter_name}_*.csv"
    )
    return np.asarray(cpu_times)


def plot_times(case_times, parameter_values, parameter_name, case_dir):
    figure, axis = plt.subplots(figsize=[9.6,7.2])
    case_colors = {"bend": "tab:blue", "bump": "tab:red"}
    order_colors = {0: "tab:green", 1: "tab:orange", 2: "tab:purple"}

    for case_name, cpu_times in case_times.items():
        if np.any(~np.isfinite(cpu_times)) or np.any(cpu_times <= 0):
            raise ValueError(f"CPU times for {case_name} must be finite and positive")
        color = case_colors.get(case_name, "tab:blue")
        axis.plot(
            parameter_values,
            cpu_times,
            marker="o",
            color=color,
            label=f"{case_name} CPU time",
        )

        if parameter_name == "CFL":
            reference_time = cpu_times[-1]
            reference_cfl = parameter_values[-1]
            for order, trend_color in order_colors.items():
                trend_times = (
                    reference_time * (parameter_values / reference_cfl) ** order
                )
                axis.plot(
                    parameter_values,
                    trend_times,
                    linestyle="--",
                    color=trend_color,
                    alpha=0.8,
                    label=f"{case_name} {order}th order (CFL^{order})",
                )

    case_names = list(case_times)
    timing_file = (
        case_dir
        / f"simulation_times_{'_'.join(case_names)}_{parameter_name}.csv"
    )
    np.savetxt(
        timing_file,
        np.column_stack([parameter_values, *case_times.values()]),
        delimiter=",",
        header=parameter_name + "," + ",".join(
            f"{name}_cpu_seconds" for name in case_names
        ),
        comments="",
        fmt="%.8g",
    )

    if parameter_name == "CFL":
        axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("CFL number" if parameter_name == "CFL" else "sfac")
    axis.set_ylabel("CPU time (seconds)")
    axis.set_title(f"Solver time versus {parameter_name}")
    axis.tick_params(direction="in", which="both")
    axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    axis.legend(fontsize="small", ncol=2)
    figure.tight_layout()
    plot_file = (
        case_dir
        / f"simulation_time_{'_'.join(case_names)}_{parameter_name}.png"
    )
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
    time_command = shutil.which("time")
    if time_command is None:
        raise FileNotFoundError("The 'time' executable was not found")

    cfl_values = np.geomspace(1e-5, 0.9, 10)
    sfac_values = np.linspace(0.0, 1.0, 11)
    for parameter_name, parameter_values in (
        ("CFL", cfl_values),
        ("sfac", sfac_values),
    ):
        print(f"Starting {parameter_name} sweep")
        case_times = {
            case_name: run_case(
                case_name,
                parameter_name,
                parameter_values,
                case_dir,
                solver,
                time_command,
            )
            for case_name in case_names
        }
        plot_times(case_times, parameter_values, parameter_name, case_dir)


if __name__ == "__main__":
    main()
