#!/usr/bin/env python3
#
# Plot saved solver timing data without running the solver.
#
# From the directory containing the timing CSV files, run:
#     python ../Code/plot_time_from_CFL_csv.py
# Optionally specify one or more case names:
#     python ../Code/plot_time_from_CFL_csv.py bend bump

from pathlib import Path
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np


def read_timing_data(csv_file, parameter_name, case_names):
    with csv_file.open(newline="") as timing_file:
        reader = csv.DictReader(timing_file)
        if reader.fieldnames is None or len(reader.fieldnames) < 2:
            raise ValueError(f"{csv_file} must have a header and timing columns")
        parameter_column = reader.fieldnames[0]
        if parameter_column.casefold() != parameter_name.casefold():
            raise ValueError(
                f"The first column in {csv_file} must be {parameter_name}"
            )

        available_columns = set(reader.fieldnames[1:])
        case_columns = {}
        for case_name in case_names:
            cpu_column = f"{case_name}_cpu_seconds"
            legacy_column = f"{case_name}_seconds"
            if cpu_column in available_columns:
                case_columns[case_name] = cpu_column
            elif legacy_column in available_columns:
                case_columns[case_name] = legacy_column
            elif len(case_names) == 1 and "elapsed_seconds" in available_columns:
                case_columns[case_name] = "elapsed_seconds"
            else:
                raise ValueError(
                    f"{csv_file} does not contain a timing column for "
                    f"{case_name}"
                )

        parameter_values = []
        case_times = {case_name: [] for case_name in case_names}
        for row_number, row in enumerate(reader, start=2):
            try:
                parameter_values.append(float(row[parameter_column]))
                for case_name, column in case_columns.items():
                    case_times[case_name].append(float(row[column]))
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid numeric data in {csv_file}, row {row_number}"
                ) from error

    if not parameter_values:
        raise ValueError(f"{csv_file} contains no timing data")

    parameter_values = np.asarray(parameter_values)
    if np.any(~np.isfinite(parameter_values)):
        raise ValueError(f"Parameter values in {csv_file} must be finite")
    if parameter_name == "CFL" and np.any(parameter_values <= 0):
        raise ValueError(f"CFL values in {csv_file} must be positive")

    case_times = {
        case_name: np.asarray(elapsed_seconds)
        for case_name, elapsed_seconds in case_times.items()
    }
    for case_name, cpu_times in case_times.items():
        if (
            np.any(~np.isfinite(cpu_times))
            or np.any(cpu_times <= 0)
        ):
            raise ValueError(
                f"CPU times for {case_name} in {csv_file} "
                "must be finite and positive"
            )
    return parameter_values, case_times


def plot_times(case_times, parameter_values, parameter_name, case_dir):
    figure, axis = plt.subplots(figsize=[8, 6])
    case_colors = {"bend": "tab:blue", "bump": "tab:orange"}
    order_colors = {0: "#000000", 1: "#006400", 2: "#4B0082"}

    for case_name, cpu_times in case_times.items():
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
                if case_name == "bend":
                    trend_times = (
                        reference_time * (parameter_values / reference_cfl) ** order
                    )
                    axis.plot(
                        parameter_values,
                        trend_times,
                        linestyle="--",
                        color=trend_color,
                        alpha=0.8,
                        label=f"{order}th order ($CFL^{order}$)",
                    )

    case_names = list(case_times)
    if parameter_name == "CFL":
        axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("CFL number" if parameter_name == "CFL" else "sfac", fontsize=14)
    axis.set_ylabel("CPU time (seconds)", fontsize=14)
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
    print(f"Saved timing plot as {plot_file.name}")


def main():
    case_names = sys.argv[1:] or ["bend", "bump"]
    case_dir = Path.cwd()
    for parameter_name in ("CFL", "sfac"):
        timing_file = (
            case_dir
            / f"simulation_times_{'_'.join(case_names)}_{parameter_name}.csv"
        )
        if not timing_file.is_file():
            raise FileNotFoundError(f"Timing data file not found: {timing_file}")
        parameter_values, case_times = read_timing_data(
            timing_file, parameter_name, case_names
        )
        plot_times(case_times, parameter_values, parameter_name, case_dir)


if __name__ == "__main__":
    main()
