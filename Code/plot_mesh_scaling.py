#!/usr/bin/env python3
"""Plot CPU time against i-direction mesh size, colored by iteration count.

Run from the directory containing mesh_timing_bump.csv:
    python ../Code/plot_mesh_scaling.py

Optionally provide a timing CSV and output image path:
    python ../Code/plot_mesh_scaling.py mesh_timing_bump.csv --output scaling.png
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TREND_ORDERS = ((0, "0th"), (1, "1st"), (2, "2nd"), (3, "3rd"))
TREND_COLORS = ("tab:red", "tab:green", "tab:purple", "tab:brown")


def read_mesh_timings(csv_file):
    required_columns = {"ni", "iterations", "cpu_seconds"}
    with csv_file.open(newline="") as timing_stream:
        reader = csv.DictReader(timing_stream)
        if reader.fieldnames is None or not required_columns.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"{csv_file} must contain columns: "
                f"{', '.join(sorted(required_columns))}"
            )

        timings = []
        for row_number, row in enumerate(reader, start=2):
            try:
                ni_value = float(row["ni"])
                iteration_value = float(row["iterations"])
                cpu_seconds = float(row["cpu_seconds"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid numeric data in {csv_file}, row {row_number}"
                ) from error

            if (
                not np.isfinite(ni_value)
                or ni_value <= 0
                or not ni_value.is_integer()
            ):
                raise ValueError(
                    f"ni in {csv_file}, row {row_number} must be a "
                    "positive integer"
                )
            if (
                not np.isfinite(iteration_value)
                or iteration_value <= 0
                or not iteration_value.is_integer()
            ):
                raise ValueError(
                    f"iterations in {csv_file}, row {row_number} must be a "
                    "positive integer"
                )
            if not np.isfinite(cpu_seconds) or cpu_seconds <= 0:
                raise ValueError(
                    f"cpu_seconds in {csv_file}, row {row_number} must be "
                    "finite and positive"
                )
            timings.append((int(ni_value), int(iteration_value), cpu_seconds))

    if not timings:
        raise ValueError(f"{csv_file} contains no mesh timing data")

    timings.sort(key=lambda result: result[0])
    i_sizes = np.asarray([timing[0] for timing in timings])
    if len(np.unique(i_sizes)) != len(i_sizes):
        raise ValueError(f"{csv_file} contains duplicate ni values")
    iterations = np.asarray([timing[1] for timing in timings])
    cpu_times = np.asarray([timing[2] for timing in timings])
    return i_sizes, iterations, cpu_times


def plot_mesh_timings(i_sizes, iterations, cpu_times, output_file):
    figure, time_axis = plt.subplots(figsize=[8, 6])
    time_axis.loglog(
        i_sizes,
        cpu_times,
        color="0.4",
        linewidth=1.5,
        label="CPU time",
        zorder=1,
    )
    reference_ni = i_sizes[-1]
    reference_time = cpu_times[-1]
    for (order, order_name), color in zip(TREND_ORDERS, TREND_COLORS):
        trend_times = reference_time * (i_sizes / reference_ni) ** order
        time_axis.loglog(
            i_sizes,
            trend_times,
            linestyle="--",
            color=color,
            linewidth=1,
            label=f"{order_name} order trend",
            zorder=1,
        )
    markers = time_axis.scatter(
        i_sizes,
        cpu_times,
        c=iterations,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.5,
        zorder=2,
    )
    for ni, count, cpu_time in zip(i_sizes, iterations, cpu_times):
        time_axis.annotate(
            f"{count}",
            (ni, cpu_time),
            xytext=(5, 5),
            textcoords="offset points",
            color="black",
            fontsize="small",
        )

    time_axis.set_xlabel("Number of i-direction cells ($ni$)", fontsize=14)
    time_axis.set_ylabel("CPU time (seconds)", fontsize=14)
    time_axis.tick_params(direction="in", which="both")
    time_axis.grid(True, which="both", linestyle="-", linewidth=0.5, alpha=0.6)
    time_axis.legend(fontsize="small", loc="best")
    figure.colorbar(markers, ax=time_axis, label="Number of iterations")
    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved mesh scaling plot as {output_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "timing_csv",
        nargs="?",
        type=Path,
        default=Path("mesh_timing_bump.csv"),
        help="CSV with ni, iterations, and cpu_seconds columns",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output image path (default: <timing CSV stem>_scaling.png)",
    )
    arguments = parser.parse_args()

    if not arguments.timing_csv.is_file():
        raise FileNotFoundError(
            f"Mesh timing data file not found: {arguments.timing_csv}"
        )
    output_file = arguments.output or arguments.timing_csv.with_name(
        f"{arguments.timing_csv.stem}_scaling.png"
    )
    i_sizes, iterations, cpu_times = read_mesh_timings(arguments.timing_csv)
    plot_mesh_timings(i_sizes, iterations, cpu_times, output_file)


if __name__ == "__main__":
    main()
