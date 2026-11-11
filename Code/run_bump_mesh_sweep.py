#!/usr/bin/env python3
"""Run the bump case at its base mesh and two successive doubled meshes.

Run from the Cases directory:
    python ../Code/run_bump_mesh_sweep.py

After the sweep, plot the final average density residual against ni with:
    python ../Code/plot_bump_mesh_convergence.py

Plot CPU time and iteration count against ni with:
    python ../Code/plot_mesh_scaling.py

The script saves one residual history per mesh as conv_bump_mesh_NIxNJ.csv,
one solver log per mesh, and iteration counts and CPU timings in
mesh_timing_bump.csv.
"""

import csv
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np


RESIDUAL_COLUMN_COUNT = 9


def read_cpu_time(time_output):
    time_lines = time_output.read_text().splitlines()
    if len(time_lines) < 2:
        raise ValueError(f"Incomplete time output in {time_output}")
    cpu_fields = time_lines[1].split()
    if len(cpu_fields) != 2 or cpu_fields[0] != "user":
        raise ValueError(
            f"Expected the second time output line to be user CPU time "
            f"in {time_output}"
        )
    try:
        cpu_seconds = float(cpu_fields[1])
    except ValueError as error:
        raise ValueError(f"Invalid CPU time in {time_output}") from error
    if not np.isfinite(cpu_seconds) or cpu_seconds < 0:
        raise ValueError(f"CPU time in {time_output} must be finite and non-negative")
    return cpu_seconds


def read_iteration_count(log_file, residual_data):
    for line in log_file.read_text().splitlines():
        match = re.search(
            r"Calculation completed after\s+(\d+)\s+iterations", line
        )
        if match:
            return int(match.group(1))

    last_iteration = residual_data[-1, 0]
    if (
        not np.isfinite(last_iteration)
        or last_iteration < 0
        or not last_iteration.is_integer()
    ):
        raise ValueError(
            f"Could not determine the iteration count from {log_file} "
            f"or its residual history"
        )
    return int(last_iteration)


def read_bump_input(input_file):
    if not input_file.is_file():
        raise FileNotFoundError(f"Bump input file not found: {input_file}")

    input_lines = input_file.read_text().splitlines()
    if len(input_lines) < 5 or input_lines[0].strip() != "bump":
        raise ValueError(f"{input_file} is not a valid bump case input")

    mesh_fields = input_lines[4].split()
    if len(mesh_fields) != 2:
        raise ValueError(f"The fifth row of {input_file} must contain ni and nj")
    try:
        ni, nj = map(int, mesh_fields)
    except ValueError as error:
        raise ValueError(f"Invalid mesh size in {input_file}") from error
    if ni <= 0 or nj <= 0:
        raise ValueError(f"Mesh dimensions in {input_file} must be positive")

    return input_lines, ni, nj


def run_mesh_sweep(case_dir, solver, time_command):
    input_file = case_dir / "input_bump.txt"
    input_lines, base_ni, base_nj = read_bump_input(input_file)
    geometry_file = case_dir / "geom_bump.txt"
    if not geometry_file.is_file():
        raise FileNotFoundError(f"Bump geometry file not found: {geometry_file}")

    timing_file = case_dir / "mesh_timing_bump.csv"
    with timing_file.open("w", newline="") as timing_stream:
        timing_writer = csv.writer(timing_stream)
        timing_writer.writerow(("ni", "nj", "iterations", "cpu_seconds"))

        for multiplier in (1, 2, 4, 6, 8, 12, 16):
            ni = base_ni * multiplier
            nj = base_nj * multiplier
            mesh_name = f"{ni}x{nj}"
            mesh_lines = input_lines.copy()
            mesh_lines[4] = f"{ni} {nj}"
            input_text = "\n".join(mesh_lines) + "\n"

            with tempfile.TemporaryDirectory(
                prefix=f".bump_mesh_{mesh_name}_", dir=case_dir
            ) as working_directory:
                working_dir = Path(working_directory)
                shutil.copy2(geometry_file, working_dir / geometry_file.name)
                time_output = working_dir / "time.txt"
                log_file = case_dir / f"output_bump_mesh_{mesh_name}.log"
                with log_file.open("w") as log_stream:
                    print(f"Running bump mesh {mesh_name}")
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
                        stdout=log_stream,
                        stderr=subprocess.STDOUT,
                        cwd=working_dir,
                        check=True,
                    )

                cpu_seconds = read_cpu_time(time_output)
                residual_source = working_dir / "conv_bump.csv"
                if not residual_source.is_file():
                    raise FileNotFoundError(
                        f"Solver did not create residual file {residual_source}"
                    )
                residual_data = np.loadtxt(residual_source, ndmin=2)
                if (
                    residual_data.shape[0] == 0
                    or residual_data.shape[1] != RESIDUAL_COLUMN_COUNT
                ):
                    raise ValueError(
                        f"{residual_source} must contain nine convergence columns"
                    )

                iterations = read_iteration_count(log_file, residual_data)
                residual_destination = (
                    case_dir / f"conv_bump_mesh_{mesh_name}.csv"
                )
                shutil.copy2(residual_source, residual_destination)
                timing_writer.writerow(
                    (ni, nj, iterations, f"{cpu_seconds:.8g}")
                )
                timing_stream.flush()
                print(
                    f"Saved {residual_destination.name}; "
                    f"{iterations} iterations; "
                    f"user CPU time: {cpu_seconds:.6g} s"
                )

    print(f"Saved CPU timings as {timing_file.name}")


def main():
    case_dir = Path.cwd()
    solver = Path(__file__).resolve().with_name("solver.x")
    if not solver.is_file():
        raise FileNotFoundError(f"Solver executable not found: {solver}")
    time_command = shutil.which("time")
    if time_command is None:
        raise FileNotFoundError("The 'time' executable was not found")

    run_mesh_sweep(case_dir, solver, time_command)


if __name__ == "__main__":
    main()
