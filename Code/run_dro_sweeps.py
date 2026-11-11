#!/usr/bin/env python3
#
# Run five-value CFL and sfac sweeps, saving each density residual history.
#
# From the Cases directory, run:
#     python ../Code/run_dro_sweeps.py
# Choose sweep values with --cfl-values and --sfac-values. For example:
#     python ../Code/run_dro_sweeps.py --cfl-values 0.1 0.2 0.3 0.5 0.8
#     --sfac-values 0.0 0.15 0.3 0.45 0.6
# Optionally specify one or more case names:
#     python ../Code/run_dro_sweeps.py bend bump

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np


DEFAULT_CFL_VALUES = (0.001, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8)
DEFAULT_SFAC_VALUES = (0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 9.9)
PARAMETER_COLUMNS = {"CFL": 0, "sfac": 1}
RESIDUAL_COLUMN_COUNT = 9


def validate_values(parameter_name, values):
    if len(values) != 7:
        raise ValueError(f"Exactly seven {parameter_name} values are required")
    if len(set(values)) != 7:
        raise ValueError(f"The seven {parameter_name} values must be distinct")
    if any(not np.isfinite(value) for value in values):
        raise ValueError(f"{parameter_name} values must be finite")
    if parameter_name == "CFL" and any(value <= 0 for value in values):
        raise ValueError("CFL values must be positive")


def read_case_input(case_name, case_dir):
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
    try:
        default_cfl, default_sfac = map(float, parameter_fields[:2])
    except ValueError as error:
        raise ValueError(f"Invalid CFL or sfac default in {input_file}") from error
    if not np.isfinite(default_cfl) or default_cfl <= 0:
        raise ValueError(f"Default CFL in {input_file} must be finite and positive")
    if not np.isfinite(default_sfac):
        raise ValueError(f"Default sfac in {input_file} must be finite")

    mesh_fields = input_lines[4].split()
    if not mesh_fields:
        raise ValueError(f"The fifth row of {input_file} has no mesh size")
    try:
        mesh_size = int(mesh_fields[0])
    except ValueError as error:
        raise ValueError(f"Invalid mesh size in {input_file}") from error
    geometry_input = (
        case_dir / f"geom_{case_name}.txt"
        if mesh_size != -1
        else case_dir / f"mesh_{case_name}.bin"
    )
    if not geometry_input.is_file():
        raise FileNotFoundError(f"Solver input file not found: {geometry_input}")

    return input_lines, geometry_input, default_cfl, default_sfac


def run_case(case_name, parameter_name, parameter_values, case_dir, solver):
    input_lines, geometry_input, default_cfl, default_sfac = read_case_input(
        case_name, case_dir
    )
    parameter_index = PARAMETER_COLUMNS[parameter_name]
    residual_source_name = f"conv_{case_name}.csv"

    for parameter_value in parameter_values:
        value_text = f"{parameter_value:.8g}"
        row_fields = input_lines[2].split()
        row_fields[parameter_index] = value_text
        input_text = "\n".join(
            [*input_lines[:2], " ".join(row_fields), *input_lines[3:]]
        ) + "\n"

        with tempfile.TemporaryDirectory(
            prefix=f".{case_name}_{parameter_name}_", dir=case_dir
        ) as working_directory:
            working_dir = Path(working_directory)
            shutil.copy2(geometry_input, working_dir / geometry_input.name)
            log_path = (
                case_dir
                / f"output_{case_name}_{parameter_name}_{value_text}.log"
            )
            with log_path.open("w") as log_file:
                print(f"Running {case_name} with {parameter_name}={value_text}")
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
            residual_data = np.loadtxt(residual_source, ndmin=2)
            if (
                residual_data.shape[0] == 0
                or residual_data.shape[1] != RESIDUAL_COLUMN_COUNT
            ):
                raise ValueError(
                    f"{residual_source} must contain nine convergence columns"
                )
            residual_destination = (
                case_dir
                / f"conv_{case_name}_{parameter_name}_{value_text}.csv"
            )
            shutil.copy2(residual_source, residual_destination)
            print(f"Saved residual history as {residual_destination.name}")

    held_parameter = "sfac" if parameter_name == "CFL" else "CFL"
    held_value = default_sfac if parameter_name == "CFL" else default_cfl
    print(
        f"During the {parameter_name} sweep, {held_parameter} used its "
        f"input default ({held_value:g})."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", default=["bend", "bump"])
    parser.add_argument(
        "--cfl-values",
        nargs=7,
        type=float,
        default=DEFAULT_CFL_VALUES,
        metavar=("CFL1", "CFL2", "CFL3", "CFL4", "CFL5", "CFL6", "CFL7"),
    )
    parser.add_argument(
        "--sfac-values",
        nargs=7,
        type=float,
        default=DEFAULT_SFAC_VALUES,
        metavar=("SFAC1", "SFAC2", "SFAC3", "SFAC4", "SFAC5", "SFAC6", "SFAC7"),
    )
    args = parser.parse_args()

    validate_values("CFL", args.cfl_values)
    validate_values("sfac", args.sfac_values)

    case_dir = Path.cwd()
    solver = Path(__file__).resolve().parent / "solver.x"
    if not solver.is_file():
        raise FileNotFoundError(f"Solver executable not found: {solver}")

    for case_name in args.cases:
        for parameter_name, parameter_values in (
            ("CFL", args.cfl_values),
            ("sfac", args.sfac_values),
        ):
            run_case(
                case_name,
                parameter_name,
                parameter_values,
                case_dir,
                solver,
            )


if __name__ == "__main__":
    main()
