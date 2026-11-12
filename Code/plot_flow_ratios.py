#!/usr/bin/env python3
"""Plot mass-averaged stagnation-property and mass-flow ratios by i-index.

Run from the directory containing the case files:
    python ../Code/plot_flow_ratios.py bend bump

Optionally choose an output image path:
    python ../Code/plot_flow_ratios.py bend bump --output flow_ratios.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from routines import calc_secondary, cut_i, face_av, mass_av, read_case, read_settings


def calculate_flow_ratios(settings, flowfield):
    mach_squared = flowfield["mach"] ** 2
    stagnation_factor = 1 + 0.5 * (settings["gam"] - 1) * mach_squared
    flowfield["tstag"] = flowfield["t"] * stagnation_factor
    flowfield["pstag"] = flowfield["p"] * stagnation_factor ** (
        settings["gam"] / (settings["gam"] - 1)
    )

    averages = []
    for i_index in range(flowfield["ni"]):
        cut = cut_i(flowfield, i_index)
        mass_flux = (
            face_av(cut["rovx"]) * cut["lx"]
            + face_av(cut["rovy"]) * cut["ly"]
        )
        mass_flow = np.sum(mass_flux)
        if not np.isfinite(mass_flow) or mass_flow == 0:
            raise ValueError(
                f"Mass flow at i-index {i_index} must be finite and nonzero"
            )

        tstag_average, _ = mass_av(cut, "tstag")
        pstag_average, _ = mass_av(cut, "pstag")
        averages.append((tstag_average, pstag_average, mass_flow))

    averages = np.asarray(averages)
    if not np.all(np.isfinite(averages)):
        raise ValueError("Mass-averaged flow properties must all be finite")
    if np.any(averages[0] == 0):
        raise ValueError("Inlet reference values must be nonzero")

    ratios = averages / averages[0]
    return np.arange(flowfield["ni"]), ratios


def plot_flow_ratios(case_profiles, output_file):
    figure, axis = plt.subplots(figsize=[8, 6])
    labels = (
        r"$T_{0}/T_{0,\mathrm{in}}$",
        r"$p_{0}/p_{0,\mathrm{in}}$",
        r"$\dot{m}/\dot{m}_{\mathrm{in}}$",
    )
    colors = ("blue", "purple", "black")
    for column, (label, color) in enumerate(zip(labels, colors)):
        for case_name, (i_indices, ratios) in case_profiles.items():
            axis.plot(
                i_indices,
                ratios[:, column],
                color=color,
                marker="." if case_name == "bump" else "",
                markersize=3,
                linewidth=1,
                linestyle=":" if case_name == "bend" else "-",
                label=label if case_name == "bump" else "_nolegend_",
            )

    axis.set_xlabel("Cell index", fontsize=14)
    axis.set_ylabel("Ratio to inlet value", fontsize=14)
    axis.set_ylim(0.985, 1.01)
    axis.tick_params(direction="in")
    axis.grid(True, linestyle="-", linewidth=0.5, alpha=0.6)
    if "bump" in case_profiles:
        axis.legend(fontsize="small")
    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved flow-ratio plot as {output_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "case_names",
        nargs="+",
        help="one or more case names, such as bend bump",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output image path (default: flow_ratios_<cases>.png)",
    )
    arguments = parser.parse_args()

    if len(set(arguments.case_names)) != len(arguments.case_names):
        parser.error("case names must be unique")

    case_profiles = {}
    for case_name in arguments.case_names:
        input_file = Path(f"input_{case_name}.txt")
        solution_file = Path(f"out_final_{case_name}.bin")
        if not input_file.is_file():
            raise FileNotFoundError(f"Case input file not found: {input_file}")
        if not solution_file.is_file():
            raise FileNotFoundError(
                f"Final solution file not found: {solution_file}"
            )

        settings = read_settings(str(input_file))
        flowfield = calc_secondary(settings, read_case(str(solution_file)))
        case_profiles[case_name] = calculate_flow_ratios(settings, flowfield)

    output_file = arguments.output or Path(
        f"flow_ratios_{'_'.join(arguments.case_names)}.png"
    )
    plot_flow_ratios(case_profiles, output_file)


if __name__ == "__main__":
    main()
