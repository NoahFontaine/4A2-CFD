#
#   plot_density_residual
#
#   Plot the cellwise density residual from a converged 4A2 solver run.
#
#   From the directory containing the output files, run:
#       python path_to_script/plot_density_residual.py casename

from routines import *


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python plot_density_residual.py casename"
        )

    case_name = sys.argv[1]
    g = read_case(f"out_final_{case_name}.bin")

    rho_ref = np.mean(g["ro"][0, :])
    if not np.isfinite(rho_ref) or rho_ref <= 0:
        raise ValueError("Reference density must be finite and positive")

    density_residual = np.abs(g["dro"])

    fig = plt.figure(figsize=[9.6, 7.2])
    ax = plt.axes()
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    hc = ax.pcolormesh(
        g["x"],
        g["y"],
        density_residual,
        vmin=0.0,
        vmax=0.0001,
        shading="flat",
    )
    x_cell = 0.25 * (
        g["x"][:-1, :-1]
        + g["x"][1:, :-1]
        + g["x"][:-1, 1:]
        + g["x"][1:, 1:]
    )
    y_cell = 0.25 * (
        g["y"][:-1, :-1]
        + g["y"][1:, :-1]
        + g["y"][:-1, 1:]
        + g["y"][1:, 1:]
    )
    ax.contour(
        x_cell,
        y_cell,
        density_residual,
        colors="white",
        linewidths=0.5,
    )
    colorbar(hc, r"Density residual ($|\Delta \rho|$)")
    plot_wall(ax, g)

    plt.savefig(
        f"contour_plots_{case_name}_density_residual.png",
        dpi=300,
        bbox_inches="tight",
    )


main()
