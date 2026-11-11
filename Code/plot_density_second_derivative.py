#
#   plot_density_second_derivative
#
#   Plot max(|d2rho/dx2|, |d2rho/dy2|) from a converged 4A2 solver run.
#
#   From the directory containing the output files, run:
#       python path_to_script/plot_density_second_derivative.py casename

from routines import *


def second_derivative_map(x, y, density):
    if min(x.shape) < 3:
        raise ValueError("The mesh must have at least three nodes in each direction")
    if x.shape != y.shape or x.shape != density.shape:
        raise ValueError("Coordinates and density must have matching shapes")
    if not all(np.all(np.isfinite(field)) for field in (x, y, density)):
        raise ValueError("Coordinates and density must contain only finite values")

    x_xi, x_eta = np.gradient(x, axis=(0, 1), edge_order=2)
    y_xi, y_eta = np.gradient(y, axis=(0, 1), edge_order=2)
    rho_xi, rho_eta = np.gradient(density, axis=(0, 1), edge_order=2)

    jacobian = x_xi * y_eta - x_eta * y_xi
    if np.any(~np.isfinite(jacobian)) or np.any(jacobian == 0):
        raise ValueError("The mesh has a singular or invalid coordinate mapping")

    rho_x = (y_eta * rho_xi - y_xi * rho_eta) / jacobian
    rho_y = (-x_eta * rho_xi + x_xi * rho_eta) / jacobian

    rho_x_xi, rho_x_eta = np.gradient(rho_x, axis=(0, 1), edge_order=2)
    rho_y_xi, rho_y_eta = np.gradient(rho_y, axis=(0, 1), edge_order=2)
    rho_xx = (y_eta * rho_x_xi - y_xi * rho_x_eta) / jacobian
    rho_yy = (-x_eta * rho_y_xi + x_xi * rho_y_eta) / jacobian

    return np.maximum(np.abs(rho_xx), np.abs(rho_yy))


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python plot_density_second_derivative.py casename"
        )

    case_name = sys.argv[1]
    g = read_case(f"out_final_{case_name}.bin")
    derivative = second_derivative_map(g["x"], g["y"], g["ro"])

    fig = plt.figure(figsize=[9.6, 7.2])
    ax = plt.axes()
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    hc = ax.pcolormesh(
        g["x"],
        g["y"],
        derivative,
        shading="gouraud",
        vmax=40,
    )
    # positive_values = derivative[derivative > 0]
    # if positive_values.size:
    #     min_value = np.min(positive_values)
    #     max_value = np.max(positive_values)
    #     if min_value < max_value:
    #         contour_levels = np.geomspace((min_value+max_value)/200, max_value, num=4)[1:-1]
    #         ax.contour(
    #             g["x"],
    #             g["y"],
    #             derivative,
    #             levels=contour_levels,
    #             colors="white",
    #             linewidths=0.5,
    #         )
    colorbar(hc, r"$\max(|\frac{\partial^2\rho}{\partial x^2}|, |\frac{\partial^2\rho}{\partial y^2}|)$")
    plot_wall(ax, g)

    plt.savefig(
        f"contour_plots_{case_name}_density_second_derivative.png",
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    main()
