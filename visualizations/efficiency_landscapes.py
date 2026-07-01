from typing import Optional, Callable
import math


def _default_objective(h_per_nadh: float, h_per_atp: float) -> float:
    nadh_atp = (10 * h_per_nadh / h_per_atp) * 1.0
    fadh2_atp = (2 * h_per_nadh * 0.67 / h_per_atp) * 1.0
    total = 4 + nadh_atp + fadh2_atp
    return min(total / 94.0, 1.0)


def efficiency_landscape(param_x: str = "h_per_nadh",
                         param_y: str = "h_per_atp",
                         objective_func: Optional[Callable] = None,
                         x_range: tuple = (20.0, 40.0),
                         y_range: tuple = (2.5, 4.0),
                         n_grid: int = 50,
                         save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    func = objective_func or _default_objective

    x_vals = np.linspace(x_range[0], x_range[1], n_grid)
    y_vals = np.linspace(y_range[0], y_range[1], n_grid)
    X, Y = np.meshgrid(x_vals, y_vals)

    Z = np.zeros_like(X)
    for i in range(n_grid):
        for j in range(n_grid):
            Z[i, j] = func(X[i, j], Y[i, j])

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    contour = ax.contourf(X, Y, Z, levels=20, cmap="viridis", alpha=0.9)
    cs = ax.contour(X, Y, Z, levels=10, colors="white", linewidths=0.8, alpha=0.6)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%.2f")

    target_x = 30.0
    target_y = 3.0
    ax.plot(target_x, target_y, "r*", markersize=18, markeredgecolor="black",
            markeredgewidth=0.5, zorder=5)
    ax.annotate("PTR-94 Target\n(30 H⁺/NADH, H⁺/ATP=3)",
                xy=(target_x, target_y), xytext=(target_x + 3, target_y + 0.3),
                fontsize=9, fontweight="bold", color="red",
                arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    cbar = fig.colorbar(contour, ax=ax, label="Normalized ATP Yield\n(fraction of 94)")
    ax.set_xlabel(f"{param_x}", fontsize=11)
    ax.set_ylabel(f"{param_y}", fontsize=11)
    ax.set_title("PTR-94 Efficiency Landscape\n(ATP yield as function of key PCM parameters)",
                 fontsize=12, fontweight="bold")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


def efficiency_surface_3d(param_x: str = "H⁺/NADH",
                          param_y: str = "H⁺/ATP",
                          objective_func: Optional[Callable] = None,
                          x_range: tuple = (20.0, 40.0),
                          y_range: tuple = (2.5, 4.0),
                          n_grid: int = 50,
                          save_path: Optional[str] = None) -> Optional[str]:
    try:
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    func = objective_func or _default_objective

    x_vals = np.linspace(x_range[0], x_range[1], n_grid)
    y_vals = np.linspace(y_range[0], y_range[1], n_grid)
    X, Y = np.meshgrid(x_vals, y_vals)

    Z = np.zeros_like(X)
    for i in range(n_grid):
        for j in range(n_grid):
            Z[i, j] = func(X[i, j], Y[i, j])

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(X, Y, Z, cmap="viridis", alpha=0.9,
                           linewidth=0, antialiased=True)

    target_x = 30.0
    target_y = 3.0
    target_z = func(target_x, target_y)
    ax.scatter([target_x], [target_y], [target_z], color="red", s=100,
               marker="*", edgecolors="black", linewidth=0.5, zorder=10)

    ax.set_xlabel(param_x, fontsize=10, labelpad=10)
    ax.set_ylabel(param_y, fontsize=10, labelpad=10)
    ax.set_zlabel("ATP Yield Fraction", fontsize=10, labelpad=10)
    ax.set_title("PTR-94 Efficiency Surface (3D)", fontsize=13, fontweight="bold")

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Normalized ATP Yield")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    r1 = efficiency_landscape(save_path="/tmp/ptr94_landscape.png")
    print(f"Efficiency landscape saved to {r1}" if r1 else "Matplotlib not available")
    r2 = efficiency_surface_3d(save_path="/tmp/ptr94_surface.png")
    print(f"Efficiency surface saved to {r2}" if r2 else "mpl_toolkits/matplotlib not available")
