from typing import Optional, List, Dict, Tuple
import random


DEFAULT_PARETO_RESULTS = {
    "points": [
        (80.0, 0.85),
        (85.0, 0.88),
        (88.0, 0.90),
        (90.0, 0.92),
        (91.0, 0.93),
        (92.0, 0.94),
        (93.0, 0.95),
        (93.5, 0.96),
        (94.0, 0.98),
        (94.0, 0.99),
        (94.0, 1.0),
    ],
    "objectives": ["ATP Yield", "Coupling Efficiency"],
}

DEFAULT_GENERATIONS = {
    "gen1": [(70, 0.7), (75, 0.72), (78, 0.74), (80, 0.75)],
    "gen5": [(78, 0.78), (82, 0.82), (85, 0.85), (87, 0.86)],
    "gen10": [(85, 0.88), (88, 0.91), (90, 0.93), (92, 0.94)],
    "gen20": [(90, 0.94), (92, 0.96), (93, 0.97), (94, 0.99)],
}


def _compute_pareto_front(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    pareto = []
    for i, (x1, y1) in enumerate(points):
        dominated = False
        for j, (x2, y2) in enumerate(points):
            if i != j and x2 >= x1 and y2 >= y1 and (x2 > x1 or y2 > y1):
                dominated = True
                break
        if not dominated:
            pareto.append((x1, y1))
    pareto.sort()
    return pareto


def plot_pareto_front(pareto_results: Optional[dict] = None,
                      objectives: Optional[List[str]] = None,
                      save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    data = pareto_results or DEFAULT_PARETO_RESULTS
    obj_names = objectives or data.get("objectives", ["Objective 1", "Objective 2"])
    points = data["points"]

    pareto = _compute_pareto_front(points)

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    all_x = [p[0] for p in points]
    all_y = [p[1] for p in points]
    ax.scatter(all_x, all_y, c="#4477AA", s=60, alpha=0.5,
               edgecolors="black", linewidth=0.5, label="Population")

    pareto_x = [p[0] for p in pareto]
    pareto_y = [p[1] for p in pareto]
    ax.plot(pareto_x, pareto_y, color="#CC3311", linewidth=2.5,
            marker="o", markersize=8, label="Pareto Front", zorder=5)

    ax.set_xlabel(obj_names[0], fontsize=11)
    ax.set_ylabel(obj_names[1], fontsize=11)
    ax.set_title("PTR-94 Pareto Front", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(True, alpha=0.3)

    for x, y in pareto:
        ax.annotate(f"({x:.1f}, {y:.2f})", (x, y),
                    xytext=(5, 5), textcoords="offset points",
                    fontsize=7, alpha=0.8)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


def pareto_evolution(generations: Optional[dict] = None,
                     save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    gen_data = generations or DEFAULT_GENERATIONS
    n_gens = len(gen_data)
    n_cols = min(3, n_gens)
    n_rows = (n_gens + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_gens == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for ax, (gen_label, points) in zip(axes, gen_data.items()):
        pareto = _compute_pareto_front(points)
        all_x = [p[0] for p in points]
        all_y = [p[1] for p in points]

        ax.scatter(all_x, all_y, c="#4477AA", s=40, alpha=0.4,
                   edgecolors="black", linewidth=0.5)

        if pareto:
            px = [p[0] for p in pareto]
            py = [p[1] for p in pareto]
            ax.plot(px, py, color="#CC3311", linewidth=2,
                    marker="o", markersize=6)

        ax.set_title(gen_label, fontsize=10, fontweight="bold")
        ax.set_xlabel("ATP Yield", fontsize=8)
        ax.set_ylabel("Efficiency", fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(65, 100)
        ax.set_ylim(0.6, 1.05)

    for i in range(n_gens, len(axes)):
        axes[i].axis("off")

    fig.suptitle("Pareto Front Evolution Over Generations", fontsize=13, fontweight="bold")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    r1 = plot_pareto_front(save_path="/tmp/ptr94_pareto_front.png")
    print(f"Pareto front saved to {r1}" if r1 else "Matplotlib not available")
    r2 = pareto_evolution(save_path="/tmp/ptr94_pareto_evolution.png")
    print(f"Pareto evolution saved to {r2}" if r2 else "Matplotlib not available")
