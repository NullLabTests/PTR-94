from typing import Optional, List, Tuple


DEFAULT_SENSITIVITY = {
    "parameters": [
        ("H⁺ per NADH", 0.85),
        ("PCM Efficiency", 0.65),
        ("H⁺/ATP Ratio", -0.55),
        ("Proton Leak", -0.40),
        ("FADH₂ Contribution", 0.25),
        ("Membrane Integrity", 0.15),
    ],
}

DEFAULT_SOBOL = {
    "parameters": [
        ("H⁺ per NADH", 0.45, 0.52),
        ("PCM Efficiency", 0.28, 0.35),
        ("H⁺/ATP Ratio", 0.18, 0.22),
        ("Proton Leak", 0.12, 0.15),
        ("FADH₂ Contribution", 0.05, 0.08),
        ("Membrane Integrity", 0.02, 0.04),
    ],
}


def tornado_plot(sensitivity_results: Optional[dict] = None,
                 save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    data = sensitivity_results or DEFAULT_SENSITIVITY

    params = data["parameters"]
    names = [p[0] for p in params]
    values = [p[1] for p in params]

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    y_pos = np.arange(len(names))
    colors = ["#CC3311" if v < 0 else "#0077BB" for v in values]

    bars = ax.barh(y_pos, values, color=colors, edgecolor="black", linewidth=0.5, height=0.6)

    for bar, v in zip(bars, values):
        label_x = bar.get_width() + 0.02 if v >= 0 else bar.get_width() - 0.02
        ha = "left" if v >= 0 else "right"
        ax.text(label_x, bar.get_y() + bar.get_height() / 2,
                f"{v:+.3f}", ha=ha, va="center", fontsize=9, fontweight="bold")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=10)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("Sensitivity Coefficient (dimensionless)", fontsize=11)
    ax.set_title("PTR-94 Parameter Sensitivity — Tornado Plot", fontsize=13, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)

    xlim = max(abs(min(values)), abs(max(values))) * 1.3
    ax.set_xlim(-xlim, xlim)

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#0077BB", label="Positive Effect"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#CC3311", label="Negative Effect"),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc="lower right")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


def sobol_plot(sobol_indices: Optional[dict] = None,
               save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    data = sobol_indices or DEFAULT_SOBOL

    params = data["parameters"]
    names = [p[0] for p in params]
    first_order = [p[1] for p in params]
    total_order = [p[2] for p in params]

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    y_pos = np.arange(len(names))
    height = 0.35

    bars1 = ax.barh(y_pos - height / 2, first_order, height,
                    label="First-Order", color="#0077BB",
                    edgecolor="black", linewidth=0.5)
    bars2 = ax.barh(y_pos + height / 2, total_order, height,
                    label="Total-Order", color="#CC3311",
                    edgecolor="black", linewidth=0.5, alpha=0.7)

    for bar, v in zip(bars1, first_order):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:.3f}", ha="left", va="center", fontsize=8, fontweight="bold")
    for bar, v in zip(bars2, total_order):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:.3f}", ha="left", va="center", fontsize=8, fontweight="bold")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel("Sobol Index", fontsize=11)
    ax.set_title("PTR-94 Global Sensitivity — Sobol Indices", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    ax.set_xlim(0, 1)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    r1 = tornado_plot(save_path="/tmp/ptr94_tornado.png")
    print(f"Tornado plot saved to {r1}" if r1 else "Matplotlib not available")
    r2 = sobol_plot(save_path="/tmp/ptr94_sobol.png")
    print(f"Sobol plot saved to {r2}" if r2 else "Matplotlib not available")
