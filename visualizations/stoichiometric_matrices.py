from typing import Optional, List


DEFAULT_MATRIX = [
    [-1,  0,  0,  0,  0],
    [ 2,  0, -1,  0,  0],
    [ 2,  2,  6,  0, -10],
    [ 0,  0,  2,  0,  -2],
    [ 2,  0,  2,  0,  90],
    [ 0,  1,  4,  0,   0],
    [ 0,  0,  0, -6,   0],
    [ 0,  0,  0,  6,   0],
    [-1,  0,  0,  0,   0],
    [ 0, -1,  0,  0,   0],
]

DEFAULT_METABOLITES = [
    "Glucose", "Pyruvate", "NADH", "FADH₂",
    "ATP", "CO₂", "O₂", "H₂O", "ADP", "Pᵢ",
]

DEFAULT_REACTIONS = [
    "Glycolysis", "PDH", "TCA", "PCM (NADH)", "PCM (FADH₂)",
]

DEFAULT_NULLSPACE = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]

DEFAULT_MOIETY_NAMES = ["Carbon Moiety", "Phosphate Moiety", "Electron Moiety"]


def stoichiometric_heatmap(matrix: Optional[List[List[float]]] = None,
                           metabolite_names: Optional[List[str]] = None,
                           reaction_names: Optional[List[str]] = None,
                           save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    mat = matrix or DEFAULT_MATRIX
    metabolites = metabolite_names or DEFAULT_METABOLITES
    reactions = reaction_names or DEFAULT_REACTIONS

    mat_array = np.array(mat, dtype=float)

    fig, ax = plt.subplots(1, 1, figsize=(max(8, len(reactions) * 1.8),
                                           max(6, len(metabolites) * 0.7)))

    vmax = max(abs(mat_array.min()), abs(mat_array.max()))
    vmax = max(vmax, 1e-10)

    im = ax.imshow(mat_array, cmap="RdBu_r", aspect="auto",
                   vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(len(reactions)))
    ax.set_xticklabels(reactions, fontsize=9, rotation=45, ha="right")
    ax.set_yticks(range(len(metabolites)))
    ax.set_yticklabels(metabolites, fontsize=9)

    for i in range(mat_array.shape[0]):
        for j in range(mat_array.shape[1]):
            val = mat_array[i, j]
            if abs(val) > 0:
                color = "white" if abs(val) > vmax * 0.5 else "black"
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=7, fontweight="bold", color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Stoichiometric Coefficient", fontsize=10)

    ax.set_title("PTR-94 Stoichiometric Matrix\n(metabolites × reactions)",
                 fontsize=12, fontweight="bold")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


def conserved_moiety_plot(nullspace: Optional[List[List[float]]] = None,
                          moiety_names: Optional[List[str]] = None,
                          save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    ns = nullspace or DEFAULT_NULLSPACE
    names = moiety_names or DEFAULT_MOIETY_NAMES

    ns_array = np.array(ns, dtype=float)

    fig, axes = plt.subplots(1, ns_array.shape[1], figsize=(5 * ns_array.shape[1], 5))
    if ns_array.shape[1] == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        col = ns_array[:, i]
        max_abs = max(abs(col).max(), 1e-10)

        ax.bar(range(len(col)), col, color="#4477AA", edgecolor="black",
               linewidth=0.5, width=0.6)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_ylim(-max_abs * 1.3, max_abs * 1.3)
        ax.set_xlabel("Metabolite Index", fontsize=9)
        ax.set_ylabel("Conserved Quantity", fontsize=9)
        ax.set_title(f"Moiety {i + 1}: {names[i] if i < len(names) else 'Unknown'}",
                     fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.3)

        for j, v in enumerate(col):
            ax.text(j, v + (0.05 if v >= 0 else -0.12),
                    f"{v:.2f}", ha="center", va="bottom" if v >= 0 else "top",
                    fontsize=7, rotation=45)

    fig.suptitle("PTR-94 Conserved Moieties (Left Nullspace)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    r1 = stoichiometric_heatmap(save_path="/tmp/ptr94_stoich_heatmap.png")
    print(f"Stoichiometric heatmap saved to {r1}" if r1 else "Matplotlib not available")
    r2 = conserved_moiety_plot(save_path="/tmp/ptr94_moiety.png")
    print(f"Conserved moiety plot saved to {r2}" if r2 else "Matplotlib not available")
