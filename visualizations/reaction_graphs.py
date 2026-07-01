from typing import Optional, List, Tuple, Dict


DEFAULT_NETWORK: Dict[str, List[Tuple[str, str, str]]] = {
    "nodes": [
        ("Glucose", "metabolite"),
        ("G3P", "metabolite"),
        ("Pyruvate", "metabolite"),
        ("Acetyl-CoA", "metabolite"),
        ("CO2", "metabolite"),
        ("NADH", "metabolite"),
        ("FADH2", "metabolite"),
        ("ATP", "metabolite"),
        ("O2", "metabolite"),
        ("H2O", "metabolite"),
        ("Glycolysis", "reaction"),
        ("PDH", "reaction"),
        ("TCA_Cycle", "reaction"),
        ("PCM", "reaction"),
    ],
    "edges": [
        ("Glucose", "Glycolysis"),
        ("Glycolysis", "G3P"),
        ("G3P", "Glycolysis"),
        ("Glycolysis", "Pyruvate"),
        ("Glycolysis", "NADH"),
        ("Glycolysis", "ATP"),
        ("Pyruvate", "PDH"),
        ("PDH", "Acetyl-CoA"),
        ("PDH", "NADH"),
        ("PDH", "CO2"),
        ("Acetyl-CoA", "TCA_Cycle"),
        ("TCA_Cycle", "CO2"),
        ("TCA_Cycle", "NADH"),
        ("TCA_Cycle", "FADH2"),
        ("TCA_Cycle", "ATP"),
        ("NADH", "PCM"),
        ("FADH2", "PCM"),
        ("O2", "PCM"),
        ("PCM", "H2O"),
        ("PCM", "ATP"),
    ],
}

STOICHIOMETRIC_MATRIX_EXAMPLE = {
    "metabolites": ["Glucose", "Pyruvate", "NADH", "FADH2", "ATP", "CO2", "O2", "H2O"],
    "reactions": ["Glycolysis", "PDH", "TCA_Cycle", "PCM"],
    "matrix": [
        [-1,  0,  0,  0],
        [ 2,  0,  0,  0],
        [ 2,  2,  6, -10],
        [ 0,  0,  2,  -2],
        [ 2,  0,  2,  90],
        [ 0,  1,  4,   0],
        [ 0,  0,  0,  -6],
        [ 0,  0,  0,   6],
    ],
}


def reaction_network_graph(network: Optional[dict] = None,
                           save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError:
        return None

    try:
        g = nx.DiGraph()

        data = network or DEFAULT_NETWORK
        for node, ntype in data["nodes"]:
            g.add_node(node, type=ntype)

        for u, v in data["edges"]:
            g.add_edge(u, v)

        pos = nx.spring_layout(g, k=1.5, iterations=50, seed=42)

        fig, ax = plt.subplots(1, 1, figsize=(14, 10))

        metabolite_nodes = [n for n, attr in g.nodes(data=True) if attr.get("type") == "metabolite"]
        reaction_nodes = [n for n, attr in g.nodes(data=True) if attr.get("type") == "reaction"]

        nx.draw_networkx_nodes(
            g, pos, nodelist=metabolite_nodes, node_color="#4477AA",
            node_size=1200, edgecolors="black", linewidths=1.5, ax=ax,
        )
        nx.draw_networkx_nodes(
            g, pos, nodelist=reaction_nodes, node_color="#CCBB44",
            node_size=1000, edgecolors="black", linewidths=1.5,
            node_shape="s", ax=ax,
        )

        nx.draw_networkx_edges(
            g, pos, arrows=True, arrowsize=20,
            edge_color="#666666", width=1.5, alpha=0.6, ax=ax,
        )

        nx.draw_networkx_labels(
            g, pos, font_size=8, font_weight="bold", ax=ax,
        )

        ax.set_title("PTR-94 Reaction Network", fontsize=14, fontweight="bold")
        ax.axis("off")

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            return save_path

        plt.close(fig)
        return None

    except ImportError:
        return None


def stoichiometric_matrix_heatmap(network: Optional[dict] = None,
                                  save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    data = network or STOICHIOMETRIC_MATRIX_EXAMPLE
    mat = np.array(data["matrix"])
    metabolites = data["metabolites"]
    reactions = data["reactions"]

    fig, ax = plt.subplots(1, 1, figsize=(max(8, len(reactions) * 1.5),
                                           max(6, len(metabolites) * 0.8)))

    vmax = max(abs(mat.min()), abs(mat.max()))
    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)

    ax.set_xticks(range(len(reactions)))
    ax.set_xticklabels(reactions, fontsize=9, rotation=45, ha="right")
    ax.set_yticks(range(len(metabolites)))
    ax.set_yticklabels(metabolites, fontsize=9)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            color = "white" if abs(val) > vmax * 0.5 else "black"
            ax.text(j, i, str(val), ha="center", va="center",
                    fontsize=8, fontweight="bold", color=color)

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Stoichiometric Coefficient")
    ax.set_title("PTR-94 Stoichiometric Matrix", fontsize=12, fontweight="bold")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    r = reaction_network_graph(save_path="/tmp/ptr94_reaction_network.png")
    print(f"Reaction network saved to {r}" if r else "networkx/matplotlib not available")
    s = stoichiometric_matrix_heatmap(save_path="/tmp/ptr94_stoich_matrix.png")
    print(f"Stoichiometric matrix saved to {s}" if s else "matplotlib not available")
