from typing import Optional


DEFAULT_ATP_YIELD = 94
DEFAULT_EFFICIENCY = 94 * 30.5 / 2870.0


def energy_flow_chart(atp_yield: int = DEFAULT_ATP_YIELD,
                      efficiency: float = DEFAULT_EFFICIENCY,
                      save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return None

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    module_labels = ["Glycolysis\n(2 ATP + 2 NADH)",
                     "PDH + TCA\n(2 ATP + 8 NADH + 2 FADH₂)",
                     "PCM\n(Redox → 90 ATP)"]
    substrate_atp = [2, 2, 0]
    redox_atp = [15, 75, 0]
    total_bars = [17, 77, 0]
    all_groups = [substrate_atp, redox_atp]
    group_names = ["Substrate-level ATP", "Redox-derived ATP"]
    colors = ["#4477AA", "#CCBB44"]

    ax = axes[0]
    x = np.arange(len(module_labels))
    bottom = np.zeros(len(module_labels))
    for i, (vals, name) in enumerate(zip(all_groups, group_names)):
        bars = ax.bar(x, vals, bottom=bottom, label=name, color=colors[i], edgecolor="black", linewidth=0.5)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_y() + v / 2,
                        str(v), ha="center", va="center", fontsize=9, fontweight="bold")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(module_labels, fontsize=8)
    ax.set_ylabel("ATP Count", fontsize=10)
    ax.set_title("Module-by-Module ATP Contribution", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 100)

    ax2 = axes[1]
    mod_names = ["Glycolysis", "PDH + TCA", "PCM", "Total"]
    mod_vals = [17, 77, 0, atp_yield]
    mod_colors = ["#66CCEE", "#228833", "#CCBB44", "#AA3377"]
    bars2 = ax2.barh(mod_names, mod_vals, color=mod_colors, edgecolor="black", linewidth=0.5)
    for bar, v in zip(bars2, mod_vals):
        if v > 0:
            ax2.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                     str(v), ha="left", va="center", fontsize=9, fontweight="bold")
    ax2.set_xlabel("ATP per Glucose", fontsize=10)
    ax2.set_title("Horizontal ATP Contribution", fontsize=11, fontweight="bold")
    ax2.set_xlim(0, 100)

    ax3 = axes[2]
    categories = ["Eukaryote", "Prokaryote", "PTR-94", "Ceiling"]
    values = [31, 37, atp_yield, 94]
    bar_colors = ["#EE6677", "#228833", "#AA3377", "#888888"]
    bars3 = ax3.bar(categories, values, color=bar_colors, edgecolor="black", linewidth=0.5)
    for bar, v in zip(bars3, values):
        ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 str(v), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax3.set_ylabel("ATP per Glucose", fontsize=10)
    ax3.set_title("Comparison: ATP Yield", fontsize=11, fontweight="bold")
    ax3.set_ylim(0, 105)

    fig.suptitle(f"PTR-94 Energy Flow — {atp_yield} ATP per Glucose (Efficiency: {efficiency * 100:.1f}%)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    result = energy_flow_chart(save_path="/tmp/ptr94_energy_flow.png")
    print(f"Energy flow chart saved to {result}" if result else "Matplotlib not available")
