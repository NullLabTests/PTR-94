from typing import Optional


DEFAULT_ENERGY_BALANCE = {
    "glucose_delta_g": 2870.0,
    "glycolysis_atp": 2 * 30.5,
    "glycolysis_nadh": 2 * 7.5 * 30.5,
    "pdh_tca_atp": 2 * 30.5,
    "pdh_tca_nadh": 8 * 7.5 * 30.5,
    "pdh_tca_fadh2": 2 * 7.5 * 30.5,
    "pcm_redox_atp": 90 * 30.5,
    "total_atp": 94 * 30.5,
    "dissipated": 2870.0 - 94 * 30.5,
}


def energy_sankey(energy_balance_data: Optional[dict] = None,
                  save_path: Optional[str] = None) -> Optional[str]:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
    except ImportError:
        return None

    data = energy_balance_data or DEFAULT_ENERGY_BALANCE

    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("PTR-94 Energy Flow (kJ/mol)", fontsize=14, fontweight="bold")

    total_in = data["glucose_delta_g"]
    glyco_energy = data["glycolysis_atp"] + data["glycolysis_nadh"]
    redox_in = data["pdh_tca_nadh"] + data["pdh_tca_fadh2"]
    atp_out = data["total_atp"]
    diss = data["dissipated"]

    module_colors = {
        "glucose": "#4477AA",
        "glycolysis": "#66CCEE",
        "pdh_tca": "#228833",
        "pcm": "#CCBB44",
        "atp": "#AA3377",
        "dissipation": "#BBBBBB",
    }

    left_positions = {
        "glucose": (0.5, 3.0),
        "glycolysis": (2.5, 3.5),
        "pdh_tca": (4.5, 3.5),
        "pcm": (6.5, 3.0),
        "atp": (8.5, 3.5),
        "dissipation": (8.5, 2.0),
    }

    module_labels = {
        "glucose": f"Glucose\nΔG°' = {total_in:.0f}",
        "glycolysis": f"Glycolysis\n{glyco_energy:.0f} kJ/mol\n(2 ATP + 2 NADH)",
        "pdh_tca": f"PDH + TCA\n{redox_in:.0f} kJ/mol\n(2 ATP + 8 NADH + 2 FADH₂)",
        "pcm": f"PCM\n{90 * 30.5:.0f} kJ/mol\n(90 ATP from redox)",
        "atp": f"ATP\n{atp_out:.0f} kJ/mol\n(94 ATP)",
        "dissipation": f"Dissipated\n{diss:.0f} kJ/mol",
    }

    box_w = 1.5
    box_h = 0.8
    for key, (cx, cy) in left_positions.items():
        color = module_colors.get(key, "#AAAAAA")
        rect = mpatches.FancyBboxPatch(
            (cx - box_w / 2, cy - box_h / 2),
            box_w, box_h,
            boxstyle="round,pad=0.1",
            facecolor=color, edgecolor="black", alpha=0.8,
        )
        ax.add_patch(rect)
        ax.text(cx, cy, module_labels.get(key, key),
                ha="center", va="center", fontsize=8, fontweight="bold")

    def draw_flow(x1, y1, x2, y2, width, label="", color="#888888"):
        ax.annotate(
            "",
            xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="->",
                color=color,
                lw=width * 3,
                connectionstyle="arc3,rad=0.1",
                alpha=0.6,
            ),
        )
        if label:
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 + 0.3
            ax.text(mx, my, label, ha="center", va="bottom",
                    fontsize=7, fontweight="bold", color=color)

    draw_flow(
        2.0, 3.0, 3.25, 3.5, width=0.4,
        label=f"{glyco_energy:.0f} kJ",
        color=module_colors["glycolysis"],
    )
    draw_flow(
        4.0, 3.5, 5.25, 3.0, width=0.4,
        label=f"{redox_in:.0f} kJ",
        color=module_colors["pdh_tca"],
    )
    draw_flow(
        7.25, 3.0, 8.0, 3.5, width=0.4,
        label=f"{90 * 30.5:.0f} kJ",
        color=module_colors["pcm"],
    )
    draw_flow(
        2.0, 3.0, 3.25, 2.5, width=0.1,
        label=f"",
        color=module_colors["glucose"],
    )

    legend_elements = [
        mpatches.Patch(color=c, label=l, alpha=0.8)
        for l, c in [
            ("Glucose Energy", module_colors["glucose"]),
            ("Glycolysis", module_colors["glycolysis"]),
            ("PDH + TCA", module_colors["pdh_tca"]),
            ("PCM (Perfect Coupling)", module_colors["pcm"]),
            ("ATP Output", module_colors["atp"]),
            ("Dissipation", module_colors["dissipation"]),
        ]
    ]
    ax.legend(handles=legend_elements, loc="lower center",
              ncol=3, fontsize=8, framealpha=0.9)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return save_path

    plt.close(fig)
    return None


if __name__ == "__main__":
    result = energy_sankey(save_path="/tmp/ptr94_sankey.png")
    print(f"Sankey diagram saved to {result}" if result else "Matplotlib not available")
