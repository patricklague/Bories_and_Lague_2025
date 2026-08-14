from pathlib import Path

import matplotlib.pyplot as plt


N_ROWS = 8
N_COLS = 2
COLUMN_TITLES = ["0.2M", "0.1M"]
ROW_LABELS = ["VAL", "LEU", "ILE", "MET", "PHE", "TYR", "TRP", r"ARG$^0$"]
REGION_BOUNDS = [0.0, 9.5, 19.5, 30.0]


fig, axes = plt.subplots(
    N_ROWS,
    N_COLS,
    figsize=(8, 1.5 * N_ROWS),
    sharex=True,
    dpi=600,
)

vertical_guides = sorted({-bound for bound in REGION_BOUNDS} | set(REGION_BOUNDS))

for row_index, row_axes in enumerate(axes):
    for column_index, ax in enumerate(row_axes):
        for x_value in vertical_guides:
            ax.axvline(
                x=x_value,
                linestyle="--",
                alpha=0.8,
                linewidth=0.8,
                color="gray",
                zorder=0,
            )

        ax.set_xlim(-40, 40)
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="y", linestyle="--", alpha=0.4, linewidth=0.3)

        if row_index == 0:
            ax.set_title(
                COLUMN_TITLES[column_index],
                fontsize=11,
                fontweight="bold",
            )

    row_axes[0].set_ylabel(
        ROW_LABELS[row_index],
        rotation=0,
        labelpad=30,
        va="center",
        fontsize=10,
        fontweight="bold",
    )

for ax in axes[-1, :]:
    ax.set_xlabel(r"z ($\AA$)")

fig.subplots_adjust(left=0.15, hspace=0.15, wspace=0.35)

output_dir = Path(__file__).resolve().parent.parent / "plot" / "FigureSS1"
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "FigureSS1_0p2M_0p1M.png"

fig.savefig(output_path, dpi=600, bbox_inches="tight", transparent=True)
plt.show()
print(f"Saved -> {output_path}")