"""Mise en forme des figures du profil NvPD."""
import matplotlib.pyplot as plt

COLORS = ["#8b9dc3", "#c98b9a", "#7fa8a1", "#d9c28a", "#a5a5a5", "#a88fb8"]
GUIDE = dict(color="0.35", lw=0.7, ls="--", zorder=1)


def use_style():
    plt.rcParams.update({
        "figure.dpi": 150, "figure.constrained_layout.use": True, "savefig.bbox": "tight",
        "font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif",
                                                "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.prop_cycle": plt.cycler(color=COLORS), "axes.linewidth": 0.8,
        "axes.axisbelow": True, "axes.grid": True, "grid.linestyle": ":",
        "grid.color": "black", "grid.linewidth": 0.6,
        "lines.linewidth": 1.3, "lines.markersize": 4, "lines.markeredgewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "legend.fancybox": False, "legend.edgecolor": "black", "legend.framealpha": 1.0,
    })


def label_panels(axes):
    for letter, axis in zip("abcdefghijklmnopqrstuvwxyz", axes):
        axis.annotate(f"({letter})", (0, 1), xycoords="axes fraction", xytext=(0, 4),
                      textcoords="offset points", va="bottom")


def legend_above(figure, axes, ncol):
    lines = {}
    for axis in axes:
        handles, labels = axis.get_legend_handles_labels()
        lines.update(zip(labels, handles))
    figure.legend(lines.values(), lines.keys(), loc="outside upper center",
                  ncol=ncol, columnspacing=1.3, handlelength=1.8, fontsize=7)
