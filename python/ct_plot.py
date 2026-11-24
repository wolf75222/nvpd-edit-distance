"""Charte des figures du carnet CHPS1001.

La palette poudre distingue les séries. Les marqueurs gardent les courbes
lisibles en niveaux de gris. Les figures utilisent un texte serif noir, des
cadres fins, une grille pointillée et des graduations majeures rentrantes.
Les en-têtes de colonne nomment les comparaisons au-dessus des axes.
"""

import re
from contextlib import contextmanager
from uuid import uuid4

import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator, ScalarFormatter


COLORS = ["#8b9dc3", "#c98b9a", "#7fa8a1", "#d9c28a", "#a5a5a5", "#a88fb8"]
INK = "#000000"
PAPER = "#ffffff"
GUIDE = dict(color=INK, lw=0.7, ls="--", zorder=1)


class _FrenchScalarFormatter(ScalarFormatter):
    """Garder la notation scientifique de Matplotlib et sa virgule sans espacement."""

    @staticmethod
    def _comma(text):
        return re.sub(r"(?<=\d)\.(?=\d)", "{,}" if "$" in text else ",", text)

    def __call__(self, value, pos=None):
        return self._comma(super().__call__(value, pos))

    def get_offset(self):
        return self._comma(super().get_offset())


def use_style():
    """Appliquer la charte avant de construire les axes."""
    plt.rcParams.update({
        "figure.dpi": 150, "figure.constrained_layout.use": True, "savefig.bbox": "tight",
        "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
        "font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif",
                                               "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.prop_cycle": plt.cycler(color=COLORS), "axes.linewidth": 0.8, "axes.edgecolor": INK,
        "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.axisbelow": True, "axes.grid": True, "axes.grid.which": "major",
        "grid.linestyle": ":", "grid.color": INK, "grid.linewidth": 0.6,
        "lines.linewidth": 1.3, "lines.markersize": 4, "lines.markeredgewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": False, "ytick.minor.visible": False,
        "legend.fancybox": False, "legend.edgecolor": INK, "legend.framealpha": 1,
    })


def label_panels(axes):
    """Marquer les panneaux et supprimer aussi les mineures des axes logarithmiques."""
    for letter, axis in zip("abcdefghijklmnopqrstuvwxyz", axes):
        for dimension in (axis.xaxis, axis.yaxis):
            dimension.set_minor_locator(NullLocator())
            original = dimension.get_major_formatter()
            if type(original) is ScalarFormatter:
                dimension.set_major_formatter(_FrenchScalarFormatter(
                    useOffset=original.get_useOffset(), useMathText=original.get_useMathText()))
        axis.annotate(f"{letter})", (0, 1), xycoords="axes fraction", xytext=(0, 4),
                      textcoords="offset points", va="bottom", fontweight="bold")


def column_titles(axes, names):
    """Placer les en-têtes de colonne, distincts des légendes et des titres d'axe."""
    if len(axes) != len(names):
        raise ValueError("Chaque colonne doit avoir un en-tête.")
    for axis, name in zip(axes, names):
        axis.annotate(name, (.5, 1), xycoords="axes fraction", xytext=(0, 4),
                      textcoords="offset points", ha="center", va="bottom", fontweight="bold")


def legend_above(figure, axes, ncol):
    """Rassembler les catégories communes au-dessus de la figure."""
    entries = {}
    for axis in axes:
        handles, labels = axis.get_legend_handles_labels()
        entries.update(zip(labels, handles))
    figure.legend(entries.values(), entries.keys(), loc="outside upper center",
                  ncol=ncol, columnspacing=1.25, handlelength=1.7, fontsize=7)


@contextmanager
def atomic_png(destination):
    """Remplacer le PNG seulement après la fermeture d'un fichier temporaire voisin."""
    temporary = destination.with_name(f".{destination.name}.{uuid4().hex}.tmp")
    stream = temporary.open("xb")
    try:
        with stream:
            yield stream
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def montage(directory, entries):
    """Réunir les exports du carnet dans l'ordre de lecture."""
    from matplotlib.font_manager import FontProperties, findfont
    from PIL import Image, ImageDraw, ImageFont

    width, height, columns = 760, 375, 3
    nrows = (len(entries) + columns - 1) // columns
    sheet = Image.new("RGB", (width * columns, height * nrows), PAPER)
    pen = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(findfont(FontProperties(family="serif")), 17)
    for index, (filename, heading) in enumerate(entries):
        x, y = index % columns * width, index // columns * height
        with Image.open(directory / filename) as source:
            panel = source.convert("RGB")
            panel.thumbnail((width - 24, height - 38), Image.Resampling.LANCZOS)
        sheet.paste(panel, (x + (width - panel.width) // 2, y))
        pen.text((x + 12, y + height - 30), f"Figure {index+1}. {heading}", fill=INK, font=font)
    with atomic_png(directory / "montage.png") as stream:
        sheet.save(stream, format="PNG", dpi=(150, 150))
