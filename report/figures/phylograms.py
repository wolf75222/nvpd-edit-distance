"""Lay out the two checked Newick trees as editable vector graphics."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.phylogeny.tree import build_tree, read_matrix


def tree_layout(tree):
    branches, tips, depths = [], [], []

    def visit(node, depth):
        depths.append(depth)
        if not node.children:
            y = len(tips) * 3.24
            tips.append((node.label, depth, y))
            return y
        heights = []
        for child, length in node.children:
            y = visit(child, depth + length)
            heights.append(y)
            branches.append((depth, y, depth + length, y))
        branches.append((depth, min(heights), depth, max(heights)))
        return sum(heights) / len(heights)

    visit(tree, 0)
    return branches, tips, min(depths), max(depths)


def phylograms(cohort):
    labels, matrix = read_matrix(cohort / "matrix.tsv")
    if len(labels) != 64:
        raise ValueError("expected 64 accessions for the full-page tree figure")
    layouts = {}
    for method in ("upgma", "nj"):
        tree = build_tree(labels, matrix, method)
        if tree.newick() + ";\n" != (cohort / f"{method}.nwk").read_text():
            raise ValueError(f"Newick file does not match the matrix ({method})")
        layouts[method] = tree_layout(tree)

    scale = 108 / sum(high - low for _, _, low, high in layouts.values())
    drawing = [r"\begin{tikzpicture}[x=1mm,y=-1mm,line width=.3pt]"]
    offset = 0
    for method, title in [("upgma", "a) UPGMA"), ("nj", "b) Neighbor-Joining")]:
        branches, tips, low, high = layouts[method]

        def x(depth):
            return offset + scale * (depth - low)

        drawing.append(fr"\node[anchor=west,font=\small\bfseries,inner sep=0pt] at ({offset:.3f},-6) {{{title}}};")
        ink = "powderBlue!50!black" if method == "upgma" else "powderTeal!50!black"
        for x0, y0, x1, y1 in branches:
            drawing.append(fr"\draw[draw={ink}] ({x(x0):.3f},{y0:.3f}) -- ({x(x1):.3f},{y1:.3f});")
        for label, depth, y in tips:
            name = label.replace("_", r"\_")
            drawing.append(fr"\draw[densely dotted,black!36] ({x(depth):.3f},{y:.3f}) -- ({x(high) + 1:.3f},{y:.3f});")
            drawing.append(fr"\node[anchor=west,inner sep=0pt,font=\rmfamily\fontsize{{7}}{{7}}\selectfont] at ({x(high) + 2:.3f},{y:.3f}) {{{name}}};")
        drawing.append(fr"\draw[line width=.6pt] ({offset:.3f},210) -- ++({10*scale:.3f},0) node[midway,below=2pt,font=\scriptsize] {{10 éditions}};")
        offset = x(high) + 28
    drawing.append(r"\end{tikzpicture}")
    return "\n".join(drawing) + "\n"

if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    cohort = directory.parents[1] / "results" / "romeo" / "cohort-837261"
    (directory / "arbres.tex").write_text(phylograms(cohort), encoding="utf-8")
