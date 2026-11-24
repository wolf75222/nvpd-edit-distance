"""Lay out the two checked Newick trees as editable vector graphics."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.tree import build_tree, read_matrix


def phylograms(cohort):
    labels, matrix = read_matrix(cohort / "matrix.tsv")
    if len(labels) != 64:
        raise ValueError("expected 64 accessions for the full-page tree figure")
    drawing = [r"\begin{tikzpicture}[x=1mm,y=-1mm,line width=.3pt]"]
    for offset, method, title in [(0, "upgma", "UPGMA"), (82, "nj", "Neighbor-Joining")]:
        tree = build_tree(labels, matrix, method)
        if tree.newick() + ";\n" != (cohort / f"{method}.nwk").read_text():
            raise ValueError(f"Newick file does not match the matrix ({method})")
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
        scale = 54 / (max(depths) - min(depths))

        def x(depth):
            return offset + scale * (depth - min(depths))

        drawing.append(fr"\node[anchor=west,font=\small\bfseries,inner sep=0pt] at ({offset},-6) {{{title}}};")
        for x0, y0, x1, y1 in branches:
            drawing.append(fr"\draw ({x(x0):.3f},{y0:.3f}) -- ({x(x1):.3f},{y1:.3f});")
        for label, depth, y in tips:
            name = label.replace("_", r"\_")
            drawing.append(fr"\draw[densely dotted,black!36] ({x(depth):.3f},{y:.3f}) -- ({offset + 55},{y:.3f});")
            drawing.append(fr"\node[anchor=west,inner sep=0pt,font=\sffamily\fontsize{{7}}{{7}}\selectfont] at ({offset + 56},{y:.3f}) {{{name}}};")
        drawing.append(fr"\draw[line width=.6pt] ({offset},210) -- ++({10*scale:.3f},0) node[midway,below=2pt,font=\scriptsize] {{10 éditions}};")
    drawing.append(r"\end{tikzpicture}")
    return "\n".join(drawing) + "\n"
