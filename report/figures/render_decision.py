"""Export the report’s TikZ diagrams used in the notebook as SVG files."""

from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


def main() -> None:
    figures = Path(__file__).resolve().parent
    report = figures.parent
    with TemporaryDirectory() as directory:
        work = Path(directory)
        (work / "tex").symlink_to(report / "tex", target_is_directory=True)
        (work / "figures").symlink_to(figures, target_is_directory=True)
        template = r"""\documentclass[10pt,a4paper]{article}
\input{tex/preambule}
\usepackage[active,tightpage]{preview}
\setlength{\PreviewBorder}{3pt}
\begin{document}
\begin{preview}\input{figures/DIAGRAM}\end{preview}
\end{document}
"""
        for name in ("decision", "contexte-nvpd", "contexte-ewfa"):
            (work / "diagram.tex").write_text(template.replace("DIAGRAM", name), encoding="utf-8")
            subprocess.run(["tectonic", "diagram.tex"], cwd=work, check=True)
            subprocess.run(["pdftocairo", "-svg", "diagram.pdf", "diagram.svg"], cwd=work, check=True)
            (figures / f"{name}.svg").write_bytes((work / "diagram.svg").read_bytes())


if __name__ == "__main__":
    main()
