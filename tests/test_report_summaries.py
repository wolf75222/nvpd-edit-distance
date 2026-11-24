"""Check the structure and cross-references of both LaTeX reports."""
from pathlib import Path
import re
import unittest


REPORT = Path(__file__).resolve().parents[1] / "report"


def expand_sections(path):
    return re.sub(r"\\input\{(sections/[^}]+)\}",
                  lambda match: expand_sections(REPORT / (match[1] + ".tex")),
                  path.read_text(encoding="utf-8"))


class ReportStructure(unittest.TestCase):
    def test_figures_and_algorithms_have_text_references(self):
        for report in ("rapport.tex", "rapport-court.tex"):
            text = expand_sections(REPORT / report)
            label = r"\{((?:fig|alg):[\w-]+)\}"
            references = set(re.findall(r"\\ref" + label, text))
            labels = set(re.findall(label, re.sub(r"\\ref" + label, "", text)))
            with self.subTest(report=report):
                self.assertTrue(labels)
                self.assertEqual(labels, references)

    def test_one_summary_per_closing_passage(self):
        for report in ("rapport.tex", "rapport-court.tex"):
            source = expand_sections(REPORT / report)
            headings = list(re.finditer(r"^\\(section|subsection)\{([^}]+)\}", source, re.M))
            self.assertTrue(headings, report)
            for index, heading in enumerate(headings):
                following = headings[index + 1] if index + 1 < len(headings) else None
                passage = source[heading.end():following.start() if following else len(source)]
                introduces_subsections = (heading[1] == "section" and following is not None
                                          and following[1] == "subsection")
                with self.subTest(report=report, section=heading[2]):
                    summaries = re.findall(r"^\s*\\synthese(?:section)?\{", passage, re.M)
                    self.assertEqual(len(summaries), 0 if introduces_subsections else 1)
