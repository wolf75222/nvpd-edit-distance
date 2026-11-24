"""Check the structure and cross-references of both LaTeX reports."""
import json
from pathlib import Path
import re
import unittest
from urllib.parse import unquote


REPORT = Path(__file__).resolve().parents[1] / "report"
REPORTS = ("rapport.tex", "rapport-court.tex")
FOUNDATIONAL_DOIS = {
    "sadiq": "10.1007/s10586-019-02962-w",
    "wfagpu": "10.1109/access.2022.3182714",
    "hirschberg": "10.1145/360825.360861",
    "wfa": "10.1093/bioinformatics/btaa777",
    "biwfa": "10.1093/bioinformatics/btad074",
    "edlib": "10.1093/bioinformatics/btw753",
    "quicked": "10.1093/bioinformatics/btaf112",
}


def expand_sections(path):
    source = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))
    return re.sub(r"\\input\{(sections/[^}]+)\}",
                  lambda match: expand_sections(REPORT / (match[1] + ".tex")),
                  source)


class ReportStructure(unittest.TestCase):
    def test_figures_and_algorithms_have_text_references(self):
        for report in REPORTS:
            text = expand_sections(REPORT / report)
            label = r"\{((?:fig|alg):[\w-]+)\}"
            references = set(re.findall(r"\\ref" + label, text))
            labels = set(re.findall(label, re.sub(r"\\ref" + label, "", text)))
            with self.subTest(report=report):
                self.assertTrue(labels)
                self.assertEqual(labels, references)

    def test_one_summary_per_closing_passage(self):
        for report in REPORTS:
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

    def test_bibliography_identifiers_and_citations(self):
        source = (REPORT / "tex" / "references.bib").read_text(encoding="utf-8")
        entries = list(re.finditer(r"^\s*@(?!comment\b|preamble\b|string\b)\w+\s*\{\s*([^\s,]+)\s*,",
                                   source, re.I | re.M))
        keys = [entry[1] for entry in entries]
        folded_keys = [key.casefold() for key in keys]
        self.assertCountEqual(folded_keys, set(folded_keys), "duplicate bibliography keys")
        dois = {}
        for index, entry in enumerate(entries):
            end = entries[index + 1].start() if index + 1 < len(entries) else len(source)
            fields = source[entry.end():end]
            doi_fields = re.findall(r'(?:^|,)\s*doi\s*=\s*(\{[^}]*\}|"[^"]*"|[^,}\n]*)', fields, re.I)
            with self.subTest(key=entry[1]):
                if not doi_fields:
                    self.assertRegex(fields, r'\burl\s*=\s*[{"]https?://', "missing DOI or source URL")
            for value in doi_fields:
                doi = value.strip('{}" ').lower()
                with self.subTest(key=entry[1], doi=doi):
                    self.assertRegex(doi, r"^10\.\d{4,9}/[^\s<>{}]+$")
                    self.assertNotIn(entry[1], dois, "duplicate DOI field")
                    dois[entry[1]] = doi
        self.assertCountEqual(dois.values(), set(dois.values()), "duplicate DOI records")
        for key, doi in FOUNDATIONAL_DOIS.items():
            with self.subTest(key=key):
                self.assertEqual(dois.get(key), doi)

        notebook = REPORT.parent / "notebooks" / "profil_nvpd.ipynb"
        cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
        markdown = "\n".join("".join(cell["source"]) for cell in cells if cell["cell_type"] == "markdown")
        doi_link = r"https://doi\.org/((?:[^()\s<>]|\([^()\s]*\))+)"
        notebook_dois = {unquote(doi).lower() for doi in re.findall(doi_link, markdown)}
        self.assertFalse(notebook_dois - set(dois.values()), "notebook DOI missing from the catalogue")
        self.assertFalse(set(FOUNDATIONAL_DOIS.values()) - notebook_dois, "foundational DOI missing from notebook")

        cited = set()
        for report in REPORTS:
            groups = re.findall(r"\\cite(?:\[[^]]*\]){0,2}\s*\{([^{}]+)\}", expand_sections(REPORT / report))
            report_keys = {key.strip() for group in groups for key in group.split(",")}
            with self.subTest(report=report):
                self.assertFalse(report_keys - set(keys), f"undefined citations: {report_keys - set(keys)}")
                missing = FOUNDATIONAL_DOIS.keys() - report_keys
                self.assertFalse(missing, f"uncited foundational papers: {missing}")
            cited.update(report_keys)
        self.assertEqual(set(keys) - cited, set(), "uncited bibliography records")
