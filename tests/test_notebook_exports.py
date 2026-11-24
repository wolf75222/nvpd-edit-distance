"""Check that the executed notebook and the published plots stay in sync."""
import ast
import base64
from hashlib import sha256
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class NotebookExports(unittest.TestCase):
    def test_embedded_png_and_caption(self):
        cells = json.loads((ROOT / "notebooks/profil_nvpd.ipynb").read_text())["cells"]
        self.assertEqual(len({cell["id"] for cell in cells}), len(cells))
        published = set()
        vectors = set()
        for index, cell in enumerate(cells):
            if cell["cell_type"] != "code":
                continue
            self.assertIsNotNone(cell["execution_count"], cell["id"])
            outputs = cell["outputs"]
            self.assertFalse(any(o["output_type"] == "error" for o in outputs), cell["id"])
            images = [o["data"]["image/png"] for o in outputs if "image/png" in o.get("data", {})]
            if not images:
                continue
            with self.subTest(cell=cell["id"]):
                code = ast.parse("".join(cell["source"]))
                calls = [x for x in ast.walk(code) if isinstance(x, ast.Call)
                         and isinstance(x.func, ast.Name) and x.func.id == "nv_save"]
                self.assertEqual(len(calls), 1)
                stem = ast.literal_eval(calls[0].args[1])
                path = ROOT / "figures" / (stem + ".png")
                if any(k.arg == "vector" and ast.literal_eval(k.value) for k in calls[0].keywords):
                    pdf = path.with_suffix(".pdf")
                    self.assertTrue(pdf.read_bytes().startswith(b"%PDF-"), str(pdf))
                    vectors.add(pdf)
                self.assertNotIn(path, published)
                published.add(path)
                self.assertEqual(len(images), 1)
                self.assertEqual(sha256(base64.b64decode("".join(images[0]))).digest(),
                                 sha256(path.read_bytes()).digest(), str(path))
                caption = "".join(cells[index + 1]["source"])
                self.assertTrue(caption.startswith(f"*Figure {len(published)}.*"), caption[:60])
        exports = {p for p in (ROOT / "figures").rglob("*.png") if p.name != "montage.png"}
        self.assertTrue(published)
        self.assertEqual(published, exports)
        self.assertEqual(vectors, set((ROOT / "figures").rglob("*.pdf")))


if __name__ == "__main__":
    unittest.main()
