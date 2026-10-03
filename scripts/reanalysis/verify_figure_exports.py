"""Optional export QA; requires Pillow, pypdf, and Poppler's pdftoppm."""

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path

from PIL import Image
from pypdf import PdfReader


def nonwhite_fraction(path):
    with Image.open(path) as image:
        gray = image.convert("L")
        histogram = gray.histogram()
        fraction = sum(histogram[:245]) / (gray.width * gray.height)
        return gray.width, gray.height, fraction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poppler", default="pdftoppm")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    figures = root / "reanalysis" / "figures"
    scratch = root / "tmp" / "figure_qa"
    rows = []
    pdfs = sorted(figures.rglob("*.pdf"))
    pngs = sorted(figures.rglob("*.png"))
    assert pdfs and {p.with_suffix("") for p in pdfs} == {p.with_suffix("") for p in pngs}
    for pdf in pdfs:
        pages = len(PdfReader(pdf).pages)
        assert pages == 1, f"Expected a one-page figure: {pdf}"
        relative = pdf.relative_to(figures)
        prefix = (scratch / relative).with_suffix("")
        prefix.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [args.poppler, "-singlefile", "-r", "120", "-png", str(pdf), str(prefix)],
            check=True,
            capture_output=True,
        )
        pdf_width, pdf_height, pdf_ink = nonwhite_fraction(prefix.with_suffix(".png"))
        png_width, png_height, png_ink = nonwhite_fraction(pdf.with_suffix(".png"))
        assert pdf_ink > 0.01 and png_ink > 0.01, f"Blank or nearly blank export: {pdf}"
        assert abs(pdf_width / pdf_height - png_width / png_height) < 0.02
        rows.append(
            {
                "Figure": relative.as_posix(),
                "PDF_Pages": pages,
                "PDF_MD5": hashlib.md5(pdf.read_bytes()).hexdigest(),
                "PNG_Width": png_width,
                "PNG_Height": png_height,
                "PNG_Nonwhite_Fraction": png_ink,
                "PDF_Render_Nonwhite_Fraction": pdf_ink,
                "Checks": "PASS",
            }
        )
    output = root / "reanalysis" / "provenance" / "figure_checks.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Export checks passed for {len(rows)} PNG/PDF pairs.")
    print("Renders still require visual review for clipping, labels, and scientific content.")


if __name__ == "__main__":
    main()
