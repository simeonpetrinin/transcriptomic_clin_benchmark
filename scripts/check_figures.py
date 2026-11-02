#Redraws every figure of the paper from the archived results and checks it against the published figure
#Each figure script is run into a temporary folder. A redrawn figure is "identical" when its SHA-256 matches the published file
#(listed in figure_generation/paper_figures.sha256); with --reference-dir (a folder holding the published PNGs) a figure that is not
#byte-identical is also compared pixel by pixel, which separates real differences from encoder or font differences
#Byte identity needs the package versions in requirements.txt and the Inter font; Figure 1 is a drawing, not made from the results
#Usage, from the project root: python scripts/check_figures.py --results-dir results [--reference-dir paper_figures] [--keep-dir figures_check]
import argparse, hashlib, os, subprocess, sys, tempfile, time
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FIGS = [("Figure 2", "fig2.py", ["fig2.png"]), ("Figure 3", "fig3.py", ["fig3.png"]), ("Figure 4", "fig4.py", ["fig4.png"]),
        ("Figure 5", "fig5.py", ["fig5.png"]), ("Figure 6", "fig6.py", ["fig6.png"]), ("Figure 7", "fig7.py", ["fig7.png"]),
        ("Figure 8", "fig8.py", ["fig8.png"]), ("Figure S1", "gene_sweep.py", ["gene_sweep_boxplot.png"]),
        ("Figure S2", "runtime.py", ["runtime_figure.png"]), ("Figures S3-S4", "external_ibs_dcal.py", ["external_ibs.png", "external_dcal.png"])]
TOL = 16  #a pixel differs when any channel differs by more than this (out of 255)

def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

#share of pixels that differ, compared in strips of rows to keep memory low; None when the sizes differ
def pixel_diff(a, b):
    A, B = Image.open(a).convert("RGB"), Image.open(b).convert("RGB")
    if A.size != B.size: return None
    w, h = A.size; n = 0
    for y in range(0, h, 512):
        box = (0, y, w, min(h, y + 512))
        d = np.abs(np.asarray(A.crop(box), np.int16) - np.asarray(B.crop(box), np.int16)).max(axis=2)
        n += int((d > TOL).sum())
    return n / (w * h)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--reference-dir", default=None)
    ap.add_argument("--keep-dir", default=None, help="keep the redrawn figures here instead of a temporary folder")
    a = ap.parse_args()
    expected = {}
    with open(os.path.join(ROOT, "figure_generation", "paper_figures.sha256")) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"): h, name = line.split(); expected[name] = h
    out = a.keep_dir or tempfile.mkdtemp(prefix="figures_check_"); os.makedirs(out, exist_ok=True)
    rows, bad = [], 0
    for label, script, files in FIGS:
        t0 = time.time()
        p = subprocess.run([sys.executable, os.path.join(ROOT, "figure_generation", script), "--results-dir", a.results_dir, "--out-dir", out],
                           capture_output=True, text=True)
        secs = time.time() - t0
        for f in files:
            new = os.path.join(out, f)
            if p.returncode != 0 or not os.path.exists(new):
                rows.append((label, f, "SCRIPT FAILED", (p.stderr.strip().splitlines() or ["no output"])[-1][:100])); bad += 1; continue
            if sha(new) == expected.get(f):
                rows.append((label, f, "identical", f"{secs:.0f} s")); continue
            ref = os.path.join(a.reference_dir, f) if a.reference_dir else None
            if ref and os.path.exists(ref):
                d = pixel_diff(new, ref)
                if d is None: rows.append((label, f, "DIFFERS", f"size {Image.open(new).size} vs {Image.open(ref).size}")); bad += 1
                elif d == 0: rows.append((label, f, "same pixels", "file bytes differ (PNG encoder)"))
                else: rows.append((label, f, "DIFFERS", f"{100 * d:.3f}% of pixels")); bad += 1
            else:
                rows.append((label, f, "NOT IDENTICAL", "give --reference-dir to compare pixels")); bad += 1
    w = max(len(r[1]) for r in rows)
    for r in rows: print(f"{r[0]:14s} {r[1]:{w}s}  {r[2]:13s} {r[3]}")
    print(f"\nredrawn figures in {out}")
    print("all figures reproduced" if bad == 0 else f"{bad} figure(s) not reproduced")
    sys.exit(1 if bad else 0)
