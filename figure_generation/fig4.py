#Figure 4: individual clinical variables (IPF-LASSO, rescued results)
#A Shapley values, B additivity; both drawn by variable_importance.py
#Usage, from the project root: python figure_generation/fig4.py --results-dir results --out-dir figures
import argparse, os
import variable_importance as vi

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results"); ap.add_argument("--out-dir", default="figures")
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    R, curves, T, P, A, B = vi.collect(a.results_dir)
    vi.main_figure(R, P, A, B, os.path.join(a.out_dir, "fig4.png"))
