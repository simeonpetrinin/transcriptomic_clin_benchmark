#Clinical and fused model failures by number of clinical variables (main TCGA run, before rescue)
import sys, glob, os, re
import pandas as pd

cols = ["n_vars", "is_anchor", "clin_separated", "fused_error", "fused_path_short", "fused_at_boundary"]
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*_results.csv"))):
    if not re.fullmatch(r"[A-Z]+_results\.csv", os.path.basename(f)):
        continue  #main run only, not robust or published files
    cancer = os.path.basename(f).split("_")[0]
    d = pd.read_csv(f, usecols=cols, dtype=str, low_memory=False)
    d = d[d["is_anchor"] != "TRUE"]
    d["n_vars"] = d["n_vars"].astype(int)
    d["clin_sep"] = d["clin_separated"] == "TRUE"
    d["fused_fail"] = d["fused_error"].notna()
    d["fused_trunc"] = (d["fused_path_short"] == "TRUE") | (d["fused_at_boundary"] == "TRUE")
    d["cancer"] = cancer
    rows.append(d[["cancer", "n_vars", "clin_sep", "fused_fail", "fused_trunc"]])
a = pd.concat(rows)

pd.set_option("display.width", 200)
print("pooled over cancers: share of fits (%) by number of clinical variables")
g = a.groupby("n_vars").agg(fits=("clin_sep", "size"), cancers=("cancer", "nunique"),
                            clin_separation=("clin_sep", "mean"), fused_failed=("fused_fail", "mean"),
                            fused_truncated=("fused_trunc", "mean"))
for c in ["clin_separation", "fused_failed", "fused_truncated"]:
    g[c] = (100 * g[c]).round(1)
print(g.to_string())
print()
print("per cancer: % clinical separation and % fused failed, smallest vs largest subsets")
out = []
for cancer, d in a.groupby("cancer"):
    kmax = d["n_vars"].max()
    lo, hi = d[d["n_vars"] == 1], d[d["n_vars"] == kmax]
    out.append({"cancer": cancer, "k": kmax,
                "sep_k1": round(100 * lo["clin_sep"].mean(), 1), "sep_kmax": round(100 * hi["clin_sep"].mean(), 1),
                "fail_k1": round(100 * lo["fused_fail"].mean(), 1), "fail_kmax": round(100 * hi["fused_fail"].mean(), 1),
                "fail_all": round(100 * d["fused_fail"].mean(), 1)})
print(pd.DataFrame(out).sort_values("k").to_string(index=False))
