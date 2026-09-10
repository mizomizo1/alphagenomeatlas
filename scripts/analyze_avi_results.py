#!/usr/bin/env python3
"""Analyze AlphaGenome AVI scoring results for the Evo2 482-variant cohort.

The script keeps the analysis explicit and reproducible:
  - coverage is reported overall and by label/variant type
  - failed rows are part of the denominator
  - AUROC is computed on each predictor's available rows
  - paired AVI vs Evo2 tests are computed only on matched rows
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


# ------------------------------- paired DeLong test for correlated ROC AUCs
# Adapted from the standard fast DeLong implementation by Sun and Xu.

def compute_midrank(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    order = np.argsort(x)
    sorted_x = x[order]
    n = len(x)
    ranks = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and sorted_x[j] == sorted_x[i]:
            j += 1
        ranks[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(n, dtype=float)
    out[order] = ranks
    return out


def fast_delong(predictions_sorted_transposed: np.ndarray, n_pos: int) -> tuple[np.ndarray, np.ndarray]:
    m = n_pos
    n = predictions_sorted_transposed.shape[1] - m
    if m <= 0 or n <= 0:
        raise ValueError("DeLong test requires at least one positive and one negative sample")

    pos = predictions_sorted_transposed[:, :m]
    neg = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty((k, m), dtype=float)
    ty = np.empty((k, n), dtype=float)
    tz = np.empty((k, m + n), dtype=float)
    for r in range(k):
        tx[r, :] = compute_midrank(pos[r, :])
        ty[r, :] = compute_midrank(neg[r, :])
        tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])

    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    cov = sx / m + sy / n
    cov = np.atleast_2d(cov)
    return aucs, cov


def delong_two_model_test(y_true: np.ndarray, score_a: np.ndarray, score_b: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    score_a = np.asarray(score_a).astype(float)
    score_b = np.asarray(score_b).astype(float)
    ok = np.isfinite(score_a) & np.isfinite(score_b) & np.isin(y_true, [0, 1])
    y_true, score_a, score_b = y_true[ok], score_a[ok], score_b[ok]

    order = np.argsort(-y_true)  # positives first
    y_sorted = y_true[order]
    preds = np.vstack([score_a[order], score_b[order]])
    n_pos = int(y_sorted.sum())
    n_neg = int(len(y_sorted) - n_pos)
    aucs, cov = fast_delong(preds, n_pos)
    contrast = np.array([[1.0, -1.0]])
    var = float(contrast @ cov @ contrast.T)
    if var <= 0:
        z = np.nan
        p = np.nan
    else:
        z = float(abs(aucs[0] - aucs[1]) / np.sqrt(var))
        p = float(2 * stats.norm.sf(z))
    return {
        "n": int(len(y_true)),
        "n_pathogenic": n_pos,
        "n_benign": n_neg,
        "auc_a": float(aucs[0]),
        "auc_b": float(aucs[1]),
        "delta_auc_a_minus_b": float(aucs[0] - aucs[1]),
        "z": z,
        "p_two_sided": p,
    }


# --------------------------------------------------------------- utilities

def safe_auc(y: pd.Series, score: pd.Series) -> float:
    y = pd.to_numeric(y, errors="coerce")
    s = pd.to_numeric(score, errors="coerce")
    ok = y.isin([0, 1]) & s.notna()
    if ok.sum() < 3 or y[ok].nunique() < 2:
        return np.nan
    return float(roc_auc_score(y[ok].astype(int), s[ok].astype(float)))


def predictor_row(df: pd.DataFrame, name: str, score_col: str) -> dict[str, Any]:
    y = pd.to_numeric(df["label"], errors="coerce")
    s = pd.to_numeric(df[score_col], errors="coerce")
    ok = y.isin([0, 1]) & s.notna()
    d = df.loc[ok]
    return {
        "predictor": name,
        "score_column": score_col,
        "n_available": int(ok.sum()),
        "n_pathogenic": int((d["label"].astype(int) == 1).sum()),
        "n_benign": int((d["label"].astype(int) == 0).sum()),
        "auroc": safe_auc(df["label"], df[score_col]),
    }


def coverage_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label in [0, 1]:
        for vt in ["ALL"] + sorted(df["variant_type"].dropna().astype(str).unique().tolist()):
            sub = df[df["label"].astype(int).eq(label)]
            if vt != "ALL":
                sub = sub[sub["variant_type"].astype(str).eq(vt)]
            n = len(sub)
            avail = int(sub["avi_available"].sum()) if n else 0
            rows.append({
                "label": label,
                "variant_type": vt,
                "n": n,
                "n_avi_available": avail,
                "n_avi_missing_or_failed": n - avail,
                "avi_coverage_rate": avail / n if n else np.nan,
            })
    return pd.DataFrame(rows)


def fisher_coverage_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = ["ALL"] + [f"variant_type={v}" for v in sorted(df["variant_type"].dropna().astype(str).unique())]
    for subset in subsets:
        sub = df.copy()
        if subset != "ALL":
            vt = subset.split("=", 1)[1]
            sub = sub[sub["variant_type"].astype(str).eq(vt)]
        pth = sub[sub["label"].astype(int).eq(1)]
        ben = sub[sub["label"].astype(int).eq(0)]
        pa = int(pth["avi_available"].sum())
        pm = int(len(pth) - pa)
        ba = int(ben["avi_available"].sum())
        bm = int(len(ben) - ba)
        if len(pth) == 0 or len(ben) == 0:
            p = np.nan
        else:
            p = float(stats.fisher_exact([[pa, pm], [ba, bm]], alternative="two-sided").pvalue)
        rows.append({
            "subset": subset,
            "pathogenic_available": pa,
            "pathogenic_missing_or_failed": pm,
            "benign_available": ba,
            "benign_missing_or_failed": bm,
            "pathogenic_coverage_rate": pa / len(pth) if len(pth) else np.nan,
            "benign_coverage_rate": ba / len(ben) if len(ben) else np.nan,
            "fisher_exact_two_sided_p": p,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored-tsv", required=True, help="01_all482_with_avi.tsv from score_avi_all482.py")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.scored_tsv, sep="\t")

    df["label"] = pd.to_numeric(df["label"], errors="coerce").astype("Int64")
    df["avi_score"] = pd.to_numeric(df["avi_score"], errors="coerce")
    df["evo2_delta_7b"] = pd.to_numeric(df["evo2_delta_7b"], errors="coerce")
    df["evo2_neg_delta_7b"] = -df["evo2_delta_7b"]
    df["evo2_abs_delta_7b"] = df["evo2_delta_7b"].abs()
    df["avi_available"] = df["avi_query_status"].eq("success") & df["avi_score"].notna()

    coverage = coverage_summary(df)
    coverage.to_csv(outdir / "05_coverage_by_label_and_variant_type.tsv", sep="\t", index=False)

    tests = fisher_coverage_tests(df)
    tests.to_csv(outdir / "06_label_coverage_test.tsv", sep="\t", index=False)

    pred_rows = [
        predictor_row(df[df["avi_available"]].copy(), "AlphaGenome AVI", "avi_score"),
        predictor_row(df, "Evo2 -delta", "evo2_neg_delta_7b"),
        predictor_row(df, "Evo2 abs(delta)", "evo2_abs_delta_7b"),
    ]
    pred = pd.DataFrame(pred_rows)
    pred.to_csv(outdir / "07_predictor_auroc_self_coverage.tsv", sep="\t", index=False)

    matched = df[df["avi_available"] & df["evo2_delta_7b"].notna() & df["label"].isin([0, 1])].copy()
    delong_rows = []
    if matched["label"].nunique() == 2:
        for evo_name, evo_col in [("Evo2 -delta", "evo2_neg_delta_7b"), ("Evo2 abs(delta)", "evo2_abs_delta_7b")]:
            r = delong_two_model_test(
                matched["label"].to_numpy(),
                matched["avi_score"].to_numpy(),
                matched[evo_col].to_numpy(),
            )
            r.update({"model_a": "AlphaGenome AVI", "model_b": evo_name, "matched_subset": "AVI_available_and_Evo2_available"})
            delong_rows.append(r)
    delong = pd.DataFrame(delong_rows)
    delong.to_csv(outdir / "08_delong_avi_vs_evo2.tsv", sep="\t", index=False)

    logo_rows = []
    for gene in sorted(matched["gene"].dropna().astype(str).unique()):
        sub = matched[~matched["gene"].astype(str).eq(gene)].copy()
        if sub["label"].nunique() < 2:
            continue
        logo_rows.append({
            "left_out_gene": gene,
            "n": len(sub),
            "n_pathogenic": int((sub["label"].astype(int) == 1).sum()),
            "n_benign": int((sub["label"].astype(int) == 0).sum()),
            "avi_auroc": safe_auc(sub["label"], sub["avi_score"]),
            "evo2_neg_delta_auroc": safe_auc(sub["label"], sub["evo2_neg_delta_7b"]),
            "evo2_abs_delta_auroc": safe_auc(sub["label"], sub["evo2_abs_delta_7b"]),
        })
    logo = pd.DataFrame(logo_rows)
    logo.to_csv(outdir / "09_leave_one_gene_out_auc_stability.tsv", sep="\t", index=False)

    summary = {
        "analysis_finished_utc": datetime.now(timezone.utc).isoformat(),
        "input_scored_tsv": str(Path(args.scored_tsv).resolve()),
        "n_total_rows": int(len(df)),
        "avi_status_counts": df["avi_query_status"].value_counts(dropna=False).to_dict(),
        "variant_type_counts": df["variant_type"].value_counts(dropna=False).to_dict(),
        "n_avi_available": int(df["avi_available"].sum()),
        "n_matched_avi_evo2": int(len(matched)),
        "main_outputs": {
            "coverage": "05_coverage_by_label_and_variant_type.tsv",
            "coverage_tests": "06_label_coverage_test.tsv",
            "auroc": "07_predictor_auroc_self_coverage.tsv",
            "delong": "08_delong_avi_vs_evo2.tsv",
            "leave_one_gene_out": "09_leave_one_gene_out_auc_stability.tsv",
        },
        "interpretation_note": (
            "Interpret AVI as a predictor available on its successful query subset. "
            "If indels are missing, full-cohort comparisons against Evo2 must discuss coverage differences."
        ),
    }
    with open(outdir / "10_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("Saved outputs in", outdir)
    print("\nAVI status counts:")
    print(df["avi_query_status"].value_counts(dropna=False))
    print("\nCoverage tests:")
    print(tests)
    print("\nAUROC:")
    print(pred)
    if len(delong):
        print("\nPaired DeLong tests:")
        print(delong)


if __name__ == "__main__":
    main()
