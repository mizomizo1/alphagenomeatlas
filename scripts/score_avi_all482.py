#!/usr/bin/env python3
"""Score a local Evo2 482-variant TSV with AlphaGenome Atlas AVI_SCORE.

No data and no API key are stored in this repository.
The API key is read from an environment variable, and all input/output paths are CLI arguments.

Important coordinate convention:
  - Query position is `pos_hg38` from the input table.
  - `pos_hg38` is expected to be 1-based.
  - The original hg19 coordinate is kept only for audit columns.
  - No ref/alt swap rescue is performed.
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from alphagenome.atlas import atlas as alphagenome_atlas
from alphagenome.data import genome


REQUIRED_COLUMNS = [
    "record_id",
    "label",
    "chrom",
    "pos_hg19",
    "pos_hg38",
    "ref",
    "alt",
    "variant_type",
    "evo2_delta_7b",
]


def clean_chrom(x: Any) -> str:
    s = str(x).strip()
    if not s:
        raise ValueError("empty chromosome")
    return s if s.startswith("chr") else "chr" + s


def clean_allele(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    if s in {"-", ".", "NAN", "NONE"}:
        return ""
    return s


def flatten_first_numeric(x: Any) -> float:
    """Return the first finite numeric value from ndarray/sparse/list-like objects."""
    if x is None:
        return float("nan")
    try:
        if hasattr(x, "toarray"):
            arr = x.toarray()
        else:
            arr = np.asarray(x)
        arr = arr.astype(float).ravel()
        arr = arr[np.isfinite(arr)]
        return float(arr[0]) if arr.size else float("nan")
    except Exception:
        return float("nan")


def simplify_metadata_value(v: Any) -> dict[str, Any]:
    d: dict[str, Any] = {"class": type(v).__name__}
    for attr in ["is_signed", "recommended", "description", "version", "name"]:
        if hasattr(v, attr):
            try:
                val = getattr(v, attr)
                d[attr] = val if isinstance(val, (str, int, float, bool, type(None))) else str(val)
            except Exception:
                d[attr] = "<unreadable>"
    attrs = [a for a in dir(v) if not a.startswith("_")]
    d["available_attrs"] = ",".join(attrs)
    for attr in ["tracks", "track_metadata", "track_names", "output_tracks"]:
        if hasattr(v, attr):
            try:
                d["track_field"] = attr
                d["n_tracks"] = len(getattr(v, attr))
            except Exception:
                d["track_field"] = attr
                d["n_tracks"] = None
            break
    return d


def make_variant(row: pd.Series) -> genome.Variant:
    chrom = clean_chrom(row["chrom"])
    pos = int(row["pos_hg38"])
    ref = clean_allele(row["ref"])
    alt = clean_allele(row["alt"])
    return genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
    )


def query_one(client: Any, row: pd.Series, scorer: str) -> dict[str, Any]:
    v = make_variant(row)
    out = {
        "query_chrom": v.chromosome,
        "query_pos_hg38_1based": v.position,
        "query_ref": v.reference_bases,
        "query_alt": v.alternate_bases,
        "query_convention": "genome.Variant with 1-based hg38 position from pos_hg38",
        "avi_query_status": "failed",
        "avi_score": np.nan,
        "avi_quantile": np.nan,
        "avi_error_type": "",
        "avi_error_message": "",
    }
    try:
        result = client.query_variant(v, requested_scorers=[scorer])
        if scorer not in result:
            raise KeyError(f"{scorer} not returned. returned={list(result.keys())}")
        adata = result[scorer]
        out["avi_score"] = flatten_first_numeric(getattr(adata, "X", None))
        layers = getattr(adata, "layers", {})
        if layers is not None and "quantiles" in layers:
            out["avi_quantile"] = flatten_first_numeric(layers["quantiles"])
        if np.isfinite(out["avi_score"]):
            out["avi_query_status"] = "success"
        else:
            out["avi_query_status"] = "missing_score"
            out["avi_error_message"] = "query returned but no finite AVI score was found"
    except Exception as e:
        out["avi_error_type"] = type(e).__name__
        out["avi_error_message"] = str(e)
    return out


def validate_input(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"input TSV missing required columns: {missing}")
    if df["record_id"].duplicated().any():
        dup = df.loc[df["record_id"].duplicated(), "record_id"].head().tolist()
        raise ValueError(f"record_id must be unique; duplicated examples: {dup}")
    if df["pos_hg38"].isna().any():
        raise ValueError("pos_hg38 contains missing values; all variants must be pre-lifted to hg38")
    same = (pd.to_numeric(df["pos_hg19"], errors="coerce") == pd.to_numeric(df["pos_hg38"], errors="coerce")).sum()
    if same:
        print(f"WARNING: {same} rows have pos_hg19 == pos_hg38. Verify liftover provenance.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-tsv", required=True, help="local all_482_for_avi.tsv")
    ap.add_argument("--outdir", required=True, help="local output directory")
    ap.add_argument("--api-key-env", default="ALPHA_GENOME_API_KEY")
    ap.add_argument("--scorer", default="AVI_SCORE")
    ap.add_argument("--sleep-sec", type=float, default=0.0)
    ap.add_argument("--limit", type=int, default=None, help="debug only: score first N rows")
    ap.add_argument("--resume", action="store_true", help="reuse existing successful rows in output TSV")
    args = ap.parse_args()

    api_key = os.environ.get(args.api_key_env)
    if not api_key:
        raise RuntimeError(f"Set {args.api_key_env} in your local environment. Do not put API keys in code.")

    input_tsv = Path(args.input_tsv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_tsv, sep="\t")
    if args.limit:
        df = df.head(args.limit).copy()
    validate_input(df)

    output_tsv = outdir / "01_all482_with_avi.tsv"
    previous: dict[str, dict[str, Any]] = {}
    if args.resume and output_tsv.exists():
        prev = pd.read_csv(output_tsv, sep="\t")
        if "record_id" in prev.columns and "avi_query_status" in prev.columns:
            ok = prev[prev["avi_query_status"].eq("success")]
            previous = ok.set_index("record_id").to_dict(orient="index")
            print(f"Resume: reusing {len(previous)} previous successful records")

    client = alphagenome_atlas.create(api_key, timeout=60)
    scorer_metadata = client.scorer_metadata()
    if args.scorer not in scorer_metadata:
        raise ValueError(f"{args.scorer} not found in scorer metadata: {list(scorer_metadata.keys())}")

    meta_rows = []
    for k, v in scorer_metadata.items():
        row = {"scorer": k}
        row.update(simplify_metadata_value(v))
        meta_rows.append(row)
    pd.DataFrame(meta_rows).to_csv(outdir / "04_scorer_metadata.tsv", sep="\t", index=False)

    run_metadata = {
        "run_started_utc": datetime.now(timezone.utc).isoformat(),
        "input_tsv_name": input_tsv.name,
        "input_tsv_absolute_path": str(input_tsv.resolve()),
        "n_input_rows": int(len(df)),
        "scorer": args.scorer,
        "api_key_source": f"environment variable {args.api_key_env}",
        "query_coordinate": "pos_hg38",
        "query_coordinate_convention": "1-based hg38",
        "ref_alt_swap_rescue": False,
        "hg19_query_used": False,
        "failed_rows_retained": True,
    }
    try:
        import alphagenome
        run_metadata["alphagenome_package_version"] = getattr(alphagenome, "__version__", "unknown")
    except Exception:
        run_metadata["alphagenome_package_version"] = "unknown"

    rows = []
    for _, row in tqdm(df.iterrows(), total=len(df)):
        rid = str(row["record_id"])
        if rid in previous:
            q = previous[rid]
        else:
            q = query_one(client, row, args.scorer)
            if args.sleep_sec:
                time.sleep(args.sleep_sec)
        merged = row.to_dict()
        merged.update(q)
        rows.append(merged)

        # checkpoint every 25 rows, while preserving failed rows
        if len(rows) % 25 == 0:
            pd.DataFrame(rows).to_csv(output_tsv, sep="\t", index=False)

    result = pd.DataFrame(rows)
    result.to_csv(output_tsv, sep="\t", index=False)
    result.loc[~result["avi_query_status"].eq("success")].to_csv(
        outdir / "02_failed_or_missing_avi.tsv", sep="\t", index=False
    )

    run_metadata["run_finished_utc"] = datetime.now(timezone.utc).isoformat()
    run_metadata["status_counts"] = result["avi_query_status"].value_counts(dropna=False).to_dict()
    run_metadata["variant_type_counts"] = result["variant_type"].value_counts(dropna=False).to_dict()
    with open(outdir / "03_run_metadata.json", "w") as f:
        json.dump(run_metadata, f, indent=2, ensure_ascii=False)

    print("Saved:")
    print(" ", output_tsv)
    print(" ", outdir / "02_failed_or_missing_avi.tsv")
    print(" ", outdir / "03_run_metadata.json")
    print(" ", outdir / "04_scorer_metadata.tsv")
    print(result["avi_query_status"].value_counts(dropna=False))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise
