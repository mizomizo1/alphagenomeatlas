#!/usr/bin/env python3
"""Validate AlphaGenome Atlas AVI queries for INDEL representations.

This script is intentionally data/API-key free. Run it locally or in Colab with:
  export ALPHA_GENOME_API_KEY="..."
  python scripts/validate_indel_query_representations.py \
    --input-tsv /path/to/all_482_for_avi.tsv \
    --outdir /path/to/indel_validation_outputs \
    --use-ucsc-reference

Why this exists:
  Atlas AVI scoring of SNVs may work while INDELs fail if the query
  representation does not match the Atlas' canonical representation.
  This script tests multiple representations of each INDEL and records every
  success/failure without deleting failed rows.

Official AlphaGenome conventions:
  - genome.Variant.position is 1-based.
  - Interval coordinates are 0-based half-open.
  - Variant reference/alternate strings can be empty, but INDEL normalization
    requires reference sequence context.
"""

from __future__ import annotations

import argparse
import json
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from alphagenome.atlas import atlas as alphagenome_atlas
from alphagenome.data import genome


REQUIRED_COLUMNS = [
    "record_id", "label", "chrom", "pos_hg19", "pos_hg38", "ref", "alt",
    "variant_type", "spdi_hg38_0based"
]


def clean_chrom(x: Any) -> str:
    s = str(x).strip()
    return s if s.startswith("chr") else f"chr{s}"


def clean_allele(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    if s in {"-", ".", "NAN", "NONE"}:
        return ""
    return s


def variant_key(v: genome.Variant) -> str:
    return f"{v.chromosome}:{v.position}:{v.reference_bases}>{v.alternate_bases}"


def flatten_first_numeric(x: Any) -> float:
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


@dataclass
class Candidate:
    strategy: str
    variant: genome.Variant
    note: str = ""


class ReferenceExtractor:
    """Small reference extractor backed by local FASTA or UCSC REST.

    Local FASTA is preferable for strict reproducibility.
    UCSC REST is convenient for quick validation of the 23 INDELs.
    """

    def __init__(self, hg38_fasta: str | None = None, use_ucsc: bool = False):
        self.hg38_fasta = hg38_fasta
        self.use_ucsc = use_ucsc
        self._cache: dict[str, str] = {}
        self._fasta = None
        if hg38_fasta:
            try:
                from pyfaidx import Fasta
            except ImportError as e:
                raise RuntimeError("Install pyfaidx to use --hg38-fasta") from e
            self._fasta = Fasta(hg38_fasta, rebuild=False)

    def available(self) -> bool:
        return self._fasta is not None or self.use_ucsc

    def extract(self, interval: genome.Interval) -> str:
        if interval.end < interval.start:
            raise ValueError(f"bad interval: {interval}")
        if interval.end == interval.start:
            return ""

        chrom = clean_chrom(interval.chromosome)
        start = int(interval.start)
        end = int(interval.end)
        key = f"{chrom}:{start}-{end}"

        if key in self._cache:
            return self._cache[key]

        if self._fasta is not None:
            if chrom in self._fasta:
                seq = str(self._fasta[chrom][start:end]).upper()
            else:
                seq = str(self._fasta[chrom.replace("chr", "")][start:end]).upper()
        elif self.use_ucsc:
            params = urllib.parse.urlencode({
                "genome": "hg38",
                "chrom": chrom,
                "start": start,
                "end": end,
            })
            url = f"https://api.genome.ucsc.edu/getData/sequence?{params}"
            with urllib.request.urlopen(url, timeout=30) as fh:
                data = json.load(fh)
            seq = str(data.get("dna", "")).upper()
            if len(seq) != end - start:
                raise ValueError(f"UCSC returned length {len(seq)} for {key}")
        else:
            raise RuntimeError("No reference source is available.")

        self._cache[key] = seq
        return seq


def parse_spdi_to_variant(spdi: str) -> genome.Variant | None:
    """Parse '<chrom>:<0based_pos>:<deleted>:<inserted>' into genome.Variant."""
    if not isinstance(spdi, str) or not spdi.strip():
        return None
    parts = spdi.strip().split(":")
    if len(parts) != 4:
        return None
    chrom, pos0, deleted, inserted = parts
    return genome.Variant(
        chromosome=clean_chrom(chrom),
        position=int(pos0) + 1,
        reference_bases=deleted.upper(),
        alternate_bases=inserted.upper(),
    )


def make_candidates(row: pd.Series, refx: ReferenceExtractor | None) -> list[Candidate]:
    chrom = clean_chrom(row["chrom"])
    pos = int(row["pos_hg38"])
    ref = clean_allele(row["ref"])
    alt = clean_allele(row["alt"])

    candidates: list[Candidate] = []

    direct = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
    candidates.append(Candidate("table_direct_pos_hg38_1based", direct, "pos_hg38 with '-' converted to empty allele"))

    spdi_v = parse_spdi_to_variant(str(row.get("spdi_hg38_0based", "")))
    if spdi_v is not None:
        candidates.append(Candidate("spdi_direct_pos0_plus1", spdi_v, "SPDI pos is 0-based; query position = pos0 + 1"))

    if pos > 1:
        candidates.append(Candidate("table_direct_pos_minus1_sanity", genome.Variant(chrom, pos - 1, ref, alt), "diagnostic only"))
    candidates.append(Candidate("table_direct_pos_plus1_sanity", genome.Variant(chrom, pos + 1, ref, alt), "diagnostic only"))

    if refx is not None and refx.available() and pos > 1:
        try:
            prev_base = refx.extract(genome.Interval(chrom, pos - 2, pos - 1))
            if len(prev_base) == 1 and set(prev_base).issubset(set("ACGTN")):
                padded = genome.Variant(
                    chromosome=chrom,
                    position=pos - 1,
                    reference_bases=prev_base + ref,
                    alternate_bases=prev_base + alt,
                )
                candidates.append(Candidate("left_padded_previous_base", padded, f"prev_base={prev_base}"))
        except Exception as e:
            candidates.append(Candidate(
                "left_padded_previous_base_ERROR",
                direct,
                f"candidate construction failed: {type(e).__name__}: {e}",
            ))

    if refx is not None and refx.available():
        for base_name, base_v in [
            ("normalize_table_direct", direct),
            ("normalize_spdi_direct", spdi_v),
        ]:
            if base_v is None:
                continue
            try:
                norm = genome.normalize_variant(base_v, refx)
                candidates.append(Candidate(base_name, norm, f"normalized from {variant_key(base_v)}"))
            except Exception as e:
                candidates.append(Candidate(
                    base_name + "_ERROR",
                    base_v,
                    f"normalization failed: {type(e).__name__}: {e}",
                ))

    seen: set[str] = set()
    unique: list[Candidate] = []
    for c in candidates:
        key = variant_key(c.variant)
        if c.strategy.endswith("_ERROR") or key not in seen:
            unique.append(c)
            seen.add(key)
    return unique


def query_avi(client: Any, variant: genome.Variant, scorer: str) -> dict[str, Any]:
    out = {
        "query_variant": variant_key(variant),
        "query_chrom": variant.chromosome,
        "query_pos_1based": int(variant.position),
        "query_ref": variant.reference_bases,
        "query_alt": variant.alternate_bases,
        "query_is_snv": bool(variant.is_snv),
        "query_is_indel": bool(variant.is_indel),
        "status": "failed",
        "avi_score": np.nan,
        "avi_quantile": np.nan,
        "error_type": "",
        "error_message": "",
    }
    try:
        res = client.query_variant(variant, requested_scorers=[scorer])
        if scorer not in res:
            raise KeyError(f"{scorer} not returned; returned={list(res)}")
        ad = res[scorer]
        out["avi_score"] = flatten_first_numeric(getattr(ad, "X", None))
        layers = getattr(ad, "layers", {})
        if layers is not None and "quantiles" in layers:
            out["avi_quantile"] = flatten_first_numeric(layers["quantiles"])
        out["status"] = "success" if np.isfinite(out["avi_score"]) else "missing_score"
    except Exception as e:
        out["error_type"] = type(e).__name__
        out["error_message"] = str(e)
    return out


def _variant_str_is_snv(s: str) -> bool:
    try:
        ref_alt = s.split(":")[-1]
        ref, alt = ref_alt.split(">")
        return len(ref) == 1 and len(alt) == 1
    except Exception:
        return False


def interval_probe(client: Any, row: pd.Series, scorer: str, flank: int) -> dict[str, Any]:
    chrom = clean_chrom(row["chrom"])
    pos0 = int(row["pos_hg38"]) - 1
    interval = genome.Interval(chrom, max(0, pos0 - flank), pos0 + flank + 1)
    out = {
        "record_id": row["record_id"],
        "interval": str(interval),
        "status": "failed",
        "n_returned": 0,
        "n_indel_returned": 0,
        "returned_examples": "",
        "error_type": "",
        "error_message": "",
    }
    try:
        res = client.query_interval(interval, requested_scorers=[scorer], progress_bar=False, max_workers=1)
        if scorer not in res:
            raise KeyError(f"{scorer} not returned; returned={list(res)}")
        obs = getattr(res[scorer], "obs", pd.DataFrame())
        variants = [str(v) for v in obs.get("variant", [])] if obs is not None and "variant" in obs else []
        out["status"] = "success"
        out["n_returned"] = len(variants)
        out["n_indel_returned"] = sum(1 for v in variants if not _variant_str_is_snv(v))
        out["returned_examples"] = ";".join(variants[:20])
    except Exception as e:
        out["error_type"] = type(e).__name__
        out["error_message"] = str(e)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-tsv", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--api-key-env", default="ALPHA_GENOME_API_KEY")
    p.add_argument("--scorer", default="AVI_SCORE")
    p.add_argument("--hg38-fasta", default=None, help="Optional local GRCh38/hg38 FASTA for padding/normalization.")
    p.add_argument("--use-ucsc-reference", action="store_true", help="Fetch small hg38 reference snippets from UCSC REST.")
    p.add_argument("--sleep-sec", type=float, default=0.0)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--interval-probe", action="store_true", help="Also query a small interval around each INDEL.")
    p.add_argument("--interval-flank", type=int, default=2)
    args = p.parse_args()

    api_key = os.environ.get(args.api_key_env)
    if not api_key:
        raise RuntimeError(f"Set {args.api_key_env}. Do not put API keys in code.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_tsv, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    indels = df[df["variant_type"].astype(str).str.upper().eq("INDEL")].copy()
    if args.limit:
        indels = indels.head(args.limit).copy()

    refx = ReferenceExtractor(args.hg38_fasta, args.use_ucsc_reference)
    if not refx.available():
        print("Reference source not set: only direct/SPDI/off-by-one candidates will be tested.")
        print("For anchored/normalized candidates, pass --hg38-fasta or --use-ucsc-reference.")

    client = alphagenome_atlas.create(api_key, timeout=60)
    meta = client.scorer_metadata()
    if args.scorer not in meta:
        raise ValueError(f"{args.scorer} not found. Available scorers: {list(meta)}")

    run_meta = {
        "run_started_utc": datetime.now(timezone.utc).isoformat(),
        "input_tsv": str(Path(args.input_tsv).resolve()),
        "n_input_rows": int(len(df)),
        "n_indels_tested": int(len(indels)),
        "scorer": args.scorer,
        "api_key_source": f"environment variable {args.api_key_env}",
        "reference_source": args.hg38_fasta or ("UCSC REST hg38" if args.use_ucsc_reference else "none"),
        "note": "Tests multiple INDEL representations against Atlas AVI_SCORE; failed rows retained.",
    }

    candidate_rows = []
    for _, row in tqdm(indels.iterrows(), total=len(indels), desc="INDEL candidate queries"):
        for cand in make_candidates(row, refx):
            q = query_avi(client, cand.variant, args.scorer)
            q.update({
                "record_id": row["record_id"],
                "label": row["label"],
                "gene": row.get("gene", ""),
                "chrom": row["chrom"],
                "pos_hg19": row["pos_hg19"],
                "pos_hg38": row["pos_hg38"],
                "ref": row["ref"],
                "alt": row["alt"],
                "spdi_hg38_0based": row.get("spdi_hg38_0based", ""),
                "vep_consequence": row.get("vep_consequence", ""),
                "candidate_strategy": cand.strategy,
                "candidate_note": cand.note,
            })
            candidate_rows.append(q)
            if args.sleep_sec:
                time.sleep(args.sleep_sec)

    cand_df = pd.DataFrame(candidate_rows)
    cand_path = outdir / "indel_candidate_query_results.tsv"
    cand_df.to_csv(cand_path, sep="\t", index=False)

    probe_df = pd.DataFrame()
    if args.interval_probe:
        probe_rows = []
        for _, row in tqdm(indels.iterrows(), total=len(indels), desc="Interval probes"):
            probe_rows.append(interval_probe(client, row, args.scorer, args.interval_flank))
            if args.sleep_sec:
                time.sleep(args.sleep_sec)
        probe_df = pd.DataFrame(probe_rows)
        probe_df.to_csv(outdir / "indel_interval_probe_results.tsv", sep="\t", index=False)

    summary = []
    if not cand_df.empty:
        for strategy, g in cand_df.groupby("candidate_strategy", dropna=False):
            summary.append({
                "candidate_strategy": strategy,
                "n_queries": int(len(g)),
                "n_success": int(g["status"].eq("success").sum()),
                "success_rate": float(g["status"].eq("success").mean()),
                "n_unique_records_success": int(g.loc[g["status"].eq("success"), "record_id"].nunique()),
            })
    summary_df = pd.DataFrame(summary).sort_values(["n_success", "candidate_strategy"], ascending=[False, True])
    summary_df.to_csv(outdir / "indel_candidate_strategy_summary.tsv", sep="\t", index=False)

    run_meta["run_finished_utc"] = datetime.now(timezone.utc).isoformat()
    run_meta["candidate_strategy_summary"] = summary
    if not probe_df.empty:
        run_meta["interval_probe_n_indel_returned_total"] = int(probe_df["n_indel_returned"].sum())
    with open(outdir / "indel_validation_run_metadata.json", "w") as f:
        json.dump(run_meta, f, indent=2, ensure_ascii=False)

    print("\nCandidate strategy summary:")
    print(summary_df.to_string(index=False))
    if not probe_df.empty:
        print("\nInterval probe summary:")
        print(probe_df[["status", "n_returned", "n_indel_returned"]].describe(include="all").to_string())
    print("\nSaved:")
    print(" ", cand_path)
    print(" ", outdir / "indel_candidate_strategy_summary.tsv")
    if args.interval_probe:
        print(" ", outdir / "indel_interval_probe_results.tsv")
    print(" ", outdir / "indel_validation_run_metadata.json")


if __name__ == "__main__":
    main()
