#!/usr/bin/env python3

import argparse
import os
import re
from typing import Dict, List, Tuple

import pandas as pd


LONG_COLUMNS = [
    "sample",
    "run_name",
    "run_key",
    "tool",
    "step",
    "benchmark_file",
    "wall_seconds",
    "wall_minutes",
    "wall_hms",
    "cpu_seconds",
    "max_rss_mb",
    "max_vms_mb",
    "max_uss_mb",
    "max_pss_mb",
    "io_in_mb",
    "io_out_mb",
    "mean_load",
]

SUMMARY_COLUMNS = [
    "sample",
    "run_name",
    "run_key",
    "tool",
    "n_steps",
    "steps",
    "wall_seconds_total",
    "wall_minutes_total",
    "cpu_seconds_total",
    "max_rss_mb_peak",
    "max_vms_mb_peak",
    "max_uss_mb_peak",
    "max_pss_mb_peak",
    "io_in_mb_total",
    "io_out_mb_total",
    "mean_load_avg",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize Snakemake benchmark files into comparable tool cost tables.")
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--run_name", required=True, help="Run folder name")
    p.add_argument("--run_key", required=True, help="Canonical run key")
    p.add_argument("--out_long", required=True, help="Output CSV path for long-format step-level costs")
    p.add_argument("--out_summary", required=True, help="Output CSV path for tool-level summarized costs")
    p.add_argument("bench_files", nargs="*", help="Benchmark TSV files from Snakemake benchmark directives")
    return p.parse_args()


def parse_tool_step(bench_path: str) -> Tuple[str, str]:
    name = os.path.basename(bench_path)
    m = re.match(r"^bench_(.+?)__(.+)\.tsv$", name)
    if m:
        return m.group(1), m.group(2)
    stem = os.path.splitext(name)[0]
    return stem, "main"


def as_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def load_one_benchmark(bench_path: str) -> Dict:
    if (not os.path.exists(bench_path)) or os.path.getsize(bench_path) == 0:
        return {}
    try:
        df = pd.read_csv(bench_path, sep="\t")
    except Exception:
        return {}
    if df.empty:
        return {}
    rec = df.iloc[-1].to_dict()
    tool, step = parse_tool_step(bench_path)
    return {
        "tool": tool,
        "step": step,
        "benchmark_file": bench_path,
        "wall_seconds": as_float(rec.get("s")),
        "wall_hms": rec.get("h:m:s", ""),
        "cpu_seconds": as_float(rec.get("cpu_time")),
        "max_rss_mb": as_float(rec.get("max_rss")),
        "max_vms_mb": as_float(rec.get("max_vms")),
        "max_uss_mb": as_float(rec.get("max_uss")),
        "max_pss_mb": as_float(rec.get("max_pss")),
        "io_in_mb": as_float(rec.get("io_in")),
        "io_out_mb": as_float(rec.get("io_out")),
        "mean_load": as_float(rec.get("mean_load")),
    }


def build_long_rows(args: argparse.Namespace) -> List[Dict]:
    rows: List[Dict] = []
    for fp in args.bench_files:
        one = load_one_benchmark(fp)
        if not one:
            continue
        one["sample"] = args.sample
        one["run_name"] = args.run_name
        one["run_key"] = args.run_key
        one["wall_minutes"] = one["wall_seconds"] / 60.0 if pd.notna(one["wall_seconds"]) else float("nan")
        rows.append(one)
    return rows


def empty_df(columns: List[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.out_long), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    rows = build_long_rows(args)
    if rows:
        long_df = pd.DataFrame(rows)
        long_df = long_df[LONG_COLUMNS].sort_values(["tool", "step", "benchmark_file"]).reset_index(drop=True)
    else:
        long_df = empty_df(LONG_COLUMNS)
    long_df.to_csv(args.out_long, index=False)

    if long_df.empty:
        summary_df = empty_df(SUMMARY_COLUMNS)
    else:
        summary_df = (
            long_df.groupby("tool", as_index=False, dropna=False)
            .agg(
                n_steps=("step", "count"),
                steps=("step", lambda s: ";".join(sorted({str(x) for x in s}))),
                wall_seconds_total=("wall_seconds", "sum"),
                cpu_seconds_total=("cpu_seconds", "sum"),
                max_rss_mb_peak=("max_rss_mb", "max"),
                max_vms_mb_peak=("max_vms_mb", "max"),
                max_uss_mb_peak=("max_uss_mb", "max"),
                max_pss_mb_peak=("max_pss_mb", "max"),
                io_in_mb_total=("io_in_mb", "sum"),
                io_out_mb_total=("io_out_mb", "sum"),
                mean_load_avg=("mean_load", "mean"),
            )
            .sort_values(["tool"])
            .reset_index(drop=True)
        )
        summary_df.insert(0, "run_key", args.run_key)
        summary_df.insert(0, "run_name", args.run_name)
        summary_df.insert(0, "sample", args.sample)
        summary_df["wall_minutes_total"] = summary_df["wall_seconds_total"] / 60.0
        summary_df = summary_df[SUMMARY_COLUMNS]
    summary_df.to_csv(args.out_summary, index=False)


if __name__ == "__main__":
    main()
