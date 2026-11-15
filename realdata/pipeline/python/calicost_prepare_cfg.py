#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calicost_prepare_cfg.py
-----------------------
A script to generate configuration files for CalicoST Step 1 and Step 2.

1. Step 1: Generates config using explicit reference paths passed from Snakemake.
2. Step 2: Reads a physical YAML template, performs variable substitution 
   (e.g., __CALICOST_DIR__), and updates run-specific paths.
"""

import os
import sys
import argparse
import datetime
import warnings
from copy import deepcopy
from typing import Any, Dict, List

# --- Dependency Check ---
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required. Please install it with `pip install pyyaml`.", file=sys.stderr)
    sys.exit(1)

# --- YAML Formatting Helpers ---
class _FlowSeq(list):
    """Marker list that dumps in flow style: [1, 2, 3]"""
    pass

def _repr_flow_seq(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(_FlowSeq, _repr_flow_seq, Dumper=yaml.SafeDumper)

# --- Utility Functions ---
def now() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def read_yaml(fp: str) -> Dict[str, Any]:
    with open(fp, "r") as f:
        return yaml.safe_load(f) or {}

def write_yaml(obj, fp, header_comment=None, flow_keys=None):
    os.makedirs(os.path.dirname(os.path.abspath(fp)), exist_ok=True)
    data = deepcopy(obj)

    # Wrap specific keys (like chromosomes) as flow-style lists [1, 2, 3]
    if flow_keys:
        def wrap_flow(x):
            if isinstance(x, dict):
                for k, v in list(x.items()):
                    if k in flow_keys and isinstance(v, list):
                        x[k] = _FlowSeq(v)
                    else:
                        wrap_flow(v)
            elif isinstance(x, list):
                for i in range(len(x)):
                    wrap_flow(x[i])
        wrap_flow(data)

    with open(fp, "w") as f:
        if header_comment:
            f.write(f"# {header_comment}\n")
            f.write(f"# Generated at: {now()}\n")
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)

def ensure_dir(path: str) -> None:
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

def eagerly_norm(path: str) -> str:
    return os.path.abspath(path) if path else path

# --- Logic for Step 1 ---
STEP1_DEFAULTS = {
    "samtools_sorting_mem": "4G",
    "UMItag": "Auto",
    "cellTAG": "CB",
    "nthreads_cellsnplite": 20,
    "chromosomes": list(range(1, 23)),
}

def generate_step1_config(
    pipeline_cfg: Dict[str, Any], 
    global_visium: str, 
    # New explicit arguments
    explicit_snp_vcf: str,
    explicit_panel: str,
    explicit_eagle_dir: str
) -> Dict[str, Any]:
    """Generates the dictionary for Step 1 based on args + main config."""
    cal = pipeline_cfg.get("calicost", {})
    
    # Base paths
    calicost_dir = cal.get("calicost_dir", "")
    
    # Use explicit arguments passed from Snakemake (Common Refs)
    region_vcf    = explicit_snp_vcf
    phasing_panel = explicit_panel
    eagledir      = explicit_eagle_dir
    
    # Fallback to Spaceranger path if not specific
    spaceranger = cal.get("spaceranger_dir", global_visium)

    # Check required Step 1 parameters
    if not all([calicost_dir, eagledir, region_vcf, phasing_panel]):
        raise ValueError("[Step1] Missing required paths (calicost_dir, eagle_dir, region_vcf, or phasing_panel). Check config or args.")

    step1 = {
        "calicost_dir": calicost_dir,
        "eagledir": eagerly_norm(eagledir),
        "samtools_sorting_mem": cal.get("samtools_sorting_mem", STEP1_DEFAULTS["samtools_sorting_mem"]),
        "UMItag": cal.get("UMItag", STEP1_DEFAULTS["UMItag"]),
        "cellTAG": cal.get("cellTAG", STEP1_DEFAULTS["cellTAG"]),
        "nthreads_cellsnplite": int(cal.get("nthreads_cellsnplite", STEP1_DEFAULTS["nthreads_cellsnplite"])),
        "region_vcf": region_vcf,
        "phasing_panel": phasing_panel,
        "chromosomes": cal.get("chromosomes", STEP1_DEFAULTS["chromosomes"]),
        "spaceranger_dir": spaceranger,
        # output_snpinfo will be set/overridden in main()
        "output_snpinfo": "", 
    }
    return step1

# --- Main CLI ---
def main():
    ap = argparse.ArgumentParser(description="Prepare CalicoST configs using a Template for Step 2.")
    
    # Inputs
    ap.add_argument("--pipeline-config", required=True, help="Main pipeline config.yaml")
    ap.add_argument("--template-file", required=True, help="Path to the physical Step 2 YAML template")
    ap.add_argument("--visium-dir", required=True, help="Path to Spaceranger output")
    
    # New Arguments for Common Refs (matching Snakemake params)
    ap.add_argument("--snp-vcf", required=True, help="Path to SNP VCF")
    ap.add_argument("--phasing-panel", required=True, help="Path to Phasing Panel Directory")
    ap.add_argument("--eagle-dir", required=True, help="Path to Eagle Directory")

    # Outputs
    ap.add_argument("--out-step1", default="calicost_step1.yaml")
    ap.add_argument("--out-step2", default="configuration_cna.yaml")
    
    # Dynamic Paths for Step 2
    ap.add_argument("--default-snp-dir", required=True, help="Where SNP info will be stored")
    ap.add_argument("--default-output-dir", required=True, help="Where CNA output will be stored")
    
    # The crucial fix for Normal Index
    ap.add_argument("--normalidx", default="", help="Path to normal ref barcodes. If empty, Template default is used.")
    
    ap.add_argument("--create-outdirs", action="store_true", help="Create output directories automatically")

    args = ap.parse_args()

    # 1. Read Main Pipeline Config
    pipeline_cfg = read_yaml(args.pipeline_config)
    
    # Extract global CalicoST directory (needed for template replacement)
    calicost_dir = pipeline_cfg.get("calicost", {}).get("calicost_dir", "")
    if not calicost_dir:
        raise ValueError("Critical: 'calicost_dir' is not defined in the main pipeline config.")

    # ---------------------------
    # Generate Step 1 Config
    # ---------------------------
    step1_data = generate_step1_config(
        pipeline_cfg, 
        args.visium_dir,
        explicit_snp_vcf=args.snp_vcf,
        explicit_panel=args.phasing_panel,
        explicit_eagle_dir=args.eagle_dir
    )
    
    # Ensure snp output matches the dynamic path passed from Snakemake
    step1_data["output_snpinfo"] = args.default_snp_dir
        
    write_yaml(step1_data, args.out_step1, 
               header_comment="CalicoST Step1 Config (Auto-generated)", 
               flow_keys={"chromosomes"})
    print(f"[{now()}] Step 1 config written to: {args.out_step1}")

    # ---------------------------
    # Generate Step 2 Config (Template Strategy)
    # ---------------------------
    if not os.path.exists(args.template_file):
        raise FileNotFoundError(f"Step 2 Template file not found: {args.template_file}")

    # A. Read Template as Raw Text
    with open(args.template_file, "r") as f:
        raw_text = f.read()

    # B. Replace Placeholder with Real Path
    # We replace __CALICOST_DIR__ with the actual path from config.yaml
    filled_text = raw_text.replace("__CALICOST_DIR__", calicost_dir)

    # C. Parse into Dictionary
    step2_data = yaml.safe_load(filled_text)

    # D. Update Dynamic Run Paths
    step2_data["spaceranger_dir"] = args.visium_dir
    step2_data["snp_dir"]         = args.default_snp_dir
    step2_data["output_dir"]      = args.default_output_dir

    # E. Handle Normal Index (Ref vs No-Ref Logic)
    if args.normalidx and args.normalidx.strip() != "":
        step2_data["normalidx_file"] = args.normalidx

    write_yaml(step2_data, args.out_step2, 
               header_comment="CalicoST Step2 Config (From Template)")
    print(f"[{now()}] Step 2 config written to: {args.out_step2}")

    # ---------------------------
    # Create Directories
    # ---------------------------
    if args.create_outdirs:
        ensure_dir(args.default_snp_dir)
        ensure_dir(args.default_output_dir)

if __name__ == "__main__":
    main()