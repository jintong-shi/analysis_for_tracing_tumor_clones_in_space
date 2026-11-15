#!/usr/bin/env python3
# Per-step tool runner: loops every Rep_* under --root_dir and runs one --step
# (prep, or a single CNV caller). Driven by run_synthetic_benchmark.sh.
import os
import glob
import subprocess
import argparse
import sys
import shutil

# ---- script paths ----
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_PREP_R   = os.path.join(SCRIPT_DIR, "distribute_inputs.R")
SCRIPT_PREP_PY  = os.path.join(SCRIPT_DIR, "prep_copyvae.py")
SCRIPT_INFERCNV = os.path.join(SCRIPT_DIR, "run_infercnv.R")
SCRIPT_NUMBAT   = os.path.join(SCRIPT_DIR, "run_numbat.R")
SCRIPT_COPYVAE  = os.path.join(SCRIPT_DIR, "run_copyvae.py")
SCRIPT_COPYKAT  = os.path.join(SCRIPT_DIR, "run_copykat.R")
SCRIPT_SCEVAN   = os.path.join(SCRIPT_DIR, "run_scevan.R")
SCRIPT_CONICSMAT = os.path.join(SCRIPT_DIR, "run_conicsmat.R")
SCRIPT_INFERCNVPY = os.path.join(SCRIPT_DIR, "run_infercnvpy.py")
STARCH_SCRIPT_DIR = "<TOOLS_ROOT>/STARCH"
SCRIPT_STARCH_PREP= os.path.join(SCRIPT_DIR, "prep_starch.py")
SCRIPT_XCLONE = os.path.join(SCRIPT_DIR, "run_xclone.py")
SCRIPT_CALICOST_CONVERT = os.path.join(SCRIPT_DIR, "prep_calicost.py")
ANALYSIS_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "analysis")
SCRIPT_POSTPROCESS = os.path.join(ANALYSIS_DIR, "extract_and_normalize.R")
SCRIPT_COPYVAE_POST = os.path.join(ANALYSIS_DIR, "postprocess_copyvae.py")
CALICOST_PYTHON = "<CONDA_ROOT>/envs/calicost_env/bin/python"
CALICOST_PARSE = "<TOOLS_ROOT>/CalicoST/src/calicost/parse_input.py"
CALICOST_MAIN = "<TOOLS_ROOT>/CalicoST/src/calicost/calicost_main.py"
DEFAULT_THREADS = "1"
CALICOST_CONFIG_TEMPLATE = """
spaceranger_dir : __SPACERANGER_DIR__
snp_dir : __SNP_DIR__
output_dir : __OUTPUT_DIR__

# supporting files and preprocessing arguments
geneticmap_file : <TOOLS_ROOT>/CalicoST/GRCh38_resources/genetic_map_GRCh38_merged.tab
hgtable_file : <TOOLS_ROOT>/CalicoST/GRCh38_resources/hgTables_hg38_gencode.txt
normalidx_file : __NORMALIDX_FILE__
tumorprop_file : None
supervision_clone_file : None
filtergenelist_file : <TOOLS_ROOT>/CalicoST/GRCh38_resources/ig_gene_list.txt
filterregion_file : <TOOLS_ROOT>/CalicoST/GRCh38_resources/HLA_regions.bed
secondary_min_umi : 200
min_snpumi_perspot : 50
min_percent_expressed_spots : 0.005
bafonly : False

# phase switch probability
nu : 1.0
logphase_shift : -2.0
npart_phasing : 3

# HMRF configurations
n_clones : __N_CLONES__
n_clones_rdr : 2
min_spots_per_clone : __MIN_SPOTS_PER_CLONE__
min_avgumi_per_clone : 10
maxspots_pooling : 19
tumorprop_threshold : 0.5
max_iter_outer : 10
nodepotential : weighted_sum
initialization_method : rectangle
num_hmrf_initialization_start : 0
num_hmrf_initialization_end : 1
spatial_weight : 1.0
construct_adjacency_method : hexagon
construct_adjacency_w : 1.0

# HMM configurations
n_states : 3
params : smp
t : 1-1e-4
t_phaseing : 0.9999
fix_NB_dispersion : False
shared_NB_dispersion : True
fix_BB_dispersion : False
shared_BB_dispersion : True
max_iter : 10
tol : 0.0001
gmm_random_state : 0
np_threshold : 1.0
np_eventminlen : 10

# integer copy number
nonbalance_bafdist : 1.0
nondiploid_rdrdist : 10.0
""".strip()

CONICSMAT_CHRARM = "<REF_ROOT>/chromosome_arm_positions_grch38.txt"
CONICSMAT_GENE_POS = "<REF_ROOT>/gene_order.csv"

# InferCNV specific
GENE_ORDER_FILE = "<DATA_ROOT>/tests/eval_B/gencode_v38_gene_pos.txt"
# =================================================

def run_cmd(cmd):
    print(f"   >>> Executing: {' '.join(cmd)}")
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
        "RCPP_PARALLEL_NUM_THREADS": "1",
        "JOBLIB_TEMP_FOLDER": "/tmp",
        "TMPDIR": "/tmp",
        "TEMP": "/tmp",
        "TMP": "/tmp",
        "MPLCONFIGDIR": "/tmp/matplotlib-cnv-pipeline",
    })
    subprocess.run(cmd, check=True, env=env)

def write_calicost_config(out_path, replacements):
    config_text = CALICOST_CONFIG_TEMPLATE
    for key, value in replacements.items():
        config_text = config_text.replace(key, value)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(config_text)
        f.write("\n")

def safe_symlink(src, dst):
    src = os.path.abspath(src)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.lexists(dst):
        if os.path.islink(dst):
            try:
                if os.readlink(dst) == src:
                    return
            except OSError:
                pass
        os.remove(dst)
    os.symlink(src, dst)

def build_calicost_tissue_positions(meta_file, out_file):
    try:
        import pandas as pd
    except Exception as exc:
        print(f"    [Warning] pandas unavailable for CalicoST tissue_positions fallback: {exc}")
        return False
    if not os.path.exists(meta_file):
        print(f"    [Warning] Missing spatial metadata: {meta_file}")
        return False
    try:
        df = pd.read_parquet(meta_file)
        if "barcode" in df.columns:
            barcode = df["barcode"].astype(str)
        elif "index" in df.columns:
            barcode = df["index"].astype(str)
        else:
            barcode = df.iloc[:, 0].astype(str)
        if "x" in df.columns:
            x = df["x"]
        elif "x_coord" in df.columns:
            x = df["x_coord"]
        else:
            x = df.iloc[:, 1]
        if "y" in df.columns:
            y = df["y"]
        elif "y_coord" in df.columns:
            y = df["y_coord"]
        else:
            y = df.iloc[:, 2]
        out = pd.DataFrame({
            "barcode": barcode,
            "in_tissue": 1,
            "x": x,
            "y": y,
            "pixel_row": y,
            "pixel_col": x,
        })
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        out.to_csv(out_file, header=False, index=False)
        return True
    except Exception as exc:
        print(f"    [Warning] Failed to build CalicoST tissue_positions_list.csv: {exc}")
        return False

def ensure_calicost_inputs(rep, std_in):
    calicost_in = os.path.join(std_in, "CalicoST")
    snpinfo_dir = os.path.join(calicost_in, "snpinfo")
    os.makedirs(snpinfo_dir, exist_ok=True)

    baf_dir = os.path.join(calicost_in, "baf")
    if not os.path.isdir(baf_dir):
        baf_dir = os.path.join(rep, "baf")

    if os.path.exists(baf_dir):
        req = [
            os.path.join(snpinfo_dir, "cell_snp_Aallele.npz"),
            os.path.join(snpinfo_dir, "cell_snp_Ballele.npz"),
            os.path.join(snpinfo_dir, "unique_snp_ids.npy"),
            os.path.join(snpinfo_dir, "barcodes.txt"),
        ]
        if not all(os.path.exists(p) for p in req):
            run_cmd([
                CALICOST_PYTHON, SCRIPT_CALICOST_CONVERT,
                "--baf_dir", baf_dir,
                "--output_dir", snpinfo_dir,
                "--snps_no_header",
            ])
    else:
        print(f"    [Skip] BAF directory missing: {baf_dir}")

    h5ad_dst = os.path.join(snpinfo_dir, "filtered_feature_bc_matrix.h5ad")
    h5ad_src = os.path.join(std_in, "CopyVAE", "st_expression.h5ad")
    if os.path.exists(h5ad_src):
        safe_symlink(h5ad_src, h5ad_dst)
    else:
        print(f"    [Warning] Missing expression h5ad: {h5ad_src}")

    tissue_src = os.path.join(calicost_in, "spatial", "tissue_positions_list.csv")
    tissue_dst = os.path.join(snpinfo_dir, "spatial", "tissue_positions_list.csv")
    if os.path.exists(tissue_src):
        safe_symlink(tissue_src, tissue_dst)
    else:
        sim_parquet = os.path.join(rep, "synthetic", "spatial_data", "st-meta.parquet")
        build_calicost_tissue_positions(sim_parquet, tissue_dst)

    normalidx_file = os.path.join(calicost_in, "normalidx_file.txt")
    if not os.path.exists(normalidx_file) or os.path.getsize(normalidx_file) == 0:
        normalidx_file = None

    required = [
        os.path.join(snpinfo_dir, "cell_snp_Aallele.npz"),
        os.path.join(snpinfo_dir, "cell_snp_Ballele.npz"),
        os.path.join(snpinfo_dir, "unique_snp_ids.npy"),
        os.path.join(snpinfo_dir, "barcodes.txt"),
        h5ad_dst,
        tissue_dst,
    ]
    missing = [p for p in required if not os.path.exists(p)]
    return snpinfo_dir, normalidx_file, missing

def read_n_tumor_spots(std_in, tumor_prop_cutoff=0.5):
    import csv
    meta_file = os.path.join(std_in, "Metadata_cellprop.csv")
    if not os.path.exists(meta_file):
        return None
    try:
        with open(meta_file, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if "Tumor_prop" not in (reader.fieldnames or []):
                return None
            count = sum(1 for row in reader if float(row["Tumor_prop"]) > tumor_prop_cutoff)
        return count if count > 0 else None
    except Exception:
        return None

def r_prep_outputs_exist(std_in):
    # Outputs produced by distribute_inputs.R only — used to detect
    # the R-script's exit-time SIGSEGV (arrow/libstdc++ teardown crash) when
    # the actual CSV/MTX outputs are already on disk.
    required = [
        os.path.join(std_in, "XClone"),
        os.path.join(std_in, "Metadata_cellprop.csv"),
        os.path.join(std_in, "CalicoST", "spatial", "tissue_positions_list.csv"),
        os.path.join(std_in, "CopyVAE", "counts.mtx"),
    ]
    return all(os.path.exists(p) for p in required)

def prep_outputs_exist(std_in):
    return (
        r_prep_outputs_exist(std_in)
        and os.path.exists(os.path.join(std_in, "CopyVAE", "st_expression.h5ad"))
    )

def main():
    parser = argparse.ArgumentParser(description="Batch Pipeline Controller")
    parser.add_argument("--root_dir", required=True, help="Path to Benchmark Batch folder (containing Rep_*)")
    parser.add_argument("--step", required=True, choices=[
        "prep", "infercnv", "numbat", "copyvae", "copykat", "scevan", "conicsmat",
        "infercnvpy", "starch", "xclone", "calicost", "postprocess"
    ],
                        help="Which tool/step to run")
    parser.add_argument("--xclone_plot", action="store_true", help="Enable XClone plotting (only for --step xclone)")
    parser.add_argument("--HMM", default="TRUE", choices=["TRUE", "FALSE", "true", "false", "1", "0"],
                        help="InferCNV HMM switch (used only for --step infercnv). Default: TRUE")
    parser.add_argument("--infercnv_out_name", default="InferCNV",
                        help="InferCNV output subfolder name under Tool_Outputs (used only for --step infercnv). Default: InferCNV")
    parser.add_argument("--tool_out_name", default="Tool_Outputs",
                        help="Tool output root folder name under each Rep_* (e.g. Tool_Outputs or Tool_Output_TME). Default: Tool_Outputs")
    parser.add_argument("--tool_out_suffix", default="",
                        help="Suffix appended to tool output folder names under Tool_Outputs (e.g. _TME). Default: empty")
    parser.add_argument("--std_in_name", default="Standardized_Inputs",
                        help="Input folder name under each Rep_* (e.g. Standardized_Inputs or Standardized_Inputs_TME). Default: Standardized_Inputs")
    args = parser.parse_args()

    # Find Replicates
    reps = sorted(
        p for p in glob.glob(os.path.join(args.root_dir, "Rep_*"))
        if os.path.isdir(p)
    )
    if not reps:
        print(f"No replicate folders (Rep_*) found in {args.root_dir}")
        sys.exit(1)

    print(f"--- Processing {len(reps)} Replicates | Step: {args.step} ---")

    failed_reps = []

    # Postprocess is a batch-level step, run once and exit
    if args.step == "postprocess":
        # Precompute CopyVAE postprocess per replicate (safe no-op if missing)
        if not os.path.exists(SCRIPT_COPYVAE_POST):
            print(f"[Postprocess] Warning: {SCRIPT_COPYVAE_POST} not found; skipping CopyVAE preprocessing.")
        else:
            for rep in reps:
                rep_name = os.path.basename(rep)
                copyvae_dir = os.path.join(
                    rep,
                    args.tool_out_name,
                    f"CopyVAE{args.tool_out_suffix}",
                )
                copyvae_out = os.path.join(copyvae_dir, "processed_log2_spots_x_genes.csv")
                if os.path.exists(copyvae_out):
                    print(f"[Postprocess] CopyVAE already processed: {rep_name}")
                    continue
                print(f"[Postprocess] CopyVAE preprocessing: {rep_name}")
                try:
                    run_cmd([sys.executable, SCRIPT_COPYVAE_POST, "--tool_dir", copyvae_dir])
                except subprocess.CalledProcessError:
                    print(f"[Postprocess] Warning: CopyVAE preprocess failed for {rep_name}; continuing.")
        post_cmd = [
            "Rscript", SCRIPT_POSTPROCESS,
            "--root_dir", args.root_dir,
            "--tool_out_name", args.tool_out_name,
            "--tool_out_suffix", args.tool_out_suffix,
            "--infercnv_dir_name", args.infercnv_out_name,
        ]
        run_cmd(post_cmd)
        print("\n[Postprocess] Done.")
        return

    for rep in reps:
        rep_name = os.path.basename(rep)
        print(f"\n[{rep_name}]")

        sim_data_dir = os.path.join(rep, "synthetic")
        if not os.path.exists(sim_data_dir):
            print(f"   [Skip] Missing synthetic folder: {sim_data_dir}")
            failed_reps.append(rep_name)
            continue
        
        # Define Standard Paths
        std_in   = os.path.join(rep, args.std_in_name)
        tool_out = os.path.join(rep, args.tool_out_name)
        tool_suffix = args.tool_out_suffix

        def tool_dir(name):
            return os.path.join(tool_out, f"{name}{tool_suffix}")

        try:
            # ---------------- STEP 1: PREP (R + Python) ----------------
            if args.step == "prep":
                # 1. Run R Distribution (Creates RDS, TSV, MTX, XClone Inputs)
                try:
                    run_cmd(["Rscript", SCRIPT_PREP_R, rep])
                except subprocess.CalledProcessError as exc:
                    # Tolerate R's exit-time SIGSEGV (arrow/libstdc++ teardown crash)
                    # when the R-side outputs are actually on disk — step 2 below
                    # then turns the MTX into the h5ad. Only raise if the R outputs
                    # really are missing.
                    if r_prep_outputs_exist(std_in):
                        print(
                            f"   [Warning] Prep R script returned {exc.returncode} "
                            f"but R-side outputs already exist. Continuing to CopyVAE prep."
                        )
                    else:
                        raise
                
                # 2. Run CopyVAE Prep (MTX -> H5AD)
                copyvae_in_dir = os.path.join(std_in, "CopyVAE")
                if os.path.exists(copyvae_in_dir):
                    try:
                        run_cmd(["python3", SCRIPT_PREP_PY, "--input_dir", copyvae_in_dir])
                    except subprocess.CalledProcessError:
                        h5ad = os.path.join(copyvae_in_dir, "st_expression.h5ad")
                        if os.path.exists(h5ad):
                            print("   [Warning] CopyVAE prep failed but st_expression.h5ad exists. Skipping error.")
                        else:
                            raise
                
                # 3. Run STARCH Prep (Raw Parquet -> CSV)
                starch_std_dir = os.path.join(std_in, "STARCH")
                print(f"   -> Running STARCH preparation...")
                # Ensures output goes to Standardized_Inputs/STARCH/gene_expression_matrix.csv
                try:
                    run_cmd([
                        "python3", SCRIPT_STARCH_PREP,
                        "--data_dir", sim_data_dir,
                        "--work_dir", starch_std_dir
                    ])
                except subprocess.CalledProcessError:
                    starch_out = os.path.join(starch_std_dir, "gene_expression_matrix.csv")
                    if os.path.exists(starch_out):
                        print("   [Warning] STARCH prep failed but output exists. Skipping error.")
                    else:
                        raise

            # ---------------- STEP 2: INFERCNV ----------------
            elif args.step == "infercnv":
                in_dir  = os.path.join(std_in, "InferCNV")
                out_dir = tool_dir(args.infercnv_out_name)
                
                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    # Note: run_infercnv.R handles the gene_order formatting internally now
                    run_cmd([
                        "Rscript", SCRIPT_INFERCNV, in_dir, out_dir, GENE_ORDER_FILE,
                        "--HMM", args.HMM
                    ])
                else:
                    print(f"   [Skip] Input dir missing: {in_dir}")

            # ---------------- STEP 3: NUMBAT ----------------
            elif args.step == "numbat":
                in_dir  = os.path.join(std_in, "Numbat")
                out_dir = tool_dir("Numbat")
                
                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    # run_numbat.R expects positional args: <Input_Dir> <Output_Dir>
                    run_cmd(["Rscript", SCRIPT_NUMBAT, in_dir, out_dir])
                else:
                    print(f"   [Skip] Input dir missing: {in_dir}")

            # ---------------- STEP 4: COPYVAE ----------------
            elif args.step == "copyvae":
                in_dir  = os.path.join(std_in, "CopyVAE") # This dir now contains st_expression.h5ad
                out_dir = tool_dir("CopyVAE")
                
                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    run_cmd(["python3", SCRIPT_COPYVAE, "--input_dir", in_dir, "--outdir", out_dir])
                else:
                    print(f"   [Skip] Input dir missing: {in_dir}")
                    
            # ---------------- STEP 5: COPYKAT ----------------
            elif args.step == "copykat":
                in_dir  = os.path.join(std_in, "CopyKAT")
                out_dir = tool_dir("CopyKAT")
                
                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    # Usage: run_copykat.R <Input_Dir> <Output_Dir> [Threads]
                    run_cmd(["Rscript", SCRIPT_COPYKAT, in_dir, out_dir, DEFAULT_THREADS])
                else:
                    print(f"    [Skip] Input dir missing: {in_dir}")

            # ---------------- STEP 6b: SCEVAN ----------------
            elif args.step == "scevan":
                in_dir  = os.path.join(std_in, "SCEVAN")
                out_dir = tool_dir("SCEVAN")

                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    run_cmd(["Rscript", SCRIPT_SCEVAN, in_dir, out_dir, DEFAULT_THREADS])
                else:
                    print(f"    [Skip] Input dir missing: {in_dir}")

            # ---------------- STEP 6c: CONICSMAT ----------------
            elif args.step == "conicsmat":
                in_dir  = os.path.join(std_in, "CONICSmat")
                out_dir = tool_dir("CONICSmat")

                if os.path.exists(in_dir):
                    os.makedirs(out_dir, exist_ok=True)
                    run_cmd(["Rscript", SCRIPT_CONICSMAT, in_dir, out_dir, CONICSMAT_CHRARM, CONICSMAT_GENE_POS])
                else:
                    print(f"    [Skip] Input dir missing: {in_dir}")
                    
            # ---------------- STEP 6: INFERCNVPY ----------------
            elif args.step == "infercnvpy":
                # Inputs
                std_in_dir  = os.path.join(std_in, "InferCNVPy")
                ref_spots   = os.path.join(std_in_dir, "ref_spots.txt")
                
                # The tool script expects the parent of 'spatial_data'
                # Assuming structure: Rep_X/synthetic/spatial_data/...
                custom_data_dir = os.path.join(std_in_dir, "data")
                sim_data_dir = custom_data_dir if os.path.isdir(custom_data_dir) else os.path.join(rep, "synthetic")

                # Output
                out_dir = tool_dir("InferCNVPy")
                
                if os.path.exists(std_in_dir) and os.path.exists(sim_data_dir):
                    # Ensure ref file exists
                    if not os.path.exists(ref_spots):
                        print(f"    [Skip] Reference spots file missing: {ref_spots}")
                        continue
                    
                    # Check if output exists to avoid re-running (optional)
                    if os.path.exists(os.path.join(out_dir, "infercnvpy_gene_cnv.parquet")):
                        print(f"    [Skip] InferCNVPy output already exists in {out_dir}")
                        continue

                    run_cmd([
                        "python3", SCRIPT_INFERCNVPY,
                        "--data_dir", sim_data_dir,
                        "--ref_spots_file", ref_spots,
                        "--gene_order_file", GENE_ORDER_FILE,
                        "--out_dir", out_dir,
                        "--n_jobs", DEFAULT_THREADS,
                        "--out_prefix", "infercnvpy"
                    ])
                else:
                    print(f"    [Skip] Inputs missing (Standardized_Inputs or Synthetic Data dir not found).")

            # ---------------- STEP 7: STARCH (UPDATED) ----------------
            elif args.step == "starch":
                # Input: The CSV created in Step 1 (Standardized_Inputs/STARCH/gene_expression_matrix.csv)
                starch_csv = os.path.join(std_in, "STARCH", "gene_expression_matrix.csv")
                out_dir    = tool_dir("STARCH")
                n_clusters = 4
                print(f"    -> STARCH n_clusters: {n_clusters} (fixed; tool may merge clones freely)")
                
                # Tool Paths
                starch_exec = os.path.join(STARCH_SCRIPT_DIR, "run_STARCH.py")
                starch_ref  = os.path.join(STARCH_SCRIPT_DIR, "hgTables_hg19.txt")

                if os.path.exists(starch_csv):
                    os.makedirs(out_dir, exist_ok=True)
                    print("    -> Running STARCH main tool...")
                    # No mode check, no prep logic here. Just run.
                    run_cmd([
                        "python3", starch_exec,
                        "-c", str(n_clusters),
                        "-t", str(DEFAULT_THREADS),
                        "-i", starch_csv,
                        "--outdir", out_dir,
                        "-m", starch_ref
                    ])
                else:
                    print(f"    [Skip] STARCH input CSV missing: {starch_csv}")
                    print(f"           (Did you run '--step prep' first?)")

            # ---------------- STEP 8: XCLONE ----------------
            elif args.step == "xclone":
                # XClone inputs are derived from the replicate directory
                out_dir = tool_dir("XClone")
                anno_file = os.path.join(std_in, "XClone", "annotation.csv")
                anno_file_alt = os.path.join(std_in, "XClone", "Annotation", "annotation.csv")
                if not os.path.exists(anno_file) and os.path.exists(anno_file_alt):
                    anno_file = anno_file_alt

                missing = [p for p in (rep, anno_file) if not os.path.exists(p)]
                if missing:
                    print(f"    [Skip] XClone inputs missing: {', '.join(missing)}")
                    print(f"           (Did you run '--step prep' first?)")
                    continue

                cmd = [
                    "python3", SCRIPT_XCLONE,
                    "--input_dir", rep,
                    "--output_dir", out_dir,
                    "--std_in_name", args.std_in_name,
                ]
                run_cmd(cmd)
            # ---------------- STEP 10: CALICOST ----------------
            elif args.step == "calicost":
                calicost_in_dir = os.path.join(std_in, "CalicoST")
                out_dir = tool_dir("CalicoST")
                if not os.path.exists(calicost_in_dir):
                    print(f"   [Skip] Input dir missing: {calicost_in_dir}")
                    continue

                os.makedirs(out_dir, exist_ok=True)
                snpinfo_dir, normalidx_file, missing_inputs = ensure_calicost_inputs(rep, std_in)
                if missing_inputs:
                    print("   [Skip] CalicoST inputs missing after preprocessing:")
                    for path in missing_inputs:
                        print(f"          - {path}")
                    continue

                config_path = os.path.join(out_dir, "configfile0")
                n_clones = 4
                print(f"    -> CalicoST n_clones: {n_clones} (fixed; tool may merge clones freely)")
                n_tumor_spots = read_n_tumor_spots(std_in)
                if n_tumor_spots is not None:
                    min_spots_per_clone = max(1, int(0.05 * n_tumor_spots))
                    print(f"    -> CalicoST min_spots_per_clone: {min_spots_per_clone} (5% of {n_tumor_spots} tumor spots)")
                else:
                    min_spots_per_clone = 20
                    print(f"    -> CalicoST min_spots_per_clone: {min_spots_per_clone} (fallback; Metadata_cellprop.csv not found)")
                replacements = {
                    "__SPACERANGER_DIR__": snpinfo_dir,
                    "__SNP_DIR__": snpinfo_dir,
                    "__OUTPUT_DIR__": out_dir,
                    "__NORMALIDX_FILE__": normalidx_file if normalidx_file else "None",
                    "__N_CLONES__": str(n_clones),
                    "__MIN_SPOTS_PER_CLONE__": str(min_spots_per_clone),
                }
                write_calicost_config(config_path, replacements)

                run_cmd([CALICOST_PYTHON, CALICOST_PARSE, "-c", config_path])
                run_cmd([CALICOST_PYTHON, CALICOST_MAIN, "-c", config_path])

        except Exception as exc:
            print(f"   [Error] {rep_name} failed: {exc}")
            failed_reps.append(rep_name)
            continue

    if failed_reps:
        print("\n--- Completed with failures ---")
        for name in failed_reps:
            print(f"   - {name}")
        fail_path = os.path.join(args.root_dir, "failed_reps.txt")
        with open(fail_path, "w", encoding="utf-8") as f:
            for name in failed_reps:
                f.write(f"{name}\n")
        print(f"--- Failed reps saved to: {fail_path} ---")
    else:
        print("\n--- Completed all replicates successfully ---")
if __name__ == "__main__":
    main()
