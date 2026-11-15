# benchmark_spatialCNV

Analysis and figure code for a benchmark of copy-number-variant (CNV) callers on
spatial transcriptomics data — Task 1 (CNV accuracy vs DNA/WES), Task 2 (subclone vs
spatial structure), Task 3 (tumour vs histology), and high-resolution ST (VisiumHD, Atera).

## Layout

```
figures/     final figure scripts + the frozen inputs they read
  envs/      figures-r4.4.yml — pinned environment for the figure scripts
  utils/     shared setup + palette (figure_setup.R, figure_settings.tsv, color.R)
  task1_cnv_accuracy/ task2_subclone/ task3_tumor_identification/ high_res_st/
    input_data/   frozen intermediates (CSV/RDS/parquet) — read directly by the scripts
    scripts/      one script per figure
realdata/
  pipeline/  Snakemake CNV pipeline + conda env specs (envs/)
  analysis/  task1/2/3 standardize + evaluate scripts + dataset registry
  highres/   high-resolution ST analysis
    visiumHD/  one folder per pipeline stage, in run order:
               binning/ segmentation/ rctd/ pseudo_spot/ infercnv/
    atera/     atera processing + inferCNV run
synthetic/
  pipeline/  run_synthetic_benchmark.sh (end-to-end driver) + master_pipeline.py
             + tool runners (run_<tool>) + tool-input prep (prep_<tool>)
  analysis/  extract / evaluate / subclone steps + their drivers
```

The synthetic-data **simulator** is released separately as the `stcnv-sim` package
(this repo contains only the analysis run on its outputs).

## Reproducing the figures

```bash
conda env create -f figures/envs/figures-r4.4.yml     # creates `stcnv-figures`
conda activate stcnv-figures

cd figures/task1_cnv_accuracy/scripts && Rscript fig2_task1_main.R
```

Each script self-locates, reads only `../input_data/`, and writes to `../main_figures/`,
`../supp_figures/` and `../data_tables/` (one CSV per panel holding the exact plotted
values). No configuration is needed.

## CNV caller versions

The versions benchmarked in the paper:

| Caller | Version | Caller | Version |
|---|---|---|---|
| inferCNV | 1.22.0 | Numbat | 1.5.1 |
| infercnvpy | 0.6.1 | XClone | 0.4.0 |
| CopyKAT | 1.1.0 | CalicoST | 1.0.0 |
| SCEVAN | 1.0.3 | copyVAE | 1.0 |
| CONICSmat | 0.0.0.1 | TUSCAN | 1.0.0 |
| STARCH | no released version — local git checkout at commit `b9c50f4` | | |

Per-tool conda environments are in `realdata/pipeline/envs/`.

### Path placeholders

Scripts reference external data through placeholders that must be substituted before
running. References to other files *inside* this repository are relative and need no
configuration.

| Placeholder | Meaning |
|---|---|
| `<DATA_ROOT>` | Public datasets, per-sample processed data, results, high-res ST data, synthetic run outputs |
| `<REF_ROOT>` | Reference annotation: gencode gene positions, Numbat/XClone reference panels, pipeline `refs/` |
| `<TOOLS_ROOT>` | External tool installs: Eagle, CalicoST, STARCH, copyVAE |
| `<CONDA_ROOT>` | Conda/miniconda installation (interpreter + envs) |

```bash
grep -rl '<DATA_ROOT>' . | xargs sed -i 's|<DATA_ROOT>|/your/data/root|g'
```

## License

MIT — see [LICENSE](LICENSE).
