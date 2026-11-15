# NovaSeq X Illumina spatial — breast cancer (high-res ST platform #4)

Fourth technology in the `high_res_ST` arm, alongside `atera` (Xenium), `cosmx`, `visiumHD`.

## Aim

Run inferCNV on the Illumina NovaSeq X spatial transcriptomics platform and manually
call subclones, so the platform can be compared against the other high-res ST
technologies in the SpatialCNV benchmark. Analysis mirrors `high_res_ST/atera`.

## Samples (3)

Raw: `<NOVASEQ_RAW>`

| Sample | BaseSpace id | barcodes | segmented cells | tissue extent (um) |
|---|---|---|---|---|
| `DCIS_IDC_Grade1` | 924422509 | 20,377,052 | 169,018 | 6161 x 12722 |
| `IDC_Grade2` | 924252345 | 27,037,184 | 388,581 | 7158 x 9926 |
| `IDC_Grade3` | 924485569 | 44,796,375 | 718,405 | 10252 x 9512 |

Per sample the download contains, and *only* contains:

- `Raw_Matrix_Files/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}` —
  **gene x spatial-barcode**, 61,906 genes, sub-cellular barcodes (~14 counts each).
  **Not** a cell-level matrix.
- `*_registered_nuclei_contour_coords_local.csv` and
  `*_registered_Expanded_5um_cell_contour_coords_local.csv` — segmentation polygons
  (`cell_id, vertex_x, vertex_y`), in micrometres.
- `*_registered.ome.tiff` — registered H&E.

No histology annotation yet — pathologists will annotate at **region level** later.

## Deliverables

1. **Prep + QC** — cell x gene matrix, `nFeature_spatial` / `nCount_spatial`
   distributions and spatial maps. **DONE 2026-09-06** (see Results below).
2. **inferCNV heatmap with subclones** — reference = benign/stromal cells; subclones
   called by cutting the inferCNV dendrogram, checked manually.
   Clustering + compartment annotation done for all three samples (see below).
3. **Subclone spatial map** — subclone label plotted on the cell coordinates.

## Key finding during prep: barcode coordinates are encoded in the barcode id

The download ships **no barcode position file**. Barcode ids look like `SBC:947855:4086433`,
and these two fields are the position in nanometres, **y first**:

```
SBC:<y_nm>:<x_nm>   ->  x_um = f3/1000,  y_um = f2/1000
```

in the same frame as the segmentation CSVs. This was inferred from range matching, not
from vendor documentation, so it is verified explicitly in
`code/novaseq_00_check_coords.py`: 100% of 20-um bins containing a segmented cell
centroid also contain barcodes, in all three samples, and the coordinate ranges agree
to within ~50 um. Barcodes cover a slightly larger area than the tissue (off-tissue
capture area), which is why the Jaccard is < 1.

## Pipeline decisions

- **Cell definition:** `Expanded_5um` cell contours, not nuclei-only — matches the
  Xenium expanded-cell convention used in `atera`, and gives more counts per cell for
  inferCNV. Nuclei contours are kept in `code/novaseq_samples.py` as an alternative.
- **Barcode -> cell assignment:** polygons rasterised into a label mask at 0.5 um/px,
  then each barcode looked up in the mask (`code/novaseq_01_prep_cells.py`).
  Overlapping expanded contours are resolved last-writer-wins.
- **Environment:** `bin2cell_env` for the heavy prep (numpy/cv2/anndata),
  `stcnv-r-shared` (Seurat 5) for anything that has to be an R object.

## Results — deliverable 1 (2026-09-06)

`processed_data/cells/<sample>_cells.h5ad`, obs = `x, y, area_um2, n_barcodes,
nCount_spatial, nFeature_spatial, sample`; var = 61,906 genes (ensembl id + symbol).

| Sample | cells | barcodes in a cell | median nCount | median nFeature | median area | counts in cells |
|---|---|---|---|---|---|---|
| DCIS_IDC_Grade1 | 169,018 | 41.5% | 1,427 | 731 | 53.8 um2 | 0.39e9 |
| IDC_Grade2 | 388,581 | 71.9% | 2,197 | 1,316 | 52.0 um2 | 1.01e9 |
| IDC_Grade3 | 718,405 | 76.6% | 3,155 | 1,588 | 46.8 um2 | 2.72e9 |

Depth is far above Xenium/CosMx — 1.4-3.2k counts and 0.7-1.6k genes per cell on a
whole-transcriptome panel — so single-cell inferCNV is realistic without pseudobulking.
At `nCount >= 50 & nFeature >= 25`, 97.0 / 99.6 / 99.6% of cells pass.

Figures in `processed_data/prep_qc/`:
`qc00_coord_decoding_check.{png,csv}`, `qc01_summary.csv`,
`qc01_distributions.{png,pdf}`, `qc01_violin.{png,pdf}`, `qc01_spatial.{png,pdf}`.

`qc01_spatial.png` is the second, independent confirmation of the coordinate decoding:
tissue morphology comes out clean — bright DCIS ducts on dim stroma in Grade1, tumour
nests vs stromal bands in Grade2/3.

Grade1's lower in-cell barcode fraction (41.5%) is the fatty/stromal, sparsely cellular
tissue, not a failure — its tissue area is largely acellular.

## Results — clustering + compartments (2026-09-07)

Seurat **v5.1.0 / SeuratObject v5.4.0** (`stcnv-r-shared`). Note there is no Seurat
5.4.0 release — 5.4.0 is the SeuratObject version.

Pipeline: LogNormalize -> 2000 vst variable features -> scale -> PCA (30 PCs) ->
UMAP (30 PCs) -> SNN graph (30 PCs) -> Louvain at resolutions 0.1-1.0.
Cells pre-filtered at `nCount >= 50 & nFeature >= 25`, genes at >= 10 cells.

Compartments are assigned by z-scoring each marker's per-cluster mean across clusters
and taking the highest-scoring marker set. `margin` = gap to the runner-up compartment;
clusters with `margin < 0.5` are flagged `compartment_ambiguous`.

**`DCIS_IDC_Grade1` — frozen at res 0.1** (Jintong's pick; the marker margins favour
0.8, where all 25 clusters are unambiguous and immune reaches its true ~7%).

| cluster | cells | % | label | margin | note |
|---|---|---|---|---|---|
| 0 | 68,440 | 41.7 | stroma | **0.058** | stromal/immune tie — mixed cluster |
| 1 | 41,597 | 25.4 | epi | 1.421 | |
| 2 | 27,147 | 16.6 | epi | 1.920 | |
| 3 | 17,495 | 10.7 | epi | 1.074 | |
| 4 | 5,638 | 3.4 | stroma | **0.459** | low-content cells |
| 5 | 1,728 | 1.1 | imm | 0.582 | T-cell-infiltrated stroma, focal aggregate |
| 6 | 972 | 0.6 | stroma | 0.828 | VWF/PECAM1 endothelial |
| 7 | 948 | 0.6 | epi | **0.465** | epithelial + stromal markers |

Totals: epi 87,187 (53.2%) / imm 1,728 (1.1%) / stroma 75,050 (45.8%).
45.8% of cells sit in an ambiguous cluster, almost all of it cluster 0.

**Consequence for deliverable 2:** at res 0.1 the immune compartment is 1.1%, versus
~7% at res 0.8 — T cells, plasma cells and mast cells are absorbed into cluster 0.
So cluster 0 is a stroma+immune mixture, and using it wholesale as the inferCNV
reference mixes two cell types with different baseline expression. Either use res 0.8
for the reference definition, or drop the flagged clusters from the reference.

**Useful for subclones:** epithelial clusters 1, 2 and 3 are *spatially segregated* into
distinct duct groups (cluster 3 = the single large duct at left, cluster 2 = the
upper-right duct cluster, cluster 1 = the remaining ducts) — a strong prior for
inferCNV subclone structure even at this coarse resolution.

### Frozen resolutions

| Sample | res | clusters | epi | imm | stroma | ambiguous cells | inferCNV reference (imm+stroma) |
|---|---|---|---|---|---|---|---|
| `DCIS_IDC_Grade1` | 0.1 | 8 | 87,187 (53.2%) | 1,728 (1.1%) | 75,050 (45.8%) | 45.8% | 76,778 (46.8%) |
| `IDC_Grade2` | 0.1 | 8 | 334,815 (86.5%) | 50,741 (13.1%) | 1,645 (0.4%) | 57.7% | 52,386 (13.5%) |
| `IDC_Grade3` | 0.3 | 12 | 525,828 (73.5%) | 74,473 (10.4%) | 115,178 (16.1%) | 38.2% | 189,651 (26.5%) |

**`IDC_Grade2` — res 0.1 caveat.** Only 8 clusters for 387k cells in an invasive
carcinoma: every cluster co-expresses keratins *and* COL1A1/COL1A2, so they are spatial
mixtures rather than cell types. Stroma comes out at 0.4%, which is not biologically
possible for this section (res 0.7-1.0 give 17-22%); real stroma is absorbed into "epi"
cluster 0 (54% of cells, margin 0.257) or labelled "imm" in cluster 2. Two clusters
contain 2 cells each. The reference is still ~52k cells and cluster 2 genuinely is the
non-epithelial compartment, so inferCNV will run — but the tumour side is
stroma-contaminated, which dilutes the CNV signal.

**`IDC_Grade3` — res 0.3.** res 0.1 gives only 5 clusters and 0% stroma (the same
failure mode, more extreme); 0.3 is the lowest resolution where the stromal compartment
appears. Cluster 3 is immune, 4/7/8 stromal, the rest epithelial (cluster 7 overridden
imm -> stroma on 2026-10-04, see `docs/IDC_Grade3_compartment_override.md`). Epithelial clusters
0, 1, 2, 5, 6, 9 occupy *spatially distinct territories* within the tumour — again a
useful subclone prior.

The general failure mode: compartment scores are z-scored across clusters, so when a
resolution yields very few clusters and all of them are mixtures, no compartment can
win on relative expression and the smaller compartments vanish. Margins below 0.5 flag
exactly this, and are recorded per cell in `compartment_ambiguous`.

## Layout

```
code/    novaseq_samples.py           shared sample paths
         novaseq_00_check_coords.py   barcode-coordinate decoding check
         novaseq_01_prep_cells.py     barcode -> cell -> cell x gene h5ad
         novaseq_02_qc_plots.py       deliverable 1 QC panels
processed_data/
         prep_qc/                     QC figures + tables
         cells/                       <sample>_cells.h5ad
```
