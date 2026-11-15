"""Shared sample paths for the NovaSeq X Illumina spatial breast-cancer set (3 samples)."""
import os

RAW = "<NOVASEQ_RAW>"
PROJ = "<DATA_ROOT>/high_res_ST/novaseq"
OUTDIR = os.path.join(PROJ, "processed_data", "prep_qc")

SAMPLES = [
    dict(
        name="DCIS_IDC_Grade1", short="G1",
        matrix=f"{RAW}/NovaSeqX_Human_Breast_DCIS_IDC_Grade1-924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/Raw_Matrix_Files/matrix.mtx.gz",
        barcodes=f"{RAW}/NovaSeqX_Human_Breast_DCIS_IDC_Grade1-924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/Raw_Matrix_Files/barcodes.tsv.gz",
        features=f"{RAW}/NovaSeqX_Human_Breast_DCIS_IDC_Grade1-924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/Raw_Matrix_Files/features.tsv.gz",
        cell_contours=f"{RAW}/924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/NovaSeqX_Human_Breast_DCIS_IDC_Grade1_registered_Expanded_5um_cell_contour_coords_local.csv",
        nuclei_contours=f"{RAW}/924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/NovaSeqX_Human_Breast_DCIS_IDC_Grade1_registered_nuclei_contour_coords_local.csv",
        image=f"{RAW}/NovaSeqX_Human_Breast_DCIS_IDC_Grade1-924422509/ds.51feb8dee0664f9cbad50ca2facb57f4/NovaSeqX_Human_Breast_DCIS_IDC_Grade1_registered.ome.tiff",
    ),
    dict(
        name="IDC_Grade2", short="G2",
        matrix=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade2-924252345/NovaSeqX_Human_Breast_IDC_Grade2-ds.6f12ba5d64644e98a1df38dd16d5d8e6/Raw_Matrix_Files/matrix.mtx.gz",
        barcodes=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade2-924252345/ds.6f12ba5d64644e98a1df38dd16d5d8e6/Raw_Matrix_Files/barcodes.tsv.gz",
        features=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade2-924252345/ds.6f12ba5d64644e98a1df38dd16d5d8e6/Raw_Matrix_Files/features.tsv.gz",
        cell_contours=f"{RAW}/924252345/ds.6f12ba5d64644e98a1df38dd16d5d8e6/NovaSeqX_Human_Breast_IDC_Grade2_registered_Expanded_5um_cell_contour_coords_local.csv",
        nuclei_contours=f"{RAW}/924252345/ds.6f12ba5d64644e98a1df38dd16d5d8e6/NovaSeqX_Human_Breast_IDC_Grade2_registered_nuclei_contour_coords_local.csv",
        image=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade2-924252345/ds.6f12ba5d64644e98a1df38dd16d5d8e6/NovaSeqX_Human_Breast_IDC_Grade2_registered.ome.tiff",
    ),
    dict(
        name="IDC_Grade3", short="G3",
        matrix=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade3-924485569/NovaSeqX_Human_Breast_IDC_Grade3-ds.42951a3bf9ca427dbbb81b66a7f7ba4c/Raw_Matrix_Files/matrix.mtx.gz",
        barcodes=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade3-924485569/ds.42951a3bf9ca427dbbb81b66a7f7ba4c/Raw_Matrix_Files/barcodes.tsv.gz",
        features=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade3-924485569/ds.42951a3bf9ca427dbbb81b66a7f7ba4c/Raw_Matrix_Files/features.tsv.gz",
        cell_contours=f"{RAW}/924485569/ds.42951a3bf9ca427dbbb81b66a7f7ba4c/NovaSeqX_Human_Breast_IDC_Grade3_registered_Expanded_5um_cell_contour_coords_local.csv",
        nuclei_contours=f"{RAW}/924485569/ds.42951a3bf9ca427dbbb81b66a7f7ba4c/NovaSeqX_Human_Breast_IDC_Grade3_registered_nuclei_contour_coords_local.csv",
        image=f"{RAW}/NovaSeqX_Human_Breast_IDC_Grade3-924485569/ds.42951a3bf9ca427dbbb81b66a7f7ba4c/NovaSeqX_Human_Breast_IDC_Grade3_registered.ome.tiff",
    ),
]
