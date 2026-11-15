import sys
import os
import argparse
import json
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import cell2location
from scipy.io import mmread
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes

def parse_args():
    parser = argparse.ArgumentParser(description="Run Cell2location Pipeline")
    # Reference inputs
    parser.add_argument("--ref_folder", required=True, help="Folder containing counts.mtx, metadata.csv, gene_name.csv")
    parser.add_argument("--cell_type_key", required=True, help="Column name for cell types")
    parser.add_argument("--batch_key", required=True, help="Column name for batch/sample ID")
    
    # Spatial inputs
    parser.add_argument("--visium_path", required=True, help="Path to Visium output folder")
    parser.add_argument("--sample_id", required=True, help="Sample ID for current run")
    
    # Output
    parser.add_argument("--outdir", required=True, help="Output directory")
    
    # Model Parameters
    parser.add_argument("--epochs_ref", type=int, default=250, help="Epochs for reference training")
    parser.add_argument("--epochs_spatial", type=int, default=30000, help="Epochs for spatial mapping")
    parser.add_argument("--use_gpu", action='store_true')
    parser.add_argument("--N_cells", type=float, default=10.0, help="Expected N cells per spot")
    parser.add_argument("--alpha", type=float, default=20.0, help="Detection alpha hyperparameter")
    
    return parser.parse_args()

def load_and_train_reference(ref_folder, batch_key, labels_key, outdir, epochs, use_gpu):
    """
    Loads reference data using manual mtx reading and trains RegressionModel.
    Checks if reference_signatures.csv already exists to avoid re-training.
    """
    ref_sig_file = os.path.join(outdir, "reference_signatures.csv")
    # Check if we already have the output
    if os.path.exists(ref_sig_file):
        print(f"Loading existing reference signatures from {ref_sig_file}")
        return pd.read_csv(ref_sig_file, index_col=0)

    print(f"Processing Reference from {ref_folder}...")
    
    # 1. Load Sparse Matrix
    X = mmread(os.path.join(ref_folder, "counts.mtx"))
    
    # 2. Create AnnData (transpose because mtx is usually Gene x Cell, AnnData needs Cell x Gene)
    # NOTE: Your script implies X.transpose().tocsr(). Check dimensions if error occurs.
    adata_ref = anndata.AnnData(X=X.transpose().tocsr())
    
    # 3. Load Metadata
    cell_meta = pd.read_csv(os.path.join(ref_folder, "metadata.csv"))
    
    # 4. Load Gene Names
    with open(os.path.join(ref_folder, "gene_name.csv"), 'r') as f:
        gene_names = f.read().splitlines()
        
    # 5. Set Observations and Variables
    adata_ref.obs = cell_meta
    
    for col in adata_ref.obs.columns:
        if adata_ref.obs[col].dtype == 'object':
            adata_ref.obs[col] = adata_ref.obs[col].astype(str)

    # Ensure index matches barcode column if it exists, else assume index 0
    if 'barcode' in adata_ref.obs.columns:
        adata_ref.obs.index = adata_ref.obs['barcode']
    
    adata_ref.var.index = gene_names
    adata_ref.var['SYMBOL'] = adata_ref.var.index
    
    # 6. QC and Filtering (As per your script)
    # Using cell2location utils
    selected = filter_genes(adata_ref, cell_count_cutoff=5, 
                            cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_ref = adata_ref[:, selected].copy()
    
    # 7. Setup AnnData
    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key=batch_key,
        labels_key=labels_key
    )
    
    # 8. Train Model
    print("Training Reference Model...")
    mod = RegressionModel(adata_ref)
    
    train_args = {"max_epochs": epochs}
        
    mod.train(**train_args)
    
    # 9. Export Posterior
    adata_ref = mod.export_posterior(
        adata_ref, 
        export_slot='mod',
        add_to_varm={'means','stds','q05','q95'},
        scale_average_detection=True
    )

    mod.save(outdir,overwrite=True)

    # 10. Extract Signatures
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
                                        
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    
    # Save results
    os.makedirs(outdir, exist_ok=True)
    inf_aver.to_csv(ref_sig_file)
    
    # Save the h5ad reference
    adata_ref.write(os.path.join(outdir, "sc_reference_model.h5ad"))
    
    return inf_aver

def _ensure_lowres_image(visium_folder):
    """
    scanpy.read_visium hard-codes spatial/tissue_lowres_image.png. Some spaceranger
    re-runs only emit the hires PNG; create a lowres surrogate by copying hires and
    patch scalefactors_json.json so tissue_lowres_scalef == tissue_hires_scalef.
    """
    spatial_dir = os.path.join(visium_folder, 'spatial')
    lowres = os.path.join(spatial_dir, 'tissue_lowres_image.png')
    hires  = os.path.join(spatial_dir, 'tissue_hires_image.png')
    sf_path = os.path.join(spatial_dir, 'scalefactors_json.json')
    if os.path.exists(lowres) or not os.path.exists(hires):
        return hires if os.path.exists(hires) else None
    import shutil
    shutil.copyfile(hires, lowres)
    if os.path.exists(sf_path):
        with open(sf_path) as f:
            sf = json.load(f)
        if 'tissue_hires_scalef' in sf:
            sf['tissue_lowres_scalef'] = sf['tissue_hires_scalef']
            with open(sf_path, 'w') as f:
                json.dump(sf, f)
    print(f"[c2l] patched missing lowres image: copied {hires} -> {lowres}")
    return hires

def load_visium_robust(visium_folder):
    """
    Robustly loads Visium data handling h5 vs mtx and various tissue position file names.
    """
    # 1. Try reading h5
    if os.path.exists(os.path.join(visium_folder, 'filtered_feature_bc_matrix.h5')):
        hires = _ensure_lowres_image(visium_folder)
        adata_vis = sc.read_visium(visium_folder, source_image_path=hires)
    # 2. Fallback to MTX
    else:
        adata_vis = sc.read_10x_mtx(
            os.path.join(visium_folder, 'filtered_feature_bc_matrix'),
            var_names='gene_symbols', cache=True
        )
        # Load positions manually
        spatial_dir = os.path.join(visium_folder, 'spatial')
        
        # Check for standard csv or list.csv
        pos_file = os.path.join(spatial_dir, 'tissue_positions.csv')
        if not os.path.exists(pos_file):
            pos_file = os.path.join(spatial_dir, 'tissue_positions_list.csv')
            
        # Read positions
        if os.path.exists(pos_file):
            # Check if it has headers (standard tissue_positions.csv usually has header, list does not)
            # Safe bet is to try reading without header if it's the list version
            try:
                positions = pd.read_csv(pos_file)
                if 'in_tissue' not in positions.columns:
                     raise ValueError("Header missing")
            except:
                positions = pd.read_csv(pos_file, header=None)
                positions.columns = [
                    'barcode', 'in_tissue', 'array_row', 'array_col',
                    'pxl_row_in_fullres', 'pxl_col_in_fullres'
                ]
            
            # Match index
            positions = positions.set_index('barcode').loc[adata_vis.obs_names]
            
            # Assign obs
            adata_vis.obs['in_tissue'] = positions['in_tissue'].values
            adata_vis.obs['array_row'] = positions['array_row'].values
            adata_vis.obs['array_col'] = positions['array_col'].values
            
            # Load scale factors
            with open(os.path.join(spatial_dir, 'scalefactors_json.json')) as f:
                scalefactors = json.load(f)
                
            adata_vis.uns['spatial'] = {
                'unknown_sample': { # Dummy key needed for scanpy structure
                    'images': {'hires': None}, # Placeholders
                    'scalefactors': scalefactors
                }
            }
            adata_vis.obsm['spatial'] = positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
        else:
            print("Warning: Spatial position files not found.")
            
    adata_vis.var_names_make_unique()
    return adata_vis

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # Hardware config
    train_kwargs_spatial = {
        "max_epochs": args.epochs_spatial,
        "batch_size": None, # IMPORTANT: From your script
        "train_size": 1     # IMPORTANT: From your script
    }
    
    if args.use_gpu:
        # For spatial mapping (Cell2location), we pass accel/devices
        train_kwargs_spatial["accelerator"] = "gpu"
        train_kwargs_spatial["devices"] = 1
        
    # ==========================
    # 1. Process Reference
    # ==========================
    # We save reference outputs to a 'ref' subdir to keep it clean, or main dir?
    # Your script saves them in `ref_run_name`. Let's put them in `outdir/ref`
    ref_outdir = os.path.join(args.outdir, "reference_signatures")
    
    inf_aver = load_and_train_reference(
        args.ref_folder, 
        args.batch_key, 
        args.cell_type_key, 
        ref_outdir, 
        args.epochs_ref, 
        args.use_gpu
    )

    # ==========================
    # 2. Process Spatial
    # ==========================
    print(f"Processing Spatial Sample {args.sample_id}...")
    
    adata_vis = load_visium_robust(args.visium_path)
    
    # Intersect Genes
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    
    # Assign Sample Key
    adata_vis.obs['sample'] = args.sample_id
    
    # Setup
    Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
    
    # Initialize Model
    print(f"Initializing Model (N_cells={args.N_cells}, alpha={args.alpha})...")
    mod_st = Cell2location(
        adata_vis, 
        cell_state_df=inf_aver,
        N_cells_per_location=args.N_cells,
        detection_alpha=args.alpha
    )
    
    # Train
    print(f"Training Spatial Model for {args.epochs_spatial} epochs...")
    mod_st.train(**train_kwargs_spatial)
    
    # Export Posterior
    adata_vis = mod_st.export_posterior(
        adata_vis, 
        sample_kwargs={'num_samples': 1000, 'batch_size': mod_st.adata.n_obs}
    )
    
    # Save Model
    mod_st.save(os.path.join(args.outdir, f"{args.sample_id}_model"), overwrite=True)
    
    # Save AnnData
    adata_vis.write(os.path.join(args.outdir, "sp.h5ad"))
    
    # Save QC Plot
    mod_st.plot_QC()
    plt.savefig(os.path.join(args.outdir, "ST_mapping_QC.pdf"))
    
    # Extract 5% quantile abundance (Standard Decon result)
    # The key is stored in obsm
    # Extract 5% quantile abundance (Standard Decon result)
    if 'q05_cell_abundance_w_sf' in adata_vis.obsm:
        # 1. Get Column Names
        try:
            cols = adata_vis.uns['mod']['factor_names']
        except KeyError:
            print("Warning: factor_names not found. Generating generic names.")
            cols = [f"Type_{i}" for i in range(adata_vis.obsm['q05_cell_abundance_w_sf'].shape[1])]

        # 2. Get Raw Data (Force Numpy Array)
        # This prevents Pandas from trying to align indices if the obsm is already a DF
        q05_data = adata_vis.obsm['q05_cell_abundance_w_sf']
        if hasattr(q05_data, "values"):
            q05_data = q05_data.values
        else:
            q05_data = np.array(q05_data)

        # 3. Create Clean DataFrame
        df_export = pd.DataFrame(
            q05_data,
            index=adata_vis.obs_names,
            columns=cols
        )
        
        # Save as Decon.csv (Raw abundance)
        df_export.to_csv(os.path.join(args.outdir, "Decon.csv"))
        print(f"Saved Decon.csv to {args.outdir}")
    else:
        print("Error: q05_cell_abundance_w_sf not found.")

if __name__ == "__main__":
    main()