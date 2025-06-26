import os
import pandas as pd
import numpy as np
import SigProfilerExtractor as sig
from umap import UMAP
from sklearn.manifold import trustworthiness
from numpy.linalg import lstsq
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Current fig options: matrix_96 or matrix_78
current_fig = "matrix_96"
# Dataset options: SBS96 or DBS78
dataset = "SBS96"
signature_count = 4
create_screeplot = False

sample_solution_path = "../../test_" + current_fig + "_output/" + dataset + "/Samples.txt"
# sample_matrix = mutation_types x samples
sample_matrix = pd.read_csv(sample_solution_path, sep="\t", index_col=0)

print("\nOriginal Samples Matrix:")
print(sample_matrix.head())

signature_solution_path = "../../test_" + current_fig + "_output/" + dataset + "/All_Solutions"
signature_folder = dataset + "_" + str(signature_count) + "_Signatures"
text_file = dataset + "_S" + str(signature_count) + "_Signatures.txt"
signatures_file = os.path.join(signature_solution_path, signature_folder, "Signatures", text_file)
# signature_matrix = mutation_types x signatures
signature_matrix = pd.read_csv(signatures_file, sep="\t", index_col=0)

print("\nOriginal Signatures Matrix:")
print(signature_matrix.head())

# exposure_matrix = signatures x samples
exposure_matrix, _, _, _ = lstsq(signature_matrix.values, sample_matrix.values, rcond=None)

# Convert exposure_matrix to a DataFrame
exposure_matrix = pd.DataFrame(exposure_matrix, index=signature_matrix.columns, columns=sample_matrix.columns)
exposure_matrix = exposure_matrix.T

print("\nExposure Matrix:")
print(exposure_matrix.head())

# PCA
if create_screeplot:
    pca = PCA(random_state=42)
else:
    pca = PCA(random_state=42, n_components=2)

pca_result = pca.fit_transform(exposure_matrix)

# Create a DataFrame for PCA results
pca_df = pd.DataFrame(data=pca_result, columns=[f'PC{i+1}' for i in range(pca_result.shape[1])], index=exposure_matrix.index)

# Add dominant mutation to the PCA DataFrame
dominant_signature = exposure_matrix.idxmax(axis=1)
pca_df['dominant_signature'] = dominant_signature

print("\nPCA Results:")
print(pca_df)

# UMAP
umap = UMAP(n_components=2, random_state=42)
umap_result = umap.fit_transform(exposure_matrix)

# Create a DataFrame for UMAP results
umap_df = pd.DataFrame(data=umap_result, columns=['UMAP1', 'UMAP2'], index=exposure_matrix.index)

# Add dominant mutation to the UMAP DataFrame
dominant_signature = exposure_matrix.idxmax(axis=1)
umap_df['dominant_signature'] = dominant_signature

print("\nUMAP Results:")
print(umap_df)

# Calculate trustworthiness of the UMAP embedding
trustworthiness_score = trustworthiness(exposure_matrix, umap_result, n_neighbors=5)
print(f"\nTrustworthiness Score: {trustworthiness_score:.4f}")

def plot_pca(pca_df, name):
    plt.figure(figsize=(12, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], s=50)
    for i, txt in enumerate(pca_df.index):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]), fontsize=6, alpha=0.7)
    plt.title(f'PCA of Samples in {name}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq2/pca_samples_{name}.png")
    plt.show()

def plot_pca_exposure(pca_df, name):
    color_map = {
        'SBS96A': 'pink',
        'SBS96B': 'red',
        'SBS96C': 'blue',
        'SBS96D': 'green',
        'DBS78A': 'pink',
        'DBS78B': 'red',
        'DBS78C': 'blue',
        'DBS78D': 'green',
    }
    
    pca_df['Color'] = pca_df['dominant_signature'].map(color_map).fillna('grey')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], c=pca_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(pca_df.index):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'PCA of Exposures Colored by Dominant Signature in {name}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    if dataset == "SBS96":
        legend_sigs = ['SBS96A', 'SBS96B', 'SBS96C', 'SBS96D']
    else:
        legend_sigs = ['DBS78A', 'DBS78B', 'DBS78C', 'DBS78D']

    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=color_map[sig], markersize=10)
               for sig in legend_sigs]
    plt.legend(handles, legend_sigs, title='Dominant Signature', loc='upper right')
    plt.tight_layout()
    plt.savefig(f"figures/rq2/pca_exposures_{name}.png")
    plt.show()

def plot_umap(umap_df, name):
    plt.figure(figsize=(12, 7))
    plt.scatter(umap_df['UMAP1'], umap_df['UMAP2'], s=50)
    for i, txt in enumerate(umap_df.index):
        plt.annotate(txt, (umap_df['UMAP1'].iloc[i], umap_df['UMAP2'].iloc[i]), fontsize=6, alpha=0.7)
    plt.title(f'UMAP of Samples in {name}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq2/umap_samples_{name}.png")
    plt.show()

def plot_umap_exposure(umap_df, name):
    # create a color for each signature
    color_map = {
        'SBS96A': 'pink',
        'SBS96B': 'red',
        'SBS96C': 'blue',
        'SBS96D': 'green',
        'DBS78A': 'pink',
        'DBS78B': 'red',
        'DBS78C': 'blue',
        'DBS78D': 'green',
    }

    umap_df['Color'] = umap_df['dominant_signature'].map(color_map).fillna('grey')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(umap_df['UMAP1'], umap_df['UMAP2'], c=umap_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(umap_df.index):
        plt.annotate(txt, (umap_df['UMAP1'].iloc[i], umap_df['UMAP2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'UMAP of Exposures Colored by Dominant Signature in {name}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(True)
    if dataset == "SBS96":
        legend_sigs = ['SBS96A', 'SBS96B', 'SBS96C', 'SBS96D']
    else:
        legend_sigs = ['DBS78A', 'DBS78B', 'DBS78C', 'DBS78D']
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=color_map[sig], markersize=10)
               for sig in legend_sigs]
    plt.legend(handles, legend_sigs, title='Dominant Signature', loc='upper right')
    plt.tight_layout()
    plt.savefig(f"figures/rq2/umap_exposures_{name}.png")
    plt.show()

def create_scree_plot_eigenvalue(pca, name):
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_, marker='o')
    plt.title(f'Scree Plot (Eigenvalues) of Exposure in {name}')
    plt.xlabel('Principal Component')
    plt.ylabel('Eigenvalue')
    plt.xticks(range(1, len(pca.explained_variance_) + 1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq2/screeplot_eigenvalue_exposure_{name}.png")
    plt.show()

def create_scree_plot_ratio(pca, name):
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_ratio_, marker='o')
    plt.title(f'Scree Plot (Eigenvalue Ratio) of Exposure in {name}')
    plt.xlabel('Principal Component')
    plt.ylabel('Eigenvalue Variance Ratio')
    plt.xticks(range(1, len(pca.explained_variance_) + 1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq2/screeplot_ratio_exposure_{name}.png")
    plt.show()

if (dataset == "SBS96" or dataset == "DBS78") and not create_screeplot:
    plot_pca(pca_df, current_fig)
    plot_pca_exposure(pca_df, current_fig)
    plot_umap(umap_df, current_fig)
    plot_umap_exposure(umap_df, current_fig)

if create_screeplot:
    # print the eigenvalues for each principal component
    print("\nEigenvalues (Variance Explained by Each PC):")
    for i, eigenvalue in enumerate(pca.singular_values_):
        print(f"PC{i+1}: {eigenvalue:.4f}")

    # print the explained variance ratio for each principal component
    print("\nExplained Variance Ratio:")
    for i, variance in enumerate(pca.explained_variance_ratio_):
        print(f"PC{i+1}: {variance:.4f}")

    create_scree_plot_eigenvalue(pca, current_fig)
    create_scree_plot_ratio(pca, current_fig)
