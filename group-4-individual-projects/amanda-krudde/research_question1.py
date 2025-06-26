import pandas as pd
import numpy as np
from umap import UMAP
from sklearn.manifold import trustworthiness
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Current fig options: matrix_96 or matrix_78
current_fig = "matrix_96"
# Dataset options: SBS96 or DBS78
dataset = "SBS96"
create_screeplot = True

solution_path_samples = "../../test_" + current_fig + "_output/" + dataset + "/Samples.txt"
sample_matrix = pd.read_csv(solution_path_samples, sep="\t", index_col=0)

print("\nOriginal Samples Matrix:")
print(sample_matrix.head())

print("\nShape of Samples Matrix:", sample_matrix.shape)

# PCA
if create_screeplot:
    pca = PCA(random_state=42)
else:
    pca = PCA(random_state=42, n_components=2)

pca_result = pca.fit_transform(sample_matrix)

# Create a DataFrame for PCA results
pca_df = pd.DataFrame(data=pca_result, columns=[f'PC{i+1}' for i in range(pca_result.shape[1])], index=sample_matrix.index)

print("\nPCA Results:")
print(pca_df)

# UMAP
umap = UMAP(n_components=2, random_state=42)
umap_result = umap.fit_transform(sample_matrix)

# Create a DataFrame for UMAP results
umap_df = pd.DataFrame(data=umap_result, columns=['UMAP1', 'UMAP2'], index=sample_matrix.index)

print("\nUMAP Results:")
print(umap_df)

# Calculate trustworthiness of the UMAP embedding
trustworthiness_score = trustworthiness(sample_matrix, umap_result, n_neighbors=5)
print(f"\nTrustworthiness Score: {trustworthiness_score:.4f}")

def plot_pca(pca_df, name):
    plt.figure(figsize=(12, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], s=50)
    for i, txt in enumerate(pca_df.index):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]), fontsize=6, alpha=0.7)
    plt.title(f'PCA of Mutation Types in {name}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/pca_mutationtype_{name}.png")
    plt.show()

def plot_pca_coloured(pca_df, name):
    color_map = {
        'C>A': 'blue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'grey',
        'T>C': 'green',
        'T>G': 'pink'
    }
    # Extract mutation class from mutation type (e.g., A[C>A]T → C>A)
    def extract_class(mutation_type):
        return mutation_type[2:5]
    
    # Add mutation class and color to the PCA dataframe
    pca_df['Mutation_Class'] = pca_df.index.map(extract_class)
    pca_df['Color'] = pca_df['Mutation_Class'].map(color_map).fillna('yellow')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], c=pca_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(pca_df.index):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'PCA of Mutation Types Colored by Mutation Class in {name}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/pca_mutationtype_{name}_coloured.png")
    plt.show()

def plot_pca_coloured_double(pca_df, name):
    color_map = {
        'AC': 'lightblue',
        'AT': 'blue',
        'CC': 'lightgreen',
        'CG': 'green',
        'CT': 'pink',
        'GC': 'red',
        'TA': 'lightsalmon',
        'TC': 'orange',
        'TG': 'purple',
        'TT': 'darkblue'
    }

    # Extract mutation class from mutation type (e.g., AC>CA -> AC)
    def extract_class(mutation_type):
        return mutation_type[0:2]

    # Add mutation class and color to the PCA dataframe
    pca_df['Mutation_Class'] = pca_df.index.map(extract_class)
    pca_df['Color'] = pca_df['Mutation_Class'].map(color_map).fillna('grey')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], c=pca_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(pca_df.index):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'PCA of Mutation Types Colored by Mutation Class in {name}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/pca_mutationtype_{name}_coloured.png")
    plt.show()

def plot_umap(umap_df, name):
    plt.figure(figsize=(12, 7))
    plt.scatter(umap_df['UMAP1'], umap_df['UMAP2'], s=50)
    for i, txt in enumerate(umap_df.index):
        plt.annotate(txt, (umap_df['UMAP1'].iloc[i], umap_df['UMAP2'].iloc[i]), fontsize=6, alpha=0.7)
    plt.title(f'UMAP of Mutation Types in {name}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/umap_mutationtype_{name}.png")
    plt.show()

def plot_umap_coloured(umap_df, name):
    color_map = {
        'C>A': 'blue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'grey',
        'T>C': 'green',
        'T>G': 'pink'
    }

    # Extract mutation class from mutation type (e.g., A[C>A]T → C>A)
    def extract_class(mutation_type):
        return mutation_type[2:5]

    # Add mutation class and color to the PCA dataframe
    umap_df['Mutation_Class'] = umap_df.index.map(extract_class)
    umap_df['Color'] = umap_df['Mutation_Class'].map(color_map).fillna('yellow')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(umap_df['UMAP1'], umap_df['UMAP2'], c=umap_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(umap_df.index):
        plt.annotate(txt, (umap_df['UMAP1'].iloc[i], umap_df['UMAP2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'UMAP of Mutation Types Colored by Mutation Class in {name}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/umap_mutationtype_{name}_coloured.png")
    plt.show()

def plot_umap_coloured_double(umap_df, name):
    color_map = {
        'AC': 'lightblue',
        'AT': 'blue',
        'CC': 'lightgreen',
        'CG': 'green',
        'CT': 'pink',
        'GC': 'red',
        'TA': 'lightsalmon',
        'TC': 'orange',
        'TG': 'purple',
        'TT': 'darkblue'
    }

    # Extract mutation class from mutation type (e.g., AC>CA -> AC)
    def extract_class(mutation_type):
        return mutation_type[0:2]

    # Add mutation class and color to the PCA dataframe
    umap_df['Mutation_Class'] = umap_df.index.map(extract_class)
    umap_df['Color'] = umap_df['Mutation_Class'].map(color_map).fillna('grey')

    # Plot with custom colors
    plt.figure(figsize=(12, 7))
    plt.scatter(umap_df['UMAP1'], umap_df['UMAP2'], c=umap_df['Color'], s=50)

    # Optionally label the points (small and light for readability)
    for i, txt in enumerate(umap_df.index):
        plt.annotate(txt, (umap_df['UMAP1'].iloc[i], umap_df['UMAP2'].iloc[i]), fontsize=6, alpha=0.6)

    plt.title(f'UMAP of Mutation Types Colored by Mutation Class in {name}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/umap_mutationtype_{name}_coloured.png")
    plt.show()

def create_scree_plot_eigenvalue(pca, name):
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_, marker='o')
    plt.title(f'Scree Plot (Eigenvalues) of Mutation Types in {name}')
    plt.xlabel('Principal Component')
    plt.ylabel('Eigenvalue')
    plt.xticks(range(1, len(pca.explained_variance_) + 1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/screeplot_eigenvalue_mutationtype_{name}.png")
    plt.show()

def create_scree_plot_ratio(pca, name):
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_ratio_, marker='o')
    plt.title(f'Scree Plot (Explained Variance Ratio) of Mutation Types in {name}')
    plt.xlabel('Principal Component')
    plt.ylabel('Explained Variance Ratio')
    plt.xticks(range(1, len(pca.explained_variance_) + 1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"figures/rq1/screeplot_ratio_mutationtype_{name}.png")
    plt.show()

def plot_pca_3d(pca_df, name):
    color_map = {
        'C>A': 'blue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'grey',
        'T>C': 'green',
        'T>G': 'pink'
    }
    
    # Extract mutation class from mutation type (e.g., A[C>A]T → C>A)
    def extract_class(mutation_type):
        return mutation_type[2:5]
    
    pca_df['Mutation_Class'] = pca_df.index.map(extract_class)
    pca_df['Color'] = pca_df['Mutation_Class'].map(color_map).fillna('yellow')

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pca_df['PC1'], pca_df['PC2'], pca_df['PC3'], c=pca_df['Color'], s=50)

    for i, txt in enumerate(pca_df.index):
        ax.text(pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i], pca_df['PC3'].iloc[i], txt, size=6, alpha=0.7)

    ax.set_title(f'3D PCA of Mutation Types in {name}')
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_zlabel('Principal Component 3')
    plt.tight_layout()
    plt.savefig(f"figures/rq1/pca_mutationtypes_{name}_3d.png")
    plt.show()

if dataset == "SBS96" and not create_screeplot:
    plot_pca(pca_df, current_fig)
    plot_pca_coloured(pca_df, current_fig)
    plot_umap(umap_df, current_fig)
    plot_umap_coloured(umap_df, current_fig)

if dataset == "DBS78" and not create_screeplot:
    # plot_pca(pca_df, current_fig)
    # plot_pca_coloured_double(pca_df, current_fig)
    # plot_umap(umap_df, current_fig)
    plot_umap_coloured_double(umap_df, current_fig)

    # find the two clusters in the UMAP plot
    # and print the mutation types in each cluster
    cluster_1 = umap_df[umap_df['UMAP1'] < -1]
    cluster_2 = umap_df[umap_df['UMAP1'] >= -1]
    print("\nCluster 1 Mutation Types:")
    print(cluster_1.index.tolist())
    print("\nCluster 2 Mutation Types:")
    print(cluster_2.index.tolist())


if create_screeplot:
    # Print the eigenvalues for each principal component
    print("\nEigenvalues (Variance Explained by Each PC):")
    for i, eigenvalue in enumerate(pca.explained_variance_):
        print(f"PC{i+1}: {eigenvalue:.4f}")

    # Print the explained variance ratio for each principal component
    print("\nExplained Variance Ratio:")
    for i, variance in enumerate(pca.explained_variance_ratio_):
        print(f"PC{i+1}: {variance:.4f}")

    # Create scree plots
    create_scree_plot_eigenvalue(pca, current_fig)
    create_scree_plot_ratio(pca, current_fig)

# plot_pca_3d(pca_df, current_fig)
