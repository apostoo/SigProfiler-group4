import pandas as pd
import umap
import hdbscan
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import csv
from scipy.stats import f_oneway
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from pathlib import Path


# Fixed parameters
UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.7
HDBSCAN_MIN_CLUSTER_SIZE = 15
NMF_COMPONENTS = 2

def load_exposure_matrix(sig_type, sig_no):
    input_dir = f"../results/exposure-matrices/results-{sig_no}-signatures/{sig_type}/Suggested_Solution/{sig_type}_De-Novo_Solution/Activities/{sig_type}_De-Novo_Activities_refit.txt"
    exposures_df = pd.read_csv(input_dir, sep='\t', index_col=0)
    norm_exposures_df = exposures_df.div(exposures_df.sum(axis=1), axis=0)
    return norm_exposures_df

def embed_and_cluster_UMAP(norm_exposure):
    umap_model = umap.UMAP(n_neighbors=UMAP_N_NEIGHBORS, min_dist=UMAP_MIN_DIST, n_components=2)
    embedding = umap_model.fit_transform(norm_exposure)
    
    clusterer = hdbscan.HDBSCAN(min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE)
    labels = clusterer.fit_predict(embedding)
    
    return embedding, labels

def embed_and_cluster_NMF(norm_exposure):
    model = NMF(n_components=NMF_COMPONENTS, init='nndsvda', max_iter=500)
    embedding = model.fit_transform(norm_exposure)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE)
    labels = clusterer.fit_predict(embedding)
    return embedding, labels


def plot_embedding(embedding, labels, sig_type, sig_no, method):
    mask = labels != -1
    embedding_plot = embedding[mask]
    labels_plot = labels[mask] + 1  # Label clusters starting from 1

    plot_df = pd.DataFrame({
        f'{method}_1': embedding_plot[:, 0],
        f'{method}_2': embedding_plot[:, 1],
        'Cluster': labels_plot
    })

    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=plot_df, x=f'{method}_1', y=f'{method}_2', hue='Cluster', palette="tab20", legend="full")
    plt.title(f"{method} + HDBSCAN Clustering (Sigs={sig_no})")
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    output_dir = Path(f"../results/clustering/clustering-{sig_no}-signatures/fixed_params")
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / f"{method}_embedding_{sig_type}_{sig_no}_sigs.jpg")
    plt.close()

    noise_pct = (labels == -1).sum() / len(labels) * 100
    noise_csv_dir = output_dir / f"{method}_HDBSCAN_noise_summary.csv"

    return noise_pct, noise_csv_dir

def save_cluster_profile_summary(df, catalog, signatures_no, method):
    output_path = Path(f"../results/summary_cluster_profiles_across_signature_numbers_{method}.csv")
    records = []

    for cluster_id in sorted(df['cluster'].unique()):
        cluster_data = df[df['cluster'] == cluster_id].drop(columns='cluster')
        mean_profile = cluster_data.mean()
        row = {
            "catalog": catalog,
            "signature_no": signatures_no,
            "cluster_id": cluster_id,
        }
        row.update(mean_profile.to_dict())
        records.append(row)

    profile_df = pd.DataFrame(records)
    profile_df.to_csv(output_path, mode='a', header=not output_path.exists(), index=False)

def plot_cluster_signature_averages(norm_exposure, labels, mutational_catalog, signatures_no, method):
    df = norm_exposure.copy()
    df['cluster'] = labels
    df = df[df['cluster'] != -1]
    df['cluster'] += 1  # Make cluster IDs 1-based

    output_dir = Path(f"../results/clustering/clustering-{signatures_no}-signatures/cluster_profiles_{method}")
    output_dir.mkdir(parents=True, exist_ok=True)

    for cluster_id in sorted(df['cluster'].unique()):
        cluster_data = df[df['cluster'] == cluster_id].drop(columns='cluster')
        mean_profile = cluster_data.mean()

        plt.figure(figsize=(10, 4))
        sns.barplot(x=mean_profile.index, y=mean_profile.values, palette="Blues_d", hue=mean_profile.index, legend=False)
        plt.xticks(rotation=90)
        plt.ylabel("Mean Exposure")
        plt.ylim(0.0, 1.0)
        plt.title(f"{method} Cluster {cluster_id} Signature Profile")
        plt.tight_layout()
        plt.savefig(output_dir / f"{method}_cluster_{cluster_id}_signature_profile_{mutational_catalog}_{signatures_no}_sigs.jpg")
        plt.close()

    cluster_means = df.groupby("cluster").mean()
    signature_variance = cluster_means.var(axis=0)

    plt.figure(figsize=(10, 4))
    sns.barplot(x=signature_variance.index, y=signature_variance.values, palette="Reds_d", hue = signature_variance.index, legend=False)
    plt.xticks(rotation=90)
    plt.ylabel("Variance Across Clusters")
    plt.title(f"{method}: Signature Variance Across Clusters ({mutational_catalog}, {signatures_no} sigs)")
    plt.tight_layout()
    plt.savefig(output_dir / f"{method}_signature_variance_across_clusters_{mutational_catalog}_{signatures_no}_sigs.jpg")
    plt.close()

    signature_variance.to_csv(output_dir / f"{method}_signature_variance_{mutational_catalog}_{signatures_no}_sigs.csv", header=["variance"])

    anova_results = {}
    for signature in norm_exposure.columns:
        groups = [df[df['cluster'] == cl][signature].values for cl in df['cluster'].unique()]
        stat, pval = f_oneway(*groups)
        anova_results[signature] = {
            "p_value": pval,
            "significant": pval < 0.05
        }

    anova_df = pd.DataFrame.from_dict(anova_results, orient='index')
    anova_df.index.name = "signature"
    anova_df.sort_values('p_value', inplace=True)
    anova_df.to_csv(output_dir / f"{method}_anova_signature_pvalues_{mutational_catalog}_{signatures_no}_sigs.csv")

    save_cluster_profile_summary(df, mutational_catalog, signatures_no, method)


def evaluate_clustering_metrics(norm_exposure, labels, sig_type, sig_no, method):
    output_dir = Path(f"../results/clustering/clustering-{sig_no}-signatures/fixed_params")
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_csv = output_dir / f"{method}_clustering_metrics_summary.csv"

    # Mask out noise points (-1)
    mask = labels != -1
    filtered_data = norm_exposure[mask]
    filtered_labels = labels[mask]

    if len(set(filtered_labels)) > 1:
        sil_score = silhouette_score(filtered_data, filtered_labels)
        dbi_score = davies_bouldin_score(filtered_data, filtered_labels)
        ch_score = calinski_harabasz_score(filtered_data, filtered_labels)
    else:
        sil_score, dbi_score, ch_score = "NA", "NA", "NA"

    with open(metrics_csv, 'a', newline='') as f:
        writer = csv.writer(f)
        if f.tell() == 0:
            writer.writerow(['Catalog', 'Signatures', 'Method', 'Silhouette', 'DaviesBouldin', 'CalinskiHarabasz'])
        writer.writerow([sig_type, sig_no, method, sil_score, dbi_score, ch_score])


def run_pipeline(method):
    if method not in ['UMAP', 'NMF']:
        raise ValueError("Method must be 'UMAP' or 'NMF'")

    embed_and_cluster = embed_and_cluster_UMAP if method == 'UMAP' else embed_and_cluster_NMF

    for mutational_sig_type in ["SBS96"]:
        for sig_no in range(6, 8):
            norm_exposure = load_exposure_matrix(mutational_sig_type, sig_no)
            embedding, labels = embed_and_cluster(norm_exposure)

            noise_pct, noise_csv = plot_embedding(embedding, labels, mutational_sig_type, sig_no, method)
            evaluate_clustering_metrics(norm_exposure, labels, mutational_sig_type, sig_no, method)

            plot_cluster_signature_averages(norm_exposure, labels, mutational_sig_type, sig_no, method)


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=FutureWarning)
    run_pipeline("UMAP")
    run_pipeline("NMF")