import pandas as pd
import umap
import hdbscan
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv

def load_exposure_matrix(sig_type, sig_no):
    input_dir = f"../results/results-{sig_no}-signatures/{sig_type}/Suggested_Solution/{sig_type}_De-Novo_Solution/Activities/{sig_type}_De-Novo_Activities_refit.txt"
    exposures_df = pd.read_csv(input_dir, sep='\t', index_col=0)
    norm_exposures_df = exposures_df.div(exposures_df.sum(axis=1), axis=0)
    return norm_exposures_df

def evaluate_clustering(X, labels):
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    if n_clusters < 2 or len(set(labels)) == 1:
        return -1, float('inf'), -1
    sil = silhouette_score(X, labels)
    db = davies_bouldin_score(X, labels)
    ch = calinski_harabasz_score(X, labels)
    return sil, db, ch

def run_grid_search(norm_exposure, n_neighbors_list, min_dist_list, min_cluster_size_list):
    results = []
    for n_neighbors in n_neighbors_list:
        for min_dist in min_dist_list:
            umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=2, random_state=42)
            embedding = umap_model.fit_transform(norm_exposure)
            for min_cluster_size in min_cluster_size_list:
                clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
                labels = clusterer.fit_predict(embedding)
                sil, db, ch = evaluate_clustering(embedding, labels)
                results.append({
                    'n_neighbors': n_neighbors,
                    'min_dist': min_dist,
                    'min_cluster_size': min_cluster_size,
                    'silhouette': sil,
                    'db': db,
                    'ch': ch,
                    'labels': labels,
                    'embedding': embedding
                })
    return results

def rank_and_select_best(results, mutational_catalog, signatures_no):
    df = pd.DataFrame([{
        'n_neighbors': r['n_neighbors'],
        'min_dist': r['min_dist'],
        'min_cluster_size': r['min_cluster_size'],
        'silhouette': r['silhouette'],
        'db': r['db'],
        'ch': r['ch']
    } for r in results])

    df['sil_rank'] = df['silhouette'].rank(ascending=False)
    df['db_rank'] = df['db'].rank(ascending=True)
    df['ch_rank'] = df['ch'].rank(ascending=False)
    df['mean_rank'] = df[['sil_rank', 'db_rank', 'ch_rank']].mean(axis=1)

    # Save CSV
    output_dir = Path(f"../results/clustering-{signatures_no}-signatures/metrics")
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"UMAP_HDBSCAN_metrics_{mutational_catalog}_{signatures_no}_sigs.csv"
    df.to_csv(csv_path, index=False)

    best_row = df.sort_values(by='mean_rank').iloc[0]
    for r in results:
        if (
            r['n_neighbors'] == best_row['n_neighbors']
            and r['min_dist'] == best_row['min_dist']
            and r['min_cluster_size'] == best_row['min_cluster_size']
        ):
            return r, df

def plot_best_result(result, method, mutational_catalog, signatures_no):
    embedding = result['embedding']
    labels = result['labels']

    # Exclude noise
    mask = labels != -1
    embedding = embedding[mask]
    labels = labels[mask]

    plt.figure(figsize=(10, 7))
    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=labels, palette="tab20", legend="full")
    plt.title(f"{method} best clustering (n_neighbors={result['n_neighbors']}, min_dist={result['min_dist']}, min_cluster_size={result['min_cluster_size']})")
    plt.xlabel(f"{method} 1")
    plt.ylabel(f"{method} 2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    output_dir = Path(f"../results/clustering-{signatures_no}-signatures/best_{method}")
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / f"{method}_for_{mutational_catalog}_{signatures_no}_sigs.jpg")
    plt.close()

    # Save noise percentage
    noise_pct = (result['labels'] == -1).sum() / len(result['labels']) * 100
    noise_csv = output_dir / f"{method}_noise_summary.csv"
    with open(noise_csv, 'a', newline='') as f:
        writer = csv.writer(f)
        if f.tell() == 0:
            writer.writerow(['Catalog', 'Signatures', 'n_neighbors', 'min_dist', 'min_cluster_size', 'noise_percent'])
        writer.writerow([mutational_catalog, signatures_no, result['n_neighbors'], result['min_dist'], result['min_cluster_size'], round(noise_pct, 2)])


def plot_cluster_signature_averages(norm_exposure, labels, mutational_catalog, signatures_no):
    df = norm_exposure.copy()
    df['cluster'] = labels
    df = df[df['cluster'] != -1]  # Exclude noise

    output_dir = Path(f"../results/clustering-{signatures_no}-signatures/cluster_profiles")
    output_dir.mkdir(parents=True, exist_ok=True)

    for cluster_id in sorted(df['cluster'].unique()):
        cluster_data = df[df['cluster'] == cluster_id].drop(columns='cluster')
        mean_profile = cluster_data.mean()

        plt.figure(figsize=(10, 4))
        sns.barplot(x=mean_profile.index, y=mean_profile.values, palette="Blues_d")
        plt.xticks(rotation=90)
        plt.ylabel("Mean Exposure")
        plt.title(f"Cluster {cluster_id} Signature Profile")
        plt.tight_layout()
        plt.savefig(output_dir / f"cluster_{cluster_id}_signature_profile_{mutational_catalog}_{signatures_no}_sigs.jpg")
        plt.close()

def get_best_embedding_clustering_params():
    n_neighbors_list = [10, 15, 30]
    min_dist_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    min_cluster_size_list = [5, 10, 15]

    all_results = []

    signature_range = range(2, 16)
    mutational_sig_type = "SBS96"

    # STEP 1: Collect results from all signature numbers
    exposures_by_sig = {}
    for sig_no in signature_range:
        norm_exposure = load_exposure_matrix(mutational_sig_type, sig_no)
        exposures_by_sig[sig_no] = norm_exposure
        results = run_grid_search(norm_exposure, n_neighbors_list, min_dist_list, min_cluster_size_list)
        for r in results:
            r['signatures_no'] = sig_no  # tag it
            all_results.append(r)

    # STEP 2: Rank across ALL results and select best global parameter set
    df = pd.DataFrame([{
        'signatures_no': r['signatures_no'],
        'n_neighbors': r['n_neighbors'],
        'min_dist': r['min_dist'],
        'min_cluster_size': r['min_cluster_size'],
        'silhouette': r['silhouette'],
        'db': r['db'],
        'ch': r['ch']
    } for r in all_results])

    df['sil_rank'] = df['silhouette'].rank(ascending=False)
    df['db_rank'] = df['db'].rank(ascending=True)
    df['ch_rank'] = df['ch'].rank(ascending=False)
    df['mean_rank'] = df[['sil_rank', 'db_rank', 'ch_rank']].mean(axis=1)

    # Save combined metrics
    output_dir = Path(f"../results/clustering-global-best/metrics")
    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_dir / "UMAP_HDBSCAN_metrics_ALL.csv", index=False)

    best_row = df.sort_values(by='mean_rank').iloc[0]
    best_n = best_row['n_neighbors']
    best_d = best_row['min_dist']
    best_m = best_row['min_cluster_size']
    print(f"Best Global Params: n_neighbors={best_n}, min_dist={best_d}, min_cluster_size={best_m}")

if __name__ == '__main__':
    get_best_embedding_clustering_params()