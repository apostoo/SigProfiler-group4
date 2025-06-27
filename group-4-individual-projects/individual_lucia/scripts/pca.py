from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cross_decomposition import CCA
from sklearn.metrics.pairwise import cosine_distances
import seaborn as sns

   

def pca_all(matrix,labels):

    # 2 components
    pca = PCA(n_components=2)
    matrix_pca = pca.fit_transform(matrix)

    # Extract group type from labels (e.g., 'WES' from 'WES_0')
    groups = np.array([label.split('_')[0] for label in labels])

    # Define colors for each group
    color_map = {'WES': 'blue', 'WGS': 'red', 'WNS': 'orange'}

    # Plot
    plt.figure(figsize=(8, 6))
    for group in np.unique(groups):
        idx = groups == group
        plt.scatter(matrix_pca[idx, 0], matrix_pca[idx, 1], label=group, color=color_map[group])
        for i in np.where(idx)[0]:

            # Manual placement of annotation to reduce overlap
            if labels[i].startswith("WNS_3"):
                offset = (-20, -10)   # below
            elif labels[i].startswith("WGS_2"):
                 offset = (10,0)    
            elif labels[i].startswith("WNS_0"):
                 offset = (-5,10)   
            elif labels[i].startswith("WGS_0"):
                 offset = (-5,-10)    
            else:
                offset = (10, 0)   
            plt.annotate(labels[i],
                         (matrix_pca[i, 0], matrix_pca[i, 1]),
                         textcoords="offset points",
                         xytext=offset,
                         ha='left',
                         fontsize=8)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
   
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
   

    plt.savefig("individual_lucia/plots/pca_all_ann.pdf")




def cosine_dist(matrix1,labels):
    
    distance_m = cosine_distances(matrix1)
    sim_df = pd.DataFrame(distance_m, index=labels, columns=labels)

    # Plot heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(sim_df, cmap='coolwarm', annot=True, xticklabels=True, yticklabels=True)
    
    plt.gca().xaxis.set_ticks_position('top')      
    plt.gca().xaxis.set_label_position('top') 
    plt.tight_layout()
    plt.savefig("individual_lucia/plots/cosine.pdf")
    
   

        
if __name__ == "__main__":

    # WES 
    df = pd.read_csv('individual_lucia\wEg2\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')
    #Ignore axis for SNV
    matrix1 = df.iloc[:, 1:].to_numpy()
    transpose_exome = matrix1.T

    # WSG data
    df = pd.read_csv('individual_lucia\wSg_4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')
    # Ignore axis for SNV
    matrix2 = df.iloc[:, 1:].to_numpy()
    transpose_wsg = matrix2.T


    # WNS data
    df = pd.read_csv('individual_lucia\wNg4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')
    # Ignore axis for SNV
    matrix3 = df.iloc[:, 1:].to_numpy()
    transpose_wng = matrix3.T


    stacked_matrix = np.vstack([transpose_exome, transpose_wsg,transpose_wng])
    labels = np.array(
    [f"WES_{i}" for i in range(transpose_exome.shape[0])] +
    [f"WGS_{i}" for i in range(transpose_wsg.shape[0])] +
    [f"WNS_{i}" for i in range(transpose_wng.shape[0])]
    )

    # pca_all(stacked_matrix,labels)

    # cosine_dist(stacked_matrix,labels)
    
    # WES_0 - WGS_0
    zero_zero=stacked_matrix[0] - stacked_matrix[2]

    # WES_0 - WNS_0
    zero_zero2=stacked_matrix[0] - stacked_matrix[6]
    print(zero_zero2)

    
    

    