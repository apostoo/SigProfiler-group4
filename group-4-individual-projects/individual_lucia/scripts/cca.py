import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cross_decomposition import CCA
import seaborn as sns

def cca(wes, wns, wgs):
    # Calculate CCA between the pairs
    group_names = ['WES-WNS', 'WES-WGS', 'WGS-WNS']
    # Store correlations for the bar chart
    first_can = []
    second_can = []

    #WES-WNS
    # 2 components
    cca = CCA(n_components=2)  
    X_c, Y_c = cca.fit_transform(wes, wns)  
    correlations = [np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1] for i in range(X_c.shape[1])]

    first_can.append(correlations[0])
    second_can.append(correlations[1])

    #WES-WGS
    cca = CCA(n_components=2)  
    X_c, Y_c = cca.fit_transform(wes, wgs)  
    correlations = [np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1] for i in range(X_c.shape[1])]

    first_can.append(correlations[0])
    second_can.append(correlations[1])

    #WGS-WNS
    cca = CCA(n_components=2)  
    X_c, Y_c = cca.fit_transform(wgs, wns)  
    correlations = [np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1] for i in range(X_c.shape[1])]

    first_can.append(correlations[0])
    second_can.append(correlations[1])

    
    x = np.arange(len(group_names))  # Group positions
    width = 0.35  # Bar width

    fig, ax = plt.subplots(figsize=(7, 4))

    bars1 = ax.bar(x - width/2, first_can, width, label='1', color='cornflowerblue')
    bars2 = ax.bar(x + width/2, second_can, width, label='2', color='lightcoral')

    # Add text labels above the bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.01,
                    f'{height:.2f}', ha='center', va='bottom', fontsize=9)

    # Formatting
    ax.set_ylabel('Canonical Correlation')
    ax.set_xticks(x)
    ax.set_xticklabels(group_names)
    ax.set_ylim(0, 1)
    ax.legend()

    plt.tight_layout()
    plt.savefig('individual_lucia/plots/cca_corr.pdf')


if __name__ == "__main__":

    # File paths
    wes_file = 'individual_lucia\wEg2\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Activities\SBS96_De-Novo_Activities_refit.txt'
    wns_file = 'individual_lucia\wNg4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Activities\SBS96_De-Novo_Activities_refit.txt'
    wgs_file = 'individual_lucia\wSg_4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Activities\SBS96_De-Novo_Activities_refit.txt'

    wes = pd.read_csv(wes_file, sep='\t', index_col=0).to_numpy()
    wgs = pd.read_csv(wgs_file, sep='\t', index_col=0).to_numpy()
    wns = pd.read_csv(wns_file, sep='\t', index_col=0).to_numpy()


    cca(wes,wns,wgs)


  