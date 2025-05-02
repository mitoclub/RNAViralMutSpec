import multiprocessing as mp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils import get_random_spectra, get_equilibrium_freqs

from warnings import filterwarnings
filterwarnings("ignore",)

NSAMPLES = 100_000
THREADS = 64

def process_random_ms(rnd_sp):
    _, rnd_eq_aa_freqs = get_equilibrium_freqs(rnd_sp, rate_col='MutSpec', gc=1)
    replic = rnd_sp['replica'].iloc[0]
    rnd_eq_aa_freqs = rnd_eq_aa_freqs.assign(replic=replic)
    return rnd_eq_aa_freqs


random_ms = get_random_spectra(NSAMPLES, 'poisson')
random_ms.to_csv('./data/random_spectra.csv', index=False, float_format='%g')
print('Random spectra saved to ./data/random_spectra.csv')

data = [random_ms.query('replica == @replic') for replic in random_ms.replica.unique()]

print('Processing random spectra...')
with mp.Pool(processes=THREADS) as pool:
    results = pool.map(process_random_ms, data)

random_eq_freqs = pd.concat(results, ignore_index=False)
X = random_eq_freqs.pivot(index=['replic'], columns='aa', values='eq_freq').fillna(0)
X.to_csv('./data/random_eq_aa_freqs.csv', float_format='%g')
print('Random equilibrium frequencies saved to ./data/random_eq_freqs.csv')


# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(X)

# Calculate explained variance ratio
explained_variance = pca.explained_variance_ratio_

# Create a biplot
plt.figure(figsize=(10, 7))
plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, label='Rnd AA eq freqs', s=10, c='blue')

# Add feature vectors (show only the most valuable ones)
threshold = 0.2  # Adjust this threshold to filter vectors
for i, (comp1, comp2) in enumerate(zip(pca.components_[0], pca.components_[1])):
    magnitude = np.sqrt(comp1**2 + comp2**2)
    if magnitude > threshold:
        plt.arrow(0, 0, comp1 , comp2 , color='red', alpha=0.25, head_width=0.01)
        plt.text(comp1 * 1.15, comp2 * 1.15, X.columns[i], alpha=0.25, color='red', ha='center', va='center')

# Update axis labels with explained variance
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f}%)')
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f}%)')
plt.title('PCA Biplot')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('./figures/aa_eq_freqs_from_rnd_ms12_pca.png', dpi=300)
plt.show()
print('PCA plot saved to ./figures/aa_eq_freqs_from_rnd_ms12_pca.png')
print('Done!')
