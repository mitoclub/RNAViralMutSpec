import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_codons

coda = CodonAnnotation(2)
alphabet = 'ACGT'


def collect_possible_changes():
    nucls = alphabet
    i = 1
    data = []
    for cdn1 in possible_codons:
        aa1 = coda.translate_codon(cdn1)
        for pic in range(3):
            nuc1 = cdn1[pic]
            for nuc2 in nucls:
                if nuc1 == nuc2:
                    continue
                cdn2 = list(cdn1)
                cdn2[pic] = nuc2
                cdn2 = ''.join(cdn2)
                aa2 = coda.translate_codon(cdn2)
                is_syn = aa1 == aa2
                sbs = f'{nuc1}>{nuc2}'
                data.append((pic, cdn1, cdn2, aa1, aa2, is_syn, sbs))
                i += 1

    df_changes = pd.DataFrame(data, columns=['pic', 'cdn1', 'cdn2', 'aa1', 'aa2', 'is_syn', 'sbs'])
    return df_changes


def nuc_spectrum_to_matrix(spec):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    M = np.zeros((4,4))
    for i1,n1 in enumerate(alphabet):
        for i2,n2 in enumerate(alphabet):
            if n1!=n2:
                M[i2,i1] = spec[f"{n1}>{n2}"]
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()
    # # will the diagonal with 'outflow' term to guarantee conservation of probability
    d = M.sum(axis=0)
    np.fill_diagonal(M,-d)

    return M


def cdn_spectrum_to_matrix(cdn_sbs):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    n = len(possible_codons)
    M = np.zeros((n, n))
    for i1,cdn1 in enumerate(possible_codons):
        for i2,cdn2 in enumerate(possible_codons):
            if cdn1!=cdn2:
                val = cdn_sbs[(cdn1, cdn2)] if (cdn1, cdn2) in cdn_sbs.index else 0.
                M[i2,i1] = val
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()
    # will the diagonal with 'outflow' term to guarantee conservation of probability
    d = M.sum(axis=0)
    np.fill_diagonal(M,-d)
    return M


def get_equilibrium_probabilities(M):
    evals, evecs = np.linalg.eig(M)
    # find zero eigenvalue
    ii = np.argmin(np.abs(evals))
    assert np.abs(evals[ii])<1e-10
    # pull out corresponding eigenvector, return normalized to sum_i p_i = 1
    p = evecs[:,ii]
    return p/p.sum()


def main():
    # derive all possible changes in the gencode 
    df_changes = collect_possible_changes()

    # load Kelley Harris spectra
    ms12 = pd.read_csv('./data/external/rates_by_clade.csv')
    ms12 = ms12[ms12.clade == '20A'].copy()
    ms12['mut_type'] = ms12['mut_type'].str.replace('to', '>')

    sbs2rate = ms12.set_index('mut_type')['rate'].to_dict()


    M = nuc_spectrum_to_matrix(sbs2rate)
    eq_prob = get_equilibrium_probabilities(M).astype(float)
    nucl_eq = pd.Series(dict(zip(alphabet, eq_prob)))
    nucl_eq.name = 'freq'
    nucl_eq.index.name = 'nucl'
    nucl_eq = nucl_eq.reset_index()
    nucl_eq['taxid'] = 0
    print(nucl_eq)


    df_changes['rate'] = df_changes['sbs'].map(sbs2rate)
    cdn_sbs = df_changes.groupby(['cdn1', 'cdn2'])['rate'].sum()
    M = cdn_spectrum_to_matrix(cdn_sbs)
    eq_prob = get_equilibrium_probabilities(M).astype(float)

    eq_freqs = pd.Series(dict(zip(possible_codons, eq_prob)))
    eq_freqs.name = 'freq'
    eq_freqs.index.name = 'cdn'
    eq_freqs = eq_freqs.reset_index()
    eq_freqs['aa'] = eq_freqs['cdn'].map(coda.translate_codon)
    eq_freqs['taxid'] = 0
    print(eq_freqs)
    # eq_freqs.to_csv('../data/equilibrium_freqs_20A.csv', index=False)


if __name__ == "__main__":
    main()
