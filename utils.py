import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chisquare, ks_2samp, pearsonr, spearmanr, uniform
from scipy.special import kl_div
from sklearn.metrics import mean_squared_error
from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_codons

alphabet = 'ACGT'
amino_acid_codes = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}


def collect_possible_changes(gc):
    coda = CodonAnnotation(gc)
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


def prepare_exp_aa_subst(spectrum: pd.DataFrame, rate_col='rate', gc=1, save_path=None):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index('Mut')[rate_col].to_dict()

    df_changes['rate'] = df_changes['sbs'].map(spectrum_dict)
    df_changes['aa1'] = df_changes['aa1'].map(amino_acid_codes)
    df_changes['aa2'] = df_changes['aa2'].map(amino_acid_codes)

    ## Calculate expected AA substitutions matrix
    exp_aa_subst = df_changes[(df_changes.aa1 != '*')&(df_changes.aa2 != '*')]\
        .groupby(['aa1', 'aa2'])['rate'].sum().reset_index()
    
    if save_path:
        exp_aa_subst.to_csv(save_path, float_format='%g', index=False)
    exp_aa_subst_matrix = exp_aa_subst.pivot(index='aa1', columns='aa2', values='rate').fillna(0.)
    return exp_aa_subst, exp_aa_subst_matrix


def prepare_aa_subst(obs_df: pd.DataFrame, exp_aa_subst: pd.DataFrame, ref_aa_freqs: dict):
    aa_subst = obs_df.groupby(['aa1', 'aa2'])['count'].sum().rename('nobs').reset_index()
    aa_subst = aa_subst[(aa_subst['aa1'] != aa_subst['aa2']) & (aa_subst['aa1'] != '*') & (aa_subst['aa2'] != '*')]
    aa_subst['aa1'] = aa_subst['aa1'].map(amino_acid_codes)
    aa_subst['aa2'] = aa_subst['aa2'].map(amino_acid_codes)
    aa_subst['ref_aa1_freq'] = aa_subst['aa1'].map(ref_aa_freqs) / sum(ref_aa_freqs.values())
    aa_subst['nobs_scaled'] = aa_subst['nobs'] / aa_subst['ref_aa1_freq']
    aa_subst['nobs_scaled'] = aa_subst['nobs_scaled'] / aa_subst['nobs_scaled'].sum() * aa_subst['nobs'].sum()
    aa_subst = aa_subst.merge(exp_aa_subst.rename(columns={'rate': 'rate_exp'}), 'right').fillna(0)
    aa_subst = aa_subst[aa_subst['aa1'] != aa_subst['aa2']]
    aa_subst['nexp'] = aa_subst['rate_exp'] / aa_subst['rate_exp'].sum() * aa_subst['nobs_scaled'].sum()
    aa_subst['diff'] = aa_subst['nobs_scaled'] - aa_subst['nexp']
    aa_subst['nobs_freqs'] = aa_subst['nobs_scaled'] / aa_subst['nobs_scaled'].sum()
    aa_subst['nexp_freqs'] = aa_subst['nexp'] / aa_subst['nexp'].sum()
    return aa_subst


def calc_accuracy(y_true, y_pred):
    '''Acc = Sum(TP1, TP2, …, TPn)/total_obs_cnt'''
    TP = np.min([y_true, y_pred], axis=0)
    total_obs_cnt = y_true.sum()
    acc = TP.sum() / total_obs_cnt
    return acc


def calc_f1(y_true, y_pred):
    '''Acc = Sum(TP1, TP2, …, TPn)/total_obs_cnt'''
    f1_weighted = 0.
    f1_macro = 0.
    total = y_true.sum()
    n = len(y_true)

    for o, e in zip(y_true, y_pred):
        if o < e:
            tp = o
            fp = e - o
            fn = 0
        elif o > e:
            tp = e
            fp = 0
            fn = o - e
        else:
            tp = o
            fp = fn = 0

        f1 = (2 * tp) / (2 * tp + fp + fn)

        f1_macro += f1 / n
        f1_weighted += f1 * (o / total)

    return f1_macro, f1_weighted


def calc_metrics(aa_subst: pd.DataFrame):
    y_true = aa_subst.nobs_scaled
    y_pred = aa_subst.nexp

    # 2. Kolmogorov-Smirnov test
    ks_stat, ks_p = ks_2samp(y_true, y_pred)

    # 3. Log-Likelihood
    log_likelihood = np.sum(aa_subst.nobs_freqs * np.log(aa_subst.nexp_freqs) + \
                                (1 - aa_subst.nobs_freqs) * np.log(1 - aa_subst.nexp_freqs))

    # 4. RMSE
    rmse = mean_squared_error(y_true, y_pred) ** 0.5

    # Spearman's rank correlation coefficient
    spearman_corr, spearman_p = spearmanr(y_true, y_pred)
    pearson_corr, pearson_p = pearsonr(y_true, y_pred)

    # 6. KL-divergence
    kl_divergence = np.sum(kl_div(aa_subst.nobs_freqs, 
                                aa_subst.nexp_freqs))
    
    # total number of ns mutations
    mut_count = np.sum(y_true)

    acc = calc_accuracy(y_true, y_pred)
    try:
        f1_macro, f1_weighted = calc_f1(y_true, y_pred)
    except ZeroDivisionError:
        print('Cannot calculate F1 score, division by zero')
        f1_macro, f1_weighted = np.nan, np.nan

    metrics = {
        'ks_stat': ks_stat,
        'ks_p': ks_p,
        'log_likelihood': log_likelihood,
        'rmse': rmse,
        'rmse_scaled': rmse / mut_count,
        'spearman_corr': spearman_corr,
        'spearman_p': spearman_p,
        'pearson_corr': pearson_corr,
        'pearson_p': pearson_p,
        'kl_divergence': kl_divergence,
        # 'bootstrap_p': bootstrap_p,
        'accuracy': acc,
        'f1_macro': f1_macro,
        'f1_weighted': f1_weighted,
        'mut_count': mut_count,
    }
    return metrics


def plot_exp_heatmap(exp_aa_subst_matrix: pd.DataFrame, save_path: str, show=True):
    """
    Plot expected amino acid substitution matrix.
    """
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()
    max_flow_value = max(freqs_to.max(), freqs_from.max())+0.1

    fig, axs = plt.subplots(2, 3, figsize=(11, 10), 
                            width_ratios=[0.1, 1, .05], height_ratios=[1, 0.1])
    sns.heatmap(exp_aa_subst_matrix, annot=False, fmt=".2f", 
                ax=axs[0, 1], cbar_ax=axs[0, 2], cmap='coolwarm', 
                cbar_kws={'label': 'Substitution rate'}, 
                mask=exp_aa_subst_matrix==0,
    )
    axs[0, 1].set_ylabel('')
    axs[0, 1].set_xlabel('')
    axs[0, 1].set_xticks([])
    axs[0, 1].set_yticks([])
    # axs[0, 1].set_title('Expected substitution rates between amino acids')

    sns.barplot(freqs_from.reset_index(), y='aa1', x=0, ax=axs[0, 0],
                color='lightgray', edgecolor='black')
    axs[0, 0].set_ylabel('From', fontsize=14)
    axs[0, 0].set_xlabel('')
    axs[0, 0].spines['top'].set_visible(False)
    axs[0, 0].spines['right'].set_visible(False)
    axs[0, 0].spines['bottom'].set_visible(False)
    axs[0, 0].spines['left'].set_visible(False)
    axs[0, 0].invert_xaxis()
    axs[1, 1].set_ylim(0, max_flow_value)
    axs[0, 0].set_xticks([])
    axs[0, 0].set_yticklabels(axs[0, 0].get_yticklabels(), fontsize=12)

    sns.barplot(freqs_to.reset_index(), x='aa2', y=0, ax=axs[1, 1],
                color='lightgray', edgecolor='black')
    axs[1, 1].set_ylabel('')
    axs[1, 1].set_xlabel('To', fontsize=14)
    axs[1, 1].spines['top'].set_visible(False)
    axs[1, 1].spines['right'].set_visible(False)
    axs[1, 1].spines['bottom'].set_visible(False)
    axs[1, 1].spines['left'].set_visible(False)
    axs[1, 1].invert_yaxis()
    axs[1, 1].set_ylim(max_flow_value, 0)
    axs[1, 1].set_yticks([])
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), fontsize=12)

    axs[1, 0].remove()
    axs[1, 2].remove()

    plt.tight_layout()
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_aa_eq_freqs(exp_aa_subst_matrix, save_path: str, show=True, figsize=(12, 3)):
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()

    _freqs_to = (freqs_to / freqs_to.sum()).copy()
    _freqs_from = (freqs_from / freqs_from.sum()).copy()
    _freqs_to.index.name = 'aa'
    _freqs_from.index.name = 'aa'

    plt.figure(figsize=figsize)
    ax = sns.barplot((_freqs_to - _freqs_from).sort_values(ascending=False).reset_index(), 
                x='aa', y=0, color='lightgray', edgecolor=".65")

    for bar in ax.patches:
        if bar.get_height() > 0.05:
            bar.set_color('orangered')
            bar.set_edgecolor('0.02')
        elif bar.get_height() > 0.02:
            bar.set_color('lightcoral')
            bar.set_edgecolor('0.05')
        elif bar.get_height() < -0.05:
            bar.set_color('royalblue')
            bar.set_edgecolor('0.02')
        elif bar.get_height() < -0.02:
            bar.set_color('lightblue')
            bar.set_edgecolor('0.45')

    plt.legend([
        plt.Rectangle((0,0),1,1,fc="orangered", edgecolor = '0.45'), 
        plt.Rectangle((0,0),1,1,fc="lightcoral", edgecolor = '0.45'), 
        plt.Rectangle((0,0),1,1,fc='lightblue', edgecolor = '0.45'),
        plt.Rectangle((0,0),1,1,fc='royalblue', edgecolor = '0.45'),
        ],[
            'Top Gainers', 'Gainers', 'Loosers', 'Top Loosers'], 
        loc='upper right')

    # plt.title(f'Expected aa subst based assignment')
    plt.ylabel('Rates [To - From]', fontsize=13)
    plt.xlabel('')
    plt.xticks(fontsize=12)
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()


def plot_obs_vs_exp(aa_subst, save_path, show=True):
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator
    aa_subst = aa_subst.copy()
    aa_subst['nobs_freqs_log'] = np.log10(aa_subst['nobs_freqs'])
    aa_subst['nexp_freqs_log'] = np.log10(aa_subst['nexp_freqs'])

    cor_res = spearmanr(aa_subst['nobs_freqs'], aa_subst['nexp_freqs'])

    plt.figure(figsize=(5, 5))
    ax = sns.regplot(aa_subst[aa_subst['nobs_freqs']>0], 
                color='blue', scatter_kws={'alpha':0.5},
                x='nobs_freqs_log', y='nexp_freqs_log', label='20A',)

    # ticks = np.log10(np.array([1e-4, 1e-3, 1e-2, 1e-1]))
    ticks = np.linspace(-4, -1, 4)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    formatter = lambda x, pos: '10$^{' + str(int(x)) + '}$'
    ax.get_xaxis().set_major_formatter(formatter)
    ax.get_yaxis().set_major_formatter(formatter)
    # ax.get_xaxis().set_minor_locator(LogLocator(-10, 'auto'))

    plt.text(-2, -4., 
            f"Spearman: {cor_res.correlation:.2f}\np-value: {cor_res.pvalue:.2g}", 
            fontsize=10)
    plt.plot([-4, -1], [-4, -1], color='black', linestyle='--', label='y=x')
    plt.xlabel('Observed AA frequencies')
    plt.ylabel('Expected AA frequencies')
    plt.legend()
    plt.savefig(save_path)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close()