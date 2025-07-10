import os
from collections import Counter

from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
from pymutspec.io import read_genbank_ref
from pymutspec.annotation import CodonAnnotation

from utils import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst, 
)

wd = 'sars-cov-2_aa_subst'
if not os.getcwd().endswith(wd):
    os.chdir(wd)

coda = CodonAnnotation(1)

ref_sites_df = read_genbank_ref('data/NC_045512.2.gb')
ref_sites_df = ref_sites_df[ref_sites_df.Codon.notna()]
ref_sites_df['AA'] = ref_sites_df['Codon']\
    .apply(coda.translate_codon).map(amino_acid_codes)


def read_aa_counts_from_gb(path: str) -> dict:
    rec = next(SeqIO.parse(path, "genbank"))
    ref_df = pd.DataFrame([f.qualifiers for f in rec.features if f.type == "CDS"])
    ref_df.drop(columns=["locus_tag", "ribosomal_slippage", "codon_start", 
                        "db_xref", 'gene_synonym'], inplace=True)
    ref_df["gene"] = ref_df["gene"].apply(lambda x: x[0])
    ref_df["product"] = ref_df["product"].apply(lambda x: x[0])
    ref_df["protein_id"] = ref_df["protein_id"].apply(lambda x: x[0])
    ref_df["translation"] = ref_df["translation"].apply(lambda x: x[0])
    ref_df = ref_df[ref_df['product'] != 'ORF1a polyprotein']

    aa_counts_df = pd.DataFrame(ref_df.set_index('gene')['translation']\
                            .apply(Counter).to_dict()).T.fillna(0).astype(int)

    aa_counts_total_dct = aa_counts_df.rename(columns=amino_acid_codes).sum(0).to_dict()
    return aa_counts_total_dct


def get_site_specific_aa_counts(sites):
    aa_counts = ref_sites_df[ref_sites_df.Pos.isin(sites)]\
        .query('AA != "*"').AA.value_counts().to_dict()
    return aa_counts


def main():
    aa_freqs_total_dct = read_aa_counts_from_gb("data/NC_045512.2.gb")

    clades_spectra = pd.read_csv('data/bloom_etal/rates_by_clade.csv').query('subset == "all"')
    clades_spectra['Mut'] = clades_spectra['mut_type'].str.replace('to', '>')

    obs_raw = pd.read_csv('data/bloom_etal/aggregated.csv').query('subset == "all" & exclude == False')
    obs = obs_raw.query('synonymous == False & noncoding == False')\
        .drop(['synonymous', 'noncoding', 'four_fold_degenerate'], axis=1)

    def _same_aa_mut(aa_mutation: str):
        """filter syn mutations? and maybe common sites for several CDS"""
        variants = aa_mutation.split(';')
        return variants.count(variants[0]) == len(variants)
    
    obs = obs[obs.aa_mutation.apply(_same_aa_mut)]
    obs['aa1'] = obs['aa_mutation'].str[0]
    obs['aa2'] = obs['aa_mutation'].str[-1]

    # uniform_aa_content = obs.query('aa1 != "*"').groupby('aa1')[['count']].sum()\
    #     .assign(count=1)['count'].rename(index=amino_acid_codes).to_dict()
    # print(uniform_aa_content)

    metrics_total = []
    for clade in tqdm.tqdm(clades_spectra.clade.unique(), desc='Clades'):
        spectrum_clade = clades_spectra[clades_spectra['clade'] == clade]
        exp_aa_subst, _ = prepare_exp_aa_subst(spectrum_clade, 'rate', 1)

        # Select clade OBS substitutions
        obs_clade = obs[obs['clade'] == clade].copy()
        obs_clade = obs_clade[(obs_clade.aa1 != '*') & (obs_clade.aa2 != '*')]

        # site_rates = obs_raw.query('noncoding == False').groupby('nt_site')['count'].sum()
        site_rates = obs_clade.groupby('nt_site')['count'].sum()
        genome_intervals = np.linspace(0, len(site_rates), 11).astype(int)

        for i, (x1, x2) in enumerate(zip(genome_intervals[:-1], genome_intervals[1:])):
            cur_sites = site_rates.sort_values().iloc[x1: x2].index.values
            df_obs = obs_clade[obs_clade.nt_site.isin(cur_sites)]
            
            cur_aa_freqs_dct = get_site_specific_aa_counts(cur_sites)
            aa_subst = prepare_aa_subst(df_obs, exp_aa_subst, cur_aa_freqs_dct)
            
            cur_metrics = calc_metrics(aa_subst)
            cur_metrics['nt_sites'] = df_obs.nt_site.nunique()
            cur_metrics['clade'] = clade
            cur_metrics['sites_sample'] = i+1
            metrics_total.append(cur_metrics)

        # for total sites set
        aa_subst = prepare_aa_subst(obs_clade, exp_aa_subst, aa_freqs_total_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['nt_sites'] = obs_clade.nt_site.nunique()
        cur_metrics['clade'] = clade
        cur_metrics['sites_sample'] = 'total'
        metrics_total.append(cur_metrics)

    metrics_total_df = pd.DataFrame(metrics_total).set_index(['clade', 'sites_sample', 'nt_sites'])\
            .drop(['mape','wape','slope','intercept','pearson_corr','pearson_p',
                   'ks_stat','ks_p','rmse','log_likelihood'], axis=1)
    metrics_total_df.to_csv('data/fit_metrics_site_rates10.csv', float_format='%g')


    # _ = metrics_total_df[['spearman_corr','r2', 'wape', 'mut_count']]\
    #     .melt(ignore_index=False, var_name='metric').reset_index()
    # g = sns.catplot(data=_, sharey=False, kind='box', col='metric', col_wrap=2,
    #                 y='value', x='sites_sample',
    #                 palette='Set2', height=3, aspect=1.25,
    # )
    # g.set_titles('{col_name}')
    # g.set_xticklabels(g.axes[2].get_xticklabels(), rotation=-45)
    # g.set_xlabels('Sites')
    # g.set_ylabels('Metric value')
    # g.axes_dict['mut_count'].set_yscale('log')
    # g.legend.set_title('Site sampling\ncutoff, %')
    # g.savefig('./figures/fit_metrics_site_rates.pdf')

    # _ = metrics_total_df[['r2']]\
    #     .melt(ignore_index=False, var_name='metric').reset_index()
    # g = sns.catplot(data=_, sharey=False, kind='box', col='metric', col_wrap=1,
    #                 y='value', x='sites_sample',
    #                 palette='Set2', height=4, aspect=1.2,
    #                 # col_order=['accuracy']
    # )
    # g.set_titles('')
    # g.set_xticklabels(g.axes[0].get_xticklabels(), rotation=-25)
    # g.set_xlabels('Sites')
    # g.set_ylabels('$R^2$')
    # g.legend.set_title('Site sampling\ncutoff, %')
    # plt.grid(axis='y')
    # g.savefig('./figures/fit_metrics_site_rates_r2.pdf')

    print(metrics_total_df.reset_index().query('sites_sample == "total"')\
        .set_index('clade')[['r2', 'spearman_corr', 'mut_count']].round(2))


if __name__ == "__main__":
    main()