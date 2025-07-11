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

    # read site rates
    site_rates_df = pd.read_csv('data/ref_sites_rates_cat4.csv')

    metrics_total = []
    for clade in tqdm.tqdm(clades_spectra.clade.unique(), desc='Clades'):
        spectrum_clade = clades_spectra[clades_spectra['clade'] == clade]
        exp_aa_subst, _ = prepare_exp_aa_subst(spectrum_clade, 'rate', 1)

        # Select clade OBS substitutions
        obs_clade = obs[obs['clade'] == clade].copy()
        obs_clade = obs_clade[(obs_clade.aa1 != '*') & (obs_clade.aa2 != '*')]

        # site_rates = obs_raw.query('noncoding == False').groupby('nt_site')['count'].sum()
        # site_rates = obs_clade.groupby('nt_site')['count'].sum()
        # genome_intervals = np.linspace(0, len(site_rates), 11).astype(int)
        # for i, (x1, x2) in enumerate(zip(genome_intervals[:-1], genome_intervals[1:])):

        for i, rcat in enumerate(range(1, site_rates_df.rate_cat.max()+1)):
            cur_sites = site_rates_df.query('rate_cat == @rcat').site.tolist()
            df_obs = obs_clade[obs_clade.nt_site.isin(cur_sites)]
            
            cur_aa_freqs_dct = get_site_specific_aa_counts(cur_sites)
            aa_subst = prepare_aa_subst(df_obs, exp_aa_subst, cur_aa_freqs_dct)
            
            cur_metrics = calc_metrics(aa_subst)
            cur_metrics['nt_sites'] = df_obs.nt_site.nunique()
            cur_metrics['clade'] = clade
            cur_metrics['rate_cat'] = rcat
            metrics_total.append(cur_metrics)

        # for total sites set
        aa_subst = prepare_aa_subst(obs_clade, exp_aa_subst, aa_freqs_total_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['nt_sites'] = obs_clade.nt_site.nunique()
        cur_metrics['clade'] = clade
        cur_metrics['rate_cat'] = 'any'
        metrics_total.append(cur_metrics)

    metrics_total_df = pd.DataFrame(metrics_total).set_index(['clade', 'rate_cat', 'nt_sites'])\
            .drop(['mape','wape','slope','intercept','pearson_corr','pearson_p',
                   'ks_stat','ks_p','rmse','log_likelihood'], axis=1)
    metrics_total_df.to_csv('data/fit_metrics_site_rate_cat4.csv', float_format='%g')


    print(metrics_total_df.reset_index().query('rate_cat == "any"')\
        .set_index('clade')[['r2', 'spearman_corr', 'mut_count']].round(2))


if __name__ == "__main__":
    main()