from collections import Counter

from Bio import SeqIO
import pandas as pd
import tqdm

from utils import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst, 
    plot_obs_vs_exp,
)


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


def main():
    # TODO read true viral counts of amino acids
    aa_freqs_total_dct = read_aa_counts_from_gb("../sars-cov-2_aa_subst/data/NC_045512.2.gb")

    viral_spectra = pd.read_csv('./data/viral_spectra_dataset.csv').query('df == "nemu"')
    viral_spectra = viral_spectra.melt(
        ['Type', 'taxname', 'virusname'], 
        viral_spectra.columns[:12].values, 'Mut', 'Rate'
    ).sort_values(['virusname', 'Mut'])

    _mut_all = pd.read_csv('./data/allmut_nemu.csv')
    obs = _mut_all[_mut_all['Label'] == 0].rename(
        columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})

    metrics_total = []
    for (vir, group), obs_vir in tqdm.tqdm(
            viral_spectra.groupby(['virusname', 'Type']), 
            desc='Viruses'):
        cur_spectrum = viral_spectra[viral_spectra['virusname'] == vir]
        exp_aa_subst, _ = prepare_exp_aa_subst(cur_spectrum, 'Rate', 1)
        # Select vir OBS substitutions
        obs_vir = obs[obs['virusname'] == vir].copy()
        obs_vir = obs_vir[(obs_vir.aa1 != '*') & (obs_vir.aa2 != '*')]

        if obs_vir['count'].sum() < 500:
            continue

        # for total sites set
        aa_subst = prepare_aa_subst(obs_vir, exp_aa_subst, aa_freqs_total_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['Type'] = group
        cur_metrics['virusname'] = vir
        metrics_total.append(cur_metrics)

        plot_obs_vs_exp(aa_subst, f'./figures/neutral_model_fit/{group}_{vir}_obs_vs_exp.pdf', show=False)

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index(['Type', 'virusname', ])
    metrics_total_df.to_csv('data/virs_fit_metrics.csv', float_format='%g')
    print(metrics_total_df)


if __name__ == "__main__":
    main()