import os
from collections import Counter

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import tqdm

from utils import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst, 
    plot_obs_vs_exp,
)

debug = True


wd = 'viral_spectra'
if not os.getcwd().endswith(wd):
    os.chdir(wd)


def read_aa_counts_from_files(viral_spectra: pd.DataFrame) -> pd.DataFrame:
    indir = './data/nemu_inputs'
    all_files = os.listdir(indir)

    data = []
    for _, row in viral_spectra.iterrows():
        taxid = row['taxid']
        virusname = row['virusname']
        fasta_file = None
        for file in all_files:
            if file.startswith(str(taxid)):
                fasta_file = os.path.join(indir, file)
                break
        
        rec = next(SeqIO.parse(fasta_file, "fasta"))
        if rec.count('V') + rec.count('L') > 0:
            print(f'{fasta_file} is a protein')
            aa_seq = str(rec.seq)
        else:
            raise NotImplementedError()
        
        data.append({
            'taxid': taxid,
            'virusname': virusname,
            'file': fasta_file,
            'aa_seq': aa_seq,
        })

    df = pd.DataFrame(data)
    aa_counts_df = pd.DataFrame(
        df.set_index('virusname')['aa_seq'].apply(Counter).to_dict()
    ).T.fillna(0).astype(int).rename(columns=amino_acid_codes)
    return aa_counts_df


def main():
    viral_spectra_raw = pd.read_csv('./data/viral_spectra_dataset.csv').query('df == "nemu"')
    viral_spectra = viral_spectra_raw.melt(
        ['Type', 'taxname', 'virusname'], 
        viral_spectra_raw.columns[:12].values, 'Mut', 'Rate'
    ).sort_values(['virusname', 'Mut'])

    # read amino acid freqs from protein files
    aa_freqs = read_aa_counts_from_files(viral_spectra_raw)

    _mut_all = pd.read_csv('./data/allmut_nemu.csv')
    obs = _mut_all[_mut_all['Label'] == 0].rename(
        columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})

    metrics_total = []
    i = 0
    nrows, ncols = 3, 3
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4-1.5))
    for (vir, group), obs_vir in tqdm.tqdm(
            viral_spectra.groupby(['virusname', 'Type']), 
            desc='Viruses'):
        cur_spectrum = viral_spectra[viral_spectra['virusname'] == vir]
        exp_aa_subst, _ = prepare_exp_aa_subst(cur_spectrum, 'Rate', 1)
        # Select vir OBS substitutions
        obs_vir = obs[obs['virusname'] == vir].copy()
        obs_vir = obs_vir[(obs_vir.aa1 != '*') & (obs_vir.aa2 != '*')]

        if obs_vir['count'].sum() < 500 or vir == 'SBV':
            continue

        # for total sites set
        cur_aa_freqs_dct = aa_freqs.loc[vir].to_dict()
        aa_subst = prepare_aa_subst(obs_vir, exp_aa_subst, cur_aa_freqs_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['Type'] = group
        cur_metrics['virusname'] = vir
        metrics_total.append(cur_metrics)

        if obs_vir['count'].sum() > 2000:
            ax = axs[i // ncols, i % ncols]
            plot_obs_vs_exp(
                aa_subst, ax=ax, show=False, 
                text=f"{vir} ({group})", text_x=-2.2, text_y=-4.)
            i += 1
    
    # # unshow empty plots
    # axs[i // ncols, i % ncols].axis('off')
    # for x in range(2):
    #     i += 1
    #     axs[i // ncols, i % ncols].axis('off')
        

    plt.tight_layout()
    plt.savefig('./figures/neutral_model_fit.pdf', dpi=300)
    plt.close()

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index(['Type', 'virusname', ])
    metrics_total_df.to_csv('data/virs_fit_metrics.csv', float_format='%g')
    print(metrics_total_df)


if __name__ == "__main__":
    main()