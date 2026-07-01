
## Analysis structure

1. derive viral spectra
2. compare them with reference ones + analyze differences and similarities
3. check spectra impact on amino acid substitutions
4. compare observed and expected amino acid frequencies differences
5. RdRp analysis


## Notes

<!-- 
Need to filter small parts using seqkit filter -->

```bash
ls fasta/ | parallel seqkit rmdup -w 80 -s '<' fasta/{} '>' flt/{}

ls fasta/ | parallel mafft fasta/{} '>' msa/{}

ls msa/*.fa | parallel -j 6  iqtree2 -s {} -T 4 --prefix trees/{/.} -n 10 --ninit 10 --fast

```
To obtain mutation spectra of viruses we utilize the NeMu pipeline (https://doi.org/10.1093/nar/gkae438).
First we indentify the protein sequence of the target gene in the selected virus using NCBI and then run the pipeline.
If NeMu fails to detect an appropriate outgroup for the chosen sequence, we provide a set of alighned sequenses obtained by tblastn. To indentify a suitable outgroup, we employ blastn with 'somewhat similar' parameter while excluded the virus of interest.
``` 



Data Flow Report for 2other_viruses
Purpose
This folder is organized around a small number of derived tables that are reused across notebooks. The main goal of cleanup should be to make these derived artifacts explicit and reduce repeated logic.

Current Backbone
1. Taxonomy and species selection
1.1get_viruses_linages.ipynb builds:
vir_linages.csv
species.csv
taxids.txt
Inputs:
df_annot_flt_en.csv
Role:
Defines virus taxonomy, lineage, and type labels used downstream.
2. Viral mutation spectrum dataset
1get_viral_spectra_dataset.ipynb builds:
viral_spectra_dataset.csv
Inputs:
ms12syn_all_virus.csv
taxid_virus_type.csv
viral_taxid_info_be.csv
species.csv indirectly through metadata joins
Role:
This is the central analysis table. Most downstream notebooks start here.
3. Observed amino-acid frequencies
get_viruses_nucl_freq_aa_freq.ipynb builds:
aminoacid_freq_all_virus.csv
Inputs:
RefSeq protein FASTA files under data/refseq
acc2taxid.csv
Role:
Provides observed amino-acid composition per virus, used later for observed-vs-equilibrium comparisons.
4. Equilibrium frequencies from spectra
get_viruses_nuc_and_aa_equilibrium_freq.ipynb builds:
aa_equilibrium_freq.csv
Legacy outputs in the same notebook:
aminoacid_eq_freq_all_virus.csv
nucl_eq_freq_all_virus.csv
aminoacid_eq_freq_sars_cov2.csv
nucl_eq_freq_sars_cov2.csv
Inputs:
viral_spectra_dataset.csv
ms12syn_all_virus.csv
Bloom_etal/rates_by_clade.csv
Role:
Converts mutation spectra into equilibrium amino-acid and nucleotide frequencies.
This notebook contains duplicated legacy logic and should be simplified.
5. Group-average equilibrium summary
2get_plot_ms12grouped_all_viruses.ipynb builds:
mean_gr_eq_freq.csv
ms12grouped.pdf
Inputs:
viral_spectra_dataset.csv
species.csv
Role:
Aggregates spectra by viral group and produces a reusable equilibrium summary.
6. Metadata assembly
prepare_table_of_vir_info.ipynb builds:
viral_spectra_long.csv
viral_meta_table.csv
viral_taxid_info_be.csv
Inputs:
viral_spectra_dataset.csv
ms12syn_all_virus.csv
taxid_virus_type.csv
vir_linages.csv
mean_gr_eq_freq.csv
aminoacid_freq_all_virus.csv
distance_to_eq.csv
Role:
Final table-building notebook; should remain a consumer of upstream artifacts, not a place for core derivations.
Downstream Consumers
Visualization and comparison notebooks
get_plot_aa_observed_freq_vs_eq_freq.ipynb
Uses aminoacid_freq_all_virus.csv
Uses aa_equilibrium_freq.csv
Produces observed-vs-equilibrium plots
4compare_aa_freqs_diff.ipynb
Uses viral_spectra_dataset.csv
Uses aminoacid_freq_all_virus.csv
Uses aa_equilibrium_freq.csv
Uses mean_gr_eq_freq.csv
Writes distance_to_eq.csv
3pca.ipynb
Uses species.csv
Uses viral_spectra_dataset.csv
2spectra_eda_plus_eq_freq_be.ipynb
Uses viral_spectra_dataset.csv
Uses vir_linages.csv
3.2model_fit_quality.ipynb
Uses viral_spectra_dataset.csv
Uses viral_taxid_info_be.csv
Uses distance_to_eq.csv
Uses virs_rnd_fit_metrics.csv
3.1apply_random_model.py
Uses viral_spectra_dataset.csv
Uses allmut_nemu.csv
Writes virs_rnd_fit_metrics.csv
Main Redundancies
The equilibrium-frequency logic is duplicated in get_viruses_nuc_and_aa_equilibrium_freq.ipynb in both the current and legacy sections.
4compare_aa_freqs_diff.ipynb appears to recompute a grouped equilibrium concept that already exists as mean_gr_eq_freq.csv.
Several notebooks are doing both transformation and plotting, which makes the data lineage harder to follow.
Some outputs are effectively archival or legacy, especially the older equilibrium CSVs in data.
Suggested Cleanup Order
Make viral_spectra_dataset.csv the single primary input for downstream analyses.
Centralize equilibrium calculations into one reusable function or script.
Keep one notebook for each major derived table:
lineage/species
spectra dataset
observed AA frequencies
equilibrium frequencies
metadata assembly
Move plotting-only logic out of data-building notebooks.
Delete or archive legacy outputs only after checking that no notebook still reads them.
Notebook State Note
The latest notebook summary for 4compare_aa_freqs_diff.ipynb shows the early cells executed successfully and later analytical cells still present but not executed. That notebook currently has in-memory variables like spectra, obs_aa_freqs, eq_aa_freqs, gr_aa_mean_content, and freqs_merged, which confirms it is a downstream comparison notebook rather than a source notebook.

If you want, I can next turn this into a cleaner dependency map with three columns: source, derived artifact, and consumer.