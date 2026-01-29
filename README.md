# RNA viruses amino acid substitutions analysis

## Environment

Create conda environment:

```bash
mamba env create -f environment.yml
mamba activate amino-acid-shift
pip install ete3 pytest PyYAML click legacy-cgi
pip install --no-deps pymutspec  # temporary solution until pymutspec is on conda-forge
```

## Structure

1. [SARS-CoV-2 amino acid substitutions ananlysis](./1sars-cov-2/)
2. [RNA viruses mutational spectra derivation and analysis](./2other_viruses/)

| Notebook | Short description | Outline |
| --- | --- | --- |
| [1sars-cov-2/1Neutral_scenario.ipynb](1sars-cov-2/1Neutral_scenario.ipynb) | Neutral-model analyses: compute equilibrium freqs and compare models. | |
| [1sars-cov-2/2.1observed_sbs_fit.ipynb](1sars-cov-2/2.1observed_sbs_fit.ipynb) | Fit observed substitution patterns, compute fit metrics and plots. | |
| [1sars-cov-2/2.2fitness_analysis.ipynb](1sars-cov-2/2.2fitness_analysis.ipynb) | Analyze Δfitness of substitutions across clades using Bloom et al. data. | |
| [1sars-cov-2/2.3changes_during_COVID-19.ipynb](1sars-cov-2/2.3changes_during_COVID-19.ipynb) | Analyze genotype data and changes during COVID-19 (genotypes2025). | |
| [1sars-cov-2/absorb_r2.ipynb](1sars-cov-2/absorb_r2.ipynb) | Analyze absorption metrics and R² relationships for model fits. | |
| [1sars-cov-2/model.ipynb](1sars-cov-2/model.ipynb) | Model linking nucleotide mutation rates to expected AA substitutions. | |
| [1sars-cov-2/verify_model.ipynb](1sars-cov-2/verify_model.ipynb) | Verify equilibrium calculations (codon/AA level) and test stop-codon handling. | |
| [2other_viruses/1.1get_viruses_linages.ipynb](2other_viruses/1.1get_viruses_linages.ipynb) | Extract viral taxonomy/lineages and select species for analyses. | |
| [2other_viruses/1get_viral_spectra_dataset.ipynb](2other_viruses/1get_viral_spectra_dataset.ipynb) | Assemble viral spectra dataset (NEMU + Bloom) and tidy data. | |
| [2other_viruses/2spectra_eda_plus_eq_freq_be.ipynb](2other_viruses/2spectra_eda_plus_eq_freq_be.ipynb) | EDA of mutation spectra and equilibrium-frequency computations. | |
| [2other_viruses/2get_plot_ms12grouped_all_viruses.ipynb](2other_viruses/2get_plot_ms12grouped_all_viruses.ipynb) | Plot ms12 mutational spectra grouped by virus type/family. | |
| [2other_viruses/2get_cossim_between_nemu_Bloom.ipynb](2other_viruses/2get_cossim_between_nemu_Bloom.ipynb) | Compute cosine similarity between NEMU dataset and Bloom et al. spectra. | |
| [2other_viruses/3pca.ipynb](2other_viruses/3pca.ipynb) | PCA and classification of viral spectra; feature importance and clustering. | |
| [2other_viruses/3.2model_fit_quality.ipynb](2other_viruses/3.2model_fit_quality.ipynb) | Assess neutral-model fit quality and mutation summaries across viruses. | |
| [2other_viruses/4compare_aa_freqs_diff.ipynb](2other_viruses/4compare_aa_freqs_diff.ipynb) | Compare AA frequency differences between viral groups and plot results. | |
| [2other_viruses/prepare_table_of_vir_info.ipynb](2other_viruses/prepare_table_of_vir_info.ipynb) | Build table of virus metadata and compute summary statistics for spectra. | |
| [2other_viruses/get_viruses_nucl_freq_aa_freq.ipynb](2other_viruses/get_viruses_nucl_freq_aa_freq.ipynb) | Compute amino-acid and nucleotide frequencies from viral reference sequences. | |
| [2other_viruses/get_viruses_nuc_and_aa_equilibrium_freq.ipynb](2other_viruses/get_viruses_nuc_and_aa_equilibrium_freq.ipynb) | Compute nucleotide and AA equilibrium frequencies across viruses. | |
| [2other_viruses/get_plot_aa_observed_freq_vs_eq_freq.ipynb](2other_viruses/get_plot_aa_observed_freq_vs_eq_freq.ipynb) | Plot observed vs equilibrium amino-acid frequencies (per virus/group). | |
| [extra/human_germline.ipynb](extra/human_germline.ipynb) | Analyze human germline mutation spectra and per-gene AA equilibrium. | |
| [extra/nuc_eq_example.ipynb](extra/nuc_eq_example.ipynb) | Examples and calculations of nucleotide-equilibrium probabilities. | |
| [extra/random_spectra_thrend.ipynb](extra/random_spectra_thrend.ipynb) | Random-spectrum simulations and their effect on equilibrium AA frequencies. | |
| [extra/rate_vs_freq.ipynb](extra/rate_vs_freq.ipynb) | Investigate relationships between mutation rates and observed frequencies. | |

<!-- 1sars-cov-2/1Neutral_scenario.ipynb:    "## Load reference and calc amino acid freqs"
1sars-cov-2/1Neutral_scenario.ipynb:    "## read external clades spectra\n",
1sars-cov-2/1Neutral_scenario.ipynb:    "## Calculate expected AA substitutions matrix \n",
1sars-cov-2/1Neutral_scenario.ipynb:    "## Estimate neutrality of the spectrum"
1sars-cov-2/1Neutral_scenario.ipynb:    "## Spectrum of amino acid substitutions"
1sars-cov-2/1Neutral_scenario.ipynb:    "## Calculate amino acid equilibrium frequencies"
1sars-cov-2/1Neutral_scenario.ipynb:    "## Calculate codon and AA equilibrium frequencies\n",
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## Imports\n",
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## Load data\n",
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## read external clades spectra\n",
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## stuff"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Plots"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### qualitative model"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## Explore OBS"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Basic features"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### terminal vs non-terminal spectrum and mutations"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## compare obs and exp with chemical distance"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## Fitting quality on mutations subsets"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Clades"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Compare with random spectra"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Codon positions"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "    ## Calculate expected AA substitutions matrix\n",
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "#### Estimate grantham distances in substitutions from different codon positions"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Genes"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "#### Check now on ORF1ab proteins"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "#### Mut distribution in clades"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "#### Mut site distribution in clades"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### How many mutations required for unbiased result?"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## Site rate categories"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Categorize sites by mutation rate"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Explore Cat-specific spectra"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Explore categories"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "#### compare rate categories with \"landscape\" categories"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Explore hotspots (5th category)"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Fitting observed substitution spectra by site-specific rates"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Distance to equilibrium VIZ"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "### Calculate aa content of most variable and constant sites"
1sars-cov-2/2.1observed_sbs_fit.ipynb:    "## view the scatterplots with negative R2"
1sars-cov-2/2.2fitness_analysis.ipynb:    "## Load data"
1sars-cov-2/2.2fitness_analysis.ipynb:    "## read external clades spectra\n",
1sars-cov-2/2.2fitness_analysis.ipynb:    "## Delta Fitness of substitutions"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### Explore variability of mutations delta fitneses in clades"
1sars-cov-2/2.2fitness_analysis.ipynb:    "## Rate cats"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### asses model"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### image for paper"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### Check that 22 muttypes are non-randomly decrease R2"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### try different cutoffs"
1sars-cov-2/2.2fitness_analysis.ipynb:    "## Fitness of distinct AA"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### Prepare site samples of different fitness landscapes"
1sars-cov-2/2.2fitness_analysis.ipynb:    "## initial classification with some assumtions\n",
1sars-cov-2/2.2fitness_analysis.ipynb:    "### Assess model on different landscapes"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### Plot model quality for landscapes"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### plot distance to equilibrium in different landscapes"
1sars-cov-2/2.2fitness_analysis.ipynb:    "### AA delta_fitness AND AA fitness "
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "## Prepare data"
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "### Load ref"
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "### Load genotypes from 2025"
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "### All proteins ???"
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "## Compare expected and observed freq changes"
1sars-cov-2/2.3changes_during_COVID-19.ipynb:    "## compare simulation equilibrium with true one\n",
2other_viruses/1.1get_viruses_linages.ipynb:    "## Select species for dataset"
2other_viruses/1get_viral_spectra_dataset.ipynb:    "### add short virus names"
2other_viruses/3.2model_fit_quality.ipynb:    "## Explore mutations features"
2other_viruses/3.2model_fit_quality.ipynb:    "## Neutral model fit results"
2other_viruses/3pca.ipynb:    "## PCA"
2other_viruses/3pca.ipynb:    "## Mutations that most affected components"
2other_viruses/3pca.ipynb:    "## PCA with complement ssRNA-"
2other_viruses/3pca.ipynb:    "## Mutations that most affected components"
2other_viruses/4compare_aa_freqs_diff.ipynb:    "### compare with spectrum"
2other_viruses/4compare_aa_freqs_diff.ipynb:    "## PCA"
2other_viruses/4compare_aa_freqs_diff.ipynb:    "## Distance to eq"
2other_viruses/get_plot_aa_observed_freq_vs_eq_freq.ipynb:    "## Observed frequency vs Equilibrium frequency"
2other_viruses/get_plot_aa_observed_freq_vs_eq_freq.ipynb:    "## Distance between observed and equilibrium aa freqs"
2other_viruses/get_plot_aa_observed_freq_vs_eq_freq.ipynb:    "## Distance between nucleotids and aminoacids"
2other_viruses/get_viruses_nuc_and_aa_equilibrium_freq.ipynb:    "## Add Sars-cov-2 (count early)"
2other_viruses/get_viruses_nucl_freq_aa_freq.ipynb:    "## For aminoacids"
2other_viruses/get_viruses_nucl_freq_aa_freq.ipynb:    "## For nucleotids"
2other_viruses/get_viruses_nucl_freq_aa_freq.ipynb:    "## For 4f codons"
extra/human_germline.ipynb:    "## Derive equilibium for each gene"
extra/human_germline.ipynb:    "## Apply neutral model genes-wise"
extra/human_germline.ipynb:    "## compare obs and eq aa freqs"
extra/human_germline.ipynb:    "## How it looks on random spectra"
extra/nuc_eq_example.ipynb:    "## Viruses check"
extra/nuc_eq_example.ipynb:    "### 12 --> 6 --> 12 on aymmetric spectrum"
extra/nuc_eq_example.ipynb:    "#### Cov20A (highly asymmetric)"
extra/nuc_eq_example.ipynb:    "#### aka spectrum (more symmetric)"
extra/random_spectra_thrend.ipynb:    "## Load rdrp freqs"
extra/random_spectra_thrend.ipynb:    "## Random spectra consequenses"
extra/random_spectra_thrend.ipynb:    "### Tests"
extra/random_spectra_thrend.ipynb:    "### Precalculated"
extra/random_spectra_thrend.ipynb:    "### CCA"
extra/random_spectra_thrend.ipynb:    "### PCA of amino acid freqs"
extra/random_spectra_thrend.ipynb:    "## Load RdRp amino acid freqs" -->

## Data

### data of `Jesse D Bloom, Annabel C Beichman, Richard A Neher, Kelley Harris, Evolution of the SARS-CoV-2 Mutational Spectrum, Molecular Biology and Evolution, Volume 40, Issue 4, April 2023, msad085, https://doi.org/10.1093/molbev/msad085`

- image of spectra https://jbloomlab.github.io/SARS2-mut-spectrum/rates-by-clade.html
- spectra table https://github.com/jbloomlab/SARS2-mut-spectrum/blob/main/results/synonymous_mut_rates/rates_by_clade.csv
- other viruses spectra https://github.com/jbloomlab/SARS2-mut-spectrum/blob/main/results/other_virus_spectra/other_virus_spectra.json


- [Serratus](https://serratus.io/trees) database

## References

- [NAR paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaa1053/5961787) about sars-cov-2 genome secondary structure
- sars-cov-2 replication [paper 1](https://www.nature.com/articles/s41579-020-00468-6) and [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7122471/)
- sars-cov-2 life circle [paper](https://www.nature.com/articles/s41579-020-00468-6)
- statistics decision tree [pingouin](https://pingouin-stats.org/guidelines.html?highlight=krus#non-parametric)
- [sars-cov-2 variants](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1036501/Technical_Briefing_29_published_26_November_2021.pdf)
- [Serratus](https://serratus.io/trees) database
- [/refseq/release/viral](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/)
- [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) db
