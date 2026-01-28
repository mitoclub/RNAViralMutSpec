# RNA viruses amino acid substitutions analysis

## Environment

Create conda environment:

```bash
mamba env create -f environment.yml
mamba activate amino-acid-shift
pip install --no-deps pymutspec  # temporary solution until pymutspec is on conda-forge
pip install ete3 pytest PyYAML click legacy-cgi
```

## Structure

1. [SARS-CoV-2 amino acid substitutions ananlysis](./1sars-cov-2/)
2. [RNA viruses mutational spectra derivation and analysis](./2other_viruses/)

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
