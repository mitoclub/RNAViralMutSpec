# Scripts for separate omicron strain processing


## Notes

1. At the filtering of omicron sequences reference seq dropped due to low quality. Consequently, need to append it to first place separately.

```bash 
> awk 'BEGIN {RS=">"} /EPI_ISL_6752027/ {printf ">"$o}' data/omicron/gisaid_omicron_13-01-22.fasta > data/external/omicron_ref.fasta
> cat data/external/omicron_ref.fasta data/omicron/sequences.filtered.fasta > data/omicron/sequences.filtered.fasta
```
