# Notes

<!-- 
Need to filter small parts using seqkit filter -->

```bash
ls fasta/ | parallel seqkit rmdup -w 80 -s '<' fasta/{} '>' flt/{}

ls fasta/ | parallel mafft fasta/{} '>' msa/{}

ls msa/*.fa | parallel -j 6  iqtree2 -s {} -T 4 --prefix trees/{/.} -n 10 --ninit 10 --fast

```
