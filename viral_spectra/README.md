# Notes

<!-- 
Need to filter small parts using seqkit filter -->

```bash
ls fasta/ | parallel seqkit rmdup -w 80 -s '<' fasta/{} '>' flt/{}

ls fasta/ | parallel mafft fasta/{} '>' msa/{}

ls msa/*.fa | parallel -j 6  iqtree2 -s {} -T 4 --prefix trees/{/.} -n 10 --ninit 10 --fast

```

txid3052310     >MG812674.1:602-7276 Mammarenavirus lassaense isolate 812285 segment L Z protein and L protein genes, complete cds
3052310(11620)
