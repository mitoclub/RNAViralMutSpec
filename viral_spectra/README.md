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

```
To obtain mutation spectra of viruses we utilize the NeMu pipeline (https://doi.org/10.1093/nar/gkae438).
First we indentify the protein sequence of the target gene in the selected virus using NCBI and then run the pipeline.
If NeMu fails to detect an appropriate outgroup for the chosen sequence, we provide a set of alighned sequenses obtained by tblastn. To indentify a suitable outgroup, we employ blastn with 'somewhat similar' parameter while excluded the virus of interest.
 

