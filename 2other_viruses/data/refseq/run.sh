
cut -d, -f1 acc2taxid.csv | tail +2 > acc.txt

#protein,gbff,genome
datasets download genome accession --inputfile acc.txt --include protein

# after that manually remove ORF1a polyprotein from CoV viruses and 
# >NP_058423.1 replicase [Transmissible gastroenteritis virus]