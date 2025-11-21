# Download SARS-CoV-2 genomes from gisaid automatically

Interface of gisaid site had changed and script `download_data.py` must be modified before run

## If you want to download specific SARS-CoV-2 records from [gisaid](https://www.epicov.org), follow the steps:

0. open gisaid site
1. set filters on Search page in EpiCoV (complete, high coverage etc.);
2. select records (just push upper checkbox button to select all records);
3. click button "Select" in the bottom (after that another subwindow will be opened and there you will see Accession Numbers of the records in the field);
4. click button "CSV" to download all records in nice format;
5. pass the csv file to script: `download_data.py -a accession_numbers.csv`
