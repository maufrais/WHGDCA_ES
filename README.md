# WHGDCA_ES
Within-Host Genomic Diversity of Candida albicans in Healthy Carriers: Implication for Microevolution Studies

## Running the test dataset
+ python src/snp_add_heterozygous_info.py -i data/snp_matrix.txt -o data/snp_matrix.withHTZ.txt
+ python src/snp_extract_microLOH.py -i data/snp_matrix.withHTZ.txt -o data/microLOH.txt

### Prerequisites
+ **Python      (2.7)**

## Authors
Corinne Maufrais
Bioinformatics and Biostatistics Hub C3BI - Institut Pasteur
Paris, France

