You may download a local copy of the latest gencode gff here to save time downloading.

You can get the latest versions with e.g.
```
cd data/annotations
gencode_release=49
#Download gencode data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencode_release}/gencode.v${gencode_release}.basic.annotation.gff3.gz
#Download GTEx data
wget https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz
cd ../../
```
Yes we know R has a cache etc. but sometimes 'AnnotationHub', 'BioMart' and 'gtexr' don't always play nicely.
