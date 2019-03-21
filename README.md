#mbf_gtf


Possibly the fastes Ensembl-GTF parser around 
(reads the 1GB human GTF in about 10s on my system).

Usage: mbf_gtf.parse_ensembl_gtf(filename, []) -> A dict of DataFrames

The second parameter may be a list of 'features' to retrieve - getting 
just a subset can greatly improve performance.

Note that this is very ensembl specific, it does not deal with any other GTF
format, and that it throws away attributes that are repeated on the sub elements - 
ie. exons have only gene_id, not gene_name, gene_version, gene_....


