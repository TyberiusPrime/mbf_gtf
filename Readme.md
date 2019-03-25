# mbf_gtf


Possibly the fastes Ensembl-GTF parser around 
(reads the 1GB human GTF in about 10s on my system).

Usage: mbf_gtf.parse_ensembl_gtf("filename.gtf", []) -> A dict of DataFrames

The file may be compressed with gzip - it must then end with ".gz".

The second parameter may be a list of 'features' to retrieve - getting 
just a subset can greatly improve performance.

Note that this is very ensembl specific, it does not deal with any other GTF
format, and that it throws away attributes that are repeated on the sub elements - 
ie. exons have only gene_id, not gene_name, gene_version, gene_....

The resulting coordinates are pythonic - ie. starting at 0 (ie. shifted -1 from
the values in the GTF).


This is part of the mbf_* family of bioinformatic libraries.
