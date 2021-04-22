# vcf2conStruct
R script to convert VCF format biallelic SNP data to the format required for the conStruct R package (https://github.com/gbradburd/conStruct)

usage example: ```Rscript vcf2conStruct.R --vcf my_file.vcf.gz --metadatafile my_meta.txt --outname my_conStruct```

OR: ```Rscript vcf2conStruct.R -v my_file.vcf.gz -m my_meta.txt -o my_conStruct```

In the above example, the output be saved as `my_conStruct.RData`. The object name will always be `conStruct_object` (after loading into R using `load("my_conStruct.RData")`.

Requires the following packages: vcfR, adegenet, fields, dplyr, argparser
run ```install.packages('vcfR', 'adegenet', 'fields', 'dplyr', 'argparser')```

Requires a metadata file with *no header* with columns in this order: `'sample population longitude latitude'` (can be tab or space delimited)
