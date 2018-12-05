# GWAS Pipeline
This is the standard GWAS pipeline. Please notice that your GWAS sites is in buid37 coordinate system, if it is not the case, you can use the UCSC liftOver tool to perform the conversion to build37 system. Please also notice it is prefer that your GWAS dataset is forward strand, otherwise the SNPs which found in both reference panel and target panel that have an incompatible alleles types will be removed during the procedure of strand alignment check in this pipeline. 
## Requirement: 
please download tools to the folder "software": 
* plink (please download the latest plink1.9, otherwise you will waste a lot of time!!!)
* shapeit
* impute2

other pre-requirements:
* a cluster or you need a server with many cores
* genotyped data(plink format, transform to bed/bim/fam before running the pipeline if your genotyped data is ped/map)
* enough storage since very large files will be created during imputation procedure
* R gawk python3.7

## Procedures



* 1.users to edit the directory in config.py
* 2.run quality_control.py

Notice: If IDs are not matched in genotype and phenotype data,you should provide "mylist.txt" file.
* 3.run imputation.py
* 4.run post_imputation.py, a manhattan plot and a qqplot shown like the below will be created. 

Notice: it may takes a long time if you want to impute the whole chromosomes, all chromosomes are separated into more than 500 chunks totally with each chunk 5MB according to the physical position, also you need to have enough storage to save the outputs!!!

![Alt text](https://github.com/verasiwei/GWAS_python/blob/master/result/manhattan_jak2_4covs.png)
![Alt text](https://github.com/verasiwei/GWAS_python/blob/master/result/qqplot_jak2.png)









