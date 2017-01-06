# PRS-on-SPARK
Polygenic risk score pipeline for large genotype data, generally imputed data

## Installation

To clone the repository, use 
```
git clone https://github.com/seriousNel/PRS-on-SPARK.git
```

## Requirements

The notebooks and scripts require the following to run :

+ spark-2.0.0 +
+ Pyhon 2.7

## What the pipeline does:
+ Calculate PRS from a genotype file (in .gen or .vcf format) and a GWAS file 


# Documentation for Command-line script PRS_run.py
## Parameters
A description of the parameters for the script can be obtained by typing: 
```
python PRS_run.py
```
which gives: 

```

positional arguments:
  GENO                  Name of the Genotype files, can be a name or path, or
                        name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  Output                The path and name for the output file

optional arguments:
  -h, --help            show this help message and exit
  --gwas_id GWAS_ID     Column number in your GWAS that contains SNP ID, with
                        first column being 0, default is 0
  --gwas_p GWAS_P       Column number in your GWAS that contains p-value, with
                        first column being 0, default is 7
  --gwas_or GWAS_OR     Column number in your GWAS that contains odds-
                        ratio/beta, with first column being 0, default is 5
  --gwas_a1 GWAS_A1     Column number in your GWAS that contains allele A1,
                        with first column being 0, default is 3. Allele A2 is
                        assumed to be at column [gwas_a1+1]
  --gwas_maf GWAS_MAF   Column number in your GWAS that contains frequency of
                        A1, with first column being 0, default is 10.
  --filetype {GEN,VCF}  The type of genotype file used as inputm choose
                        between VCF and GEN, default is VCF
  --thresholds THRESHOLDS [THRESHOLDS ...]
                        The p-value thresholds that controls which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. default is [0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001]
  --GWAS_delim GWAS_DELIM
                        Delimtier of the GWAS file, default is tab-delimiter
  --GWAS_no_header      Adding this parameter signals that there is no headers
                        for the GWAS. The default is to assume that GWAS has
                        column names
  --log_or              Adding this parameter tells the script to log the
                        effect sizes provided in the GWAS
  --check_ref           Adding this option tells the script to theck reference
                        allele when determining genoypte calls. Default is not
                        checking
  --app_name APP_NAME   Give your spark application a name. Default is PRS.
  --sample_file SAMPLE_FILE
                        path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_delim SAMPLE_DELIM
                        Delimiter of the sample file. Default is comma
  --sample_file_ID SAMPLE_FILE_ID [SAMPLE_FILE_ID ...]
                        Specify which columns in the sample file are used as
                        labels. Can use one integer to specify one column, or
                        multiple integers to specify multiple columns. Default
                        is the first column
  --sample_file_skip SAMPLE_SKIP
                        Specify how many lines to skip in the sample file,
                        i.e. which row do the labels start. Default is 1,
                        which assumes that the sample files has column names
                        and the labels start on the second line
  --use_maf             Use this paramter to tell the script to calculate MAF
                        in the provided propulation and compare it with MAF in
                        the GWAS, in order to check the reference alleles of
                        ambiguous SNPs (those whose A1 and A2 are reverese
                        complements). Not using this will result in ambiguous
                        SNPs be discarded. Default is not using MAF
  --log LOG             Specify the location of the log file. Default is no
                        log file


```
