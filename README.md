# PRS-on-SPARK
Polygenic risk score pipeline for large genotype data, generally imputed data.
Published in BMC bioinformatics
[PRS-on-Spark (PRSoS): a novel, efficient and flexible approach for generating polygenic risk scores](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2289-9)



## Installation

To clone the repository, use 
```
git clone https://github.com/seriousNel/PRS-on-SPARK.git
```

## Software Requirements

The notebooks and scripts require the following to run :
+ spark-2.0.0 +
+ Pyhon 2.7

Instruction of installing Apache Spark on linux can be found [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/)

The prerequisite to install spark are :
+ java8 (oracle)
+ scala 
+ sbt

Some extra libraries are required for regression and plotting. To install them, first make sure pip is installed on your computer, then type:
```
cd PRS-on-SPARK
pip install -r requirements.txt
```

## What the pipeline does:
+ Calculate PRS from a genotype file (in .gen or .vcf format) and a GWAS file
+ Correct the strand alignment descrpencies between genotype  and GWAS data. 

## What it cannot do :
+ Performs quality control of genotype data

## Default format :
### GWAS
By default, the GWAS should have the same format as that of a GWAS file obtained from Psychiatric Genomics Consortium (PGC). 

|     snpid|   pval|    or| a1| a2|    CEUaf|
|----------|-------|------|---|---|---------|
| rs3131972| 0.2032| 1.047|  A|  G|  0.16055|
| rs3131969|0.08597| 1.067|  A|  G| 0.133028|
| rs3131967|0.06683| 1.077|  T|  C|        .|
| rs1048488| 0.2808|0.9617|  T|  C| 0.836449|
|rs12562034| 0.8489|0.9931|  A|  G|0.0925926|



You can change your GWAS to the same format, or use optional parameter flags to let the script know about the format you are using. More details below.

### .gen file
from [www.shapeit.fr](http://www.shapeit.fr/pages/m02_formats/gensample.html) :
A .gen file is a SPACE delimited file. Each line corresponds to a single SNP. The first 5 columns are:
Chromosome number [integer]
SNP ID [string]
SNP physical position (bp) [integer]
First allele [string]
Second allele [string]

### .vcf file 
This is a default format for the genotype data returned from Sanger Institute. Details about the format can be found [here](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) 

### Output file
By default, the output format of PRS results and SNP logs is csv. 

### SNP log
If the user selects the option, the SNPs that are used for each p-value threshold will be logged into a file with the same name as the output file but with ".snplog" as the extension. This may include the SNPs that have missing calls for some samples. 

## Running command-line script PRS_run.py
### Parameters



A description of the parameters for the script can be obtained by typing: 
```
python PRS_run.py --help
```
Commandline-stype flags are used to specify how the scores are calculated. 

To run the script, use ```spark-submit```. You can add other parameters for spark before the script if desired. 

```
spark-submit PRS_run.py 
```
Followed by three positional parameters 
```

  GENO                  Name of the Genotype files, can be a name or path, or name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  Output                The path and name for the output file
```
Followed by some optional parameters. By default, the pipeline assumes the following. 

* To specify the format of genotype file :
```
  --filetype {GEN,VCF}  
```
The type of genotype file used as input,  choose between VCF and GEN, default is VCF

* To specify the p-value thresholds of the score :
```
  --thresholds  0.5 0.2 0.3
```
Enter one or more float numbers separated by comma. Default is 0.5 0.2 0.1 0.05 0.01 0.001 0.0001

Alternatively you can specify a sequence of thresholds and calculate scores at each threshold:

```
  --threshold_seq THRESHOLD_SEQ 0.1 0.5 0.01
```
After the flag, the first number is the starting point of the sequence, the second is the end point of the sequence, the third number denotes the step size. The above example would yield the sequence 0.1,0.11,0.12.....0.49,0.50. Note the interval is inclusive of the endpoints.

* By default 
                        
  --GWAS_no_header      Adding this parameter signals that there is no headers
                        for the GWAS. The default is to assume that GWAS has
                        column names
  --log_or              Adding this parameter tells the script to log the
                        effect sizes provided in the GWAS
  --check_ref           Adding this option tells the script to theck reference
                        allele when determining genoypte calls. Default is not
                        checking

  --sample_file SAMPLE_FILE
                        path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_delim SAMPLE_DELIM
                        Delimiter of the sample file. Default is comma


  --no_maf              By default, the pipeline calculated the allele
                        frequency in the genotype population. Use this flag to
                        tell the script NOT to calculate MAF in the provided
                        propulation and compare it with MAF in the GWAS, e.g,
                        when the GWAS does not provide information for allele
                        frequencies. MAF is needed to check the reference
                        alleles of ambiguous SNPs (those whose A1 and A2 are
                        reverese complements). Not using this will result in
                        ambiguous SNPs be discarded.
  --snp_log SNP_LOG     Specify the path for a log file that records the SNPs
                        that are used at each threshold. Default is no log




### Examples:
To calculate PRS from a series of .vcf files, while checking the allele allignment between the genotype and the GWAS, and take the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:
```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --sample_file samplefile.csv --sample_file_id 0 --check_ref --log_or --thresholds  0.2 0.1 0.05
```
To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and taking the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:

```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --filetype GEN --sample_file samplefile.csv --sample_file_id 0 --GWAS_no_header --log_or --thresholds  0.2 0.1 0.05
```




### Full list of Parameters when type `python PRS_run.py --help` :
```
positional arguments:
  GENO                  Name of the Genotype files, can be a name or path, or
                        name patterns with wildcard character
  GWAS                  Name of the GWAS file, can be a name or path.
  Output                The path and name stem for the output files. One name
                        will be used for the score output, the snp log and the
                        regression output. This is similar to the --out flag
                        in pLink

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --gwas_id GWAS_ID     Column number in your GWAS that contains SNP ID, with
                        first column being 0, default is 0
  --gwas_p GWAS_P       Column number in your GWAS that contains p-value, with
                        first column being 0, default is 1
  --gwas_or GWAS_OR     Column number in your GWAS that contains odds-
                        ratio/beta, with first column being 0, default is 2
  --gwas_a1 GWAS_A1     Column number in your GWAS that contains allele A1,
                        with first column being 0, default is 3.
  --gwas_a2 GWAS_A2     Column number in your GWAS that contains allele A@,
                        with first column being 0, default is 4.
  --gwas_a1f GWAS_A1F   Column number in your GWAS that contains frequency of
                        A1, with first column being 0, default is 5.
  --filetype {GEN,VCF}  The type of genotype file used as input , choose
                        between VCF and GEN, default is VCF
  --thresholds THRESHOLDS [THRESHOLDS ...]
                        The p-value thresholds that controls which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. default is [0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001]
  --threshold_seq THRESHOLD_SEQ [THRESHOLD_SEQ ...]
                        Defines a sequence that contains all the p-value
                        thresholds that controls which SNPs are used from the
                        GWAS. Input is three numbers separted by space: lower
                        bound, upper bound, step size. Default is None.
                        Defining a sequence automatically overwrites the
                        threshold list defined under --thresholds
  --GWAS_delim GWAS_DELIM
                        Delimtier of the GWAS file, default is tab-delimiter
  --GWAS_no_header      Adding this parameter signals that there is no headers
                        for the GWAS. The default is to assume that GWAS has
                        column names
  --log_or              Adding this parameter tells the script to log the
                        effect sizes provided in the GWAS
  --no_check_ref        Adding this option tells the script to theck reference
                        allele when determining genoypte calls. Default is not
                        checking
  --app_name APP_NAME   Give your spark application a name. Default is PRS.
  --sample_file SAMPLE_FILE
                        path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_file_delim SAMPLE_DELIM
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
  --no_maf              By default, the pipeline calculated the allele
                        frequency in the genotype population. Use this flag to
                        tell the script NOT to calculate MAF in the provided
                        propulation and compare it with MAF in the GWAS, e.g,
                        when the GWAS does not provide information for allele
                        frequencies. MAF is needed to check the reference
                        alleles of ambiguous SNPs (those whose A1 and A2 are
                        reverese complements). Not using this will result in
                        ambiguous SNPs be discarded.
  --snp_log             Specify whether to write the IDs of the SNPs for each
                        scores and how they were processed. If added, the SNP
                        ids will be saved to a file with the name specified in
                        the Output flag, with .snplog as suffix
  --check_dup           Add this flag if you want to check for and discard
                        SNPs that are duplicated, which will take extra time.
                        By default, the script will assume there is no
                        duplicate SNPs.
  --pheno_file PHENO_FILE
                        Sepcify the path to the data file for the phenotype.
                        It is assumed that the phenotype data is organized in
                        the same order as the samples in the genoytpe file.
  --pheno_columns PHENO_COLUMNS [PHENO_COLUMNS ...]
                        Specify which columns that the phenotype data is in
                        the provided phenotype data file. Multiple column
                        numbers can be specified to conduct regression with
                        multiple phenotypes. Default is the first column.
  --pheno_delim PHENO_DELIM
                        Specify the delimiter for the phenotype data file.
                        Default is comma
  --pheno_no_header     Sepcify whether the phenotype has a header row
  --covar_columns COVAR_COLUMNS [COVAR_COLUMNS ...]
                        Specify which columns that the phenotype data is in
                        the provided phenotype data file. Multiple column
                        numbers can be specified to conduct regression with
                        multiple phenotypes. Default is the first column.

``` 
