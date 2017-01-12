# PRS-on-SPARK
Polygenic risk score pipeline for large genotype data, generally imputed data

## Installation

To clone the repository, use 
```
git clone https://github.com/seriousNel/PRS-on-SPARK.git
```

## Software Requirements

The notebooks and scripts require the following to run :
+ spark-2.0.0 +
+ Pyhon 2.7

Instruction of installing Apache Spark can be found [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/)


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

## Running command-line script PRS_run.py
### Parameters

There are three required positional parameters: 
1. path to the genotype file(s)
2. path to the GWAS file
3. path to the output file


A description of the parameters for the script can be obtained by typing: 
```
python PRS_run.py
```
Some important parameters are:

```

positional arguments:
  GENO                  Name of the Genotype files, can be a name or path, or
                        name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  Output                The path and name for the output file

optional arguments:
 
  --filetype {GEN,VCF}  The type of genotype file used as inputm choose
                        between VCF and GEN, default is VCF
  --thresholds THRESHOLDS [THRESHOLDS ...]
                        The p-value thresholds that controls which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. default is [0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001]
                        
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


```

### Examples:
To calculate PRS from a series of .vcf files, while checking the allele allignment between the genotype and the GWAS, and take the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:
```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --sample_file samplefile.csv --sample_file_id 0 --check_ref --log_or --thresholds  0.2 0.1 0.05
```
To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and taking the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:

```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --filetype GEN --sample_file samplefile.csv --sample_file_id 0 --GWAS_no_header --log_or --thresholds  0.2 0.1 0.05
```

