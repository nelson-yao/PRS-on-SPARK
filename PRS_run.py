
 # -*- coding: utf-8 -*-
from __future__ import division


from operator import add
from math import log
import csv
import pickle
import sys
from collections import Counter
import re
import glob, os

import ntpath
import functools
import itertools

from time import time
import argparse


# Takes a line in the genotype file and return the frequency of A1 allele
def getA1f(geno):
    AA=geno[0::3]
    AB=geno[1::3]
    AA2=[x*2 for x in AA]
    A1count=map(add, AA2, AB)
    A1F=sum(A1count)/(float(len(AA2))*2)
    return A1F

# Takes an array of triplets and convert to number of A1 alleles
def makeGenotype(line):
    AA=line[0::3]
    AB= line[1::3]
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return  genotype

# Takes an array of triplets and convert to number of A1 alleles, while checking the strand alignments
def makeGenotypeCheckRef(line, checkMap, toDF=False):
    rsid=line[0]
    gen=line[1]
    if rsid in checkMap:
        if checkMap[rsid]=="keep":
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))

        elif checkMap[rsid]=="flip":
            AA=gen[2::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))


    else:
        # fail-safe, in case a SNP does not have a flag but still exist in the genotype file
        print("SNP {} was not accounted for in the alignment checking step, discarding this SNP".format(rsid))
        #genotype=[0.0]*(len(gen)/3)
    if genotype:
        if toDF:  # future work, in case want to use SQL methods on Spark DataFrame
            return [rsid]+genotype
        else:
            return (rsid, genotype)


def filterGWASByP(GWASRdd, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASRdd.filter(lambda line: (float(eval(line[pcolumn]))<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair


def filterGWASByP_DF(GWASdf, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASdf.rdd.filter(lambda line: (float(line[pcolumn])<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))

    else:

        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))  ## make tuples of the format (snpid, effect size)
    GWASoddspair=GWAS_Odds.collectAsMap()   # make a python dictionary of the format { snpid: effect_size }
    return GWASoddspair

## The following function is for checking alignment in genotype RDDs
def checkAlignment(line):
    bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
    genoA1=line[0][0][0]
    genoA2=line[0][0][1]
    genoA1F=line[0][1]
    gwasA1=line[1][0][0]
    gwasA2=line[1][0][1]
    gwasA1F=line[1][1]

    try:
        if genoA2==bpMap[genoA1]:
            if gwasA1F==".":
                flag="discard"
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5

                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff*genoA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"

    finally:
        return flag


# checking for reference allele alignment (i.e is the GWAS A1 the same as A1 in the genotype)
def checkAlignmentDF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1]
    genoA2=dataframe[2]
    genoA1F=dataframe[3]
    # no number 4 because the 5th column has the snpid from the GWAS
    gwasA1=dataframe[5]
    gwasA2=dataframe[6]
    gwasA1F=dataframe[7]
    if genoA1 in bpMap:
        if genoA1==genoA2 or gwasA1 == gwasA2:
            flag='discard'
        elif genoA2==bpMap[genoA1]:  # checking ambiguous SNP in the genotype
            if gwasA1F==".":
                flag="discard"
            else:

                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5
                # if the allele frequency in the GWAS or genotype is between 0.4 and 0.6, too close to decide, discard the SNP
                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff * gwasA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        #example if geno has (A,G) and GWAS has (G,T), then need to discard
        elif (genoA2==gwasA1 and genoA1==gwasA2) or (genoA2==bpMap[gwasA1] and genoA1==bpMap[gwasA2]) :
            flag="flip"

        elif (genoA1==gwasA1 and genoA2==gwasA2) or (genoA1==bpMap[gwasA1] and genoA2==bpMap[gwasA2]) :
            flag="keep"

        else:
            flag="discard"

    else :
        flag="discard"
        print("Invalid Genotypes for SNP {}; genotype Alleles:{},{}; GWAS alleles: {},{}".format(snpid,genoA1, genoA2, gwasA1, gwasA2))

    return (snpid,flag)

# check reference allele alignment without using A1F in the GWAS. Ambiguous SNPS are automatically discarded
def checkAlignmentDFnoMAF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1]
    genoA2=dataframe[2]
    gwasA1=dataframe[4]
    gwasA2=dataframe[5]

    if genoA1 in bpMap:
        if genoA2==bpMap[genoA1] or genoA1==genoA2:  # discard ambiguous case and the case where A1 = A2 in the genotype data
            flag="discard"
        elif (genoA2==gwasA1 and genoA1==gwasA2) or (genoA2==bpMap[gwasA1] and genoA1==bpMap[gwasA2]) :
            flag="flip"
        elif (genoA1==gwasA1 and genoA2==gwasA2) or (genoA1==bpMap[gwasA1] and genoA2==bpMap[gwasA2]) :
            flag="keep"
        else:
            flag='discard'

    else:
        flag="discard"
        print("Invalid Genotypes for SNP {}; Genotype alleles:{},{}; GWAS alleles: {},{}".format(snpid,genoA1, genoA2, gwasA1, gwasA2))

    return (snpid,flag)


def getSampleNames(sampleFileName, sampleDelim, sampleIDCol, skip=0):
    labels=[]
    with open(sampleFileName, "r") as f:
        subjList=[item.split(sampleDelim) for item in f.read().splitlines()]
        for i in range(len(sampleIDCol)):
            subjNames=[x[i] for x in subjList[skip::]]
            subjNames=[name.strip('"') for name in subjNames]
            column=["Label"+str(i)]+subjNames
            labels.append(column)
    return labels

# remove duplicates from flag list
def rmDup(checkList):
    checkCount=Counter([x[0] for x in checkList])
    nodupMap={snp:flag for snp, flag in checkList if checkCount[snp]==1}
    return nodupMap

# output each PRS for each sample, in the form of [sample, *scores],
# and a list of pvalues that are in the order of the scores from each p-value

def writePRS(prsResults, outputFile, samplenames=None, dialect=None):
    onescore=list(prsResults.values())[0][1]
    samplesize=len(onescore)
    if not samplenames:
        print("No sample names provided, generating sample names")
        samplenames=["Label"]+["Sample"+str(i+1) for i in range(samplesize)]
    labelnumber=len(samplenames[0])-1
    print("Collected {} sample labels".format(labelnumber))
    if labelnumber==samplesize:
        outputdata=samplenames
        pvaluelist=sorted(list(prsResults.keys()))
        for pvalue in pvaluelist:
            outputdata.append(["SNP_count_{}".format(pvalue)]+[prsResults[pvalue][0]]*samplesize)
            outputdata.append(["PRS_{}".format(pvalue)]+prsResults[pvalue][1])

        try:
            with open(outputFile, "w") as f:
                csvwriter=csv.writer(f, dialect=dialect)
                for row in zip(*outputdata):
                    csvwriter.writerow(row)
                print("Successfully wrote scores to "+ os.path.basename(outputFile))
        except:
            e = sys.exc_info()[0]
            print( "<p>Error: %s</p>" % e )
            print("Data output was unsuccessful.")
            print("All is not lost, final results saved as binary format in file 'PRSOutput.pk'")
            with open(os.path.dirname(outputFile)+"/PRSOutput.pk", "wb") as f:
                pickle.dump(outputdata, f)
    else:
        print("Unequal number of labels extracted from sample sheet and number of samples detected in the genotype data, saving results to PRSOutput.pk")
        outputdata=False
        with open(os.path.join(os.path.dirname(outputFile),"PRSOutput.pk"), "wb") as f:
            pickle.dump({"Samples": samplenames, "results": prsResults}, f)

    return outputdata

def writeSNPlog(snpid, outputFile, dialect=None):
    outputdata=[]
    maxT=max(snpid.keys())
    maxLen=len(snpid[maxT])
    for pvalue in snpid.keys():
        outputdata.append([str(pvalue)]+sorted(snpid[pvalue])+[""]*(maxLen-len(snpid[pvalue])))
    try:
        with open(outputFile, "w") as f:
            csvwriter=csv.writer(f, dialect=dialect)
            for row in zip(*outputdata):
                csvwriter.writerow(row)
            print("Successfully output log to "+ os.path.basename(outputFile))
    except:
        e = sys.exc_info()[0]
        print( "<p>Error: %s</p>" % e )
        print("Data output was unsuccessful.")
        print("All is not lost, final results saved as binary format in file 'SNPlog.pk'")
        with open(os.path.join(os.path.dirname(outputFile),"SNPlog.pk"), "wb") as f:
            pickle.dump(outputdata, f)
    return outputdata



if __name__=="__main__":
	import pyspark
	from pyspark.sql import SparkSession
	from pyspark.sql import Row
	from pyspark.sql.functions import udf
	from pyspark.sql.types import *
	from pyspark import SparkConf, SparkContext
    # ATTN: python index starts at 0, so if you want to specify the second column, use 1
    # define column number for contents in GWAS
	parser = argparse.ArgumentParser(description='PRS Script Parameters')
	# Mandatory positional arguments
	parser.add_argument("GENO", action="store", help="Name of the Genotype files, can be a name or path, or name patterns with wildcard character ")
	parser.add_argument("GWAS", action="store", help="Name of the GWAS file, can be a name or path.")
	parser.add_argument("Output", action="store", help="The path and name for the output file")
    # Optional arguments
	parser.add_argument("--gwas_id", action="store", default=0, dest="gwas_id",type=int, help="Column number in your GWAS that contains SNP ID, with first column being 0, default is 0")
	parser.add_argument("--gwas_p", action="store", default=1, dest="gwas_p", type=int, help="Column number in your GWAS that contains p-value, with first column being 0, default is 1")
	parser.add_argument("--gwas_or", action="store", default=2, dest="gwas_or", type=int, help="Column number in your GWAS that contains odds-ratio/beta, with first column being 0, default is 2")
	parser.add_argument("--gwas_a1", action="store", default=3, dest="gwas_a1", type=int, help="Column number in your GWAS that contains allele A1, with first column being 0, default is 3. Allele A2 is assumed to be at column [gwas_a1+1]")
	parser.add_argument("--gwas_a1f", action="store", default=5, dest="gwas_a1f", type=int, help="Column number in your GWAS that contains frequency of A1, with first column being 0, default is 5.")

    #parser.add_argument("--geno_id", action="store", default=2, dest="geno_id",type=int, help="Column number in your genotype files that contains SNP ID, with first column being 0, default is 2")
    #parser.add_argument("--geno_start", action="store", default=9,dest="geno_start", type=int, help="Column number in your genotype files that contains the first genotype,  with first column being 0, default is 9.")
    #parser.add_argument("--geno_a1", action="store", default=10, dest="geno_a1",type=int, help="Column number in your genotype files that contains SNP ID, with first column being 0, default is 10.")
    #parser.add_argument("--GENO_delim", action="store", default="\t", dest="GENO_delim", help="Delimtier of the GWAS file, default is tab-delimiter ")


	parser.add_argument("--filetype", action="store",default="VCF", dest="filetype", help="The type of genotype file used as input , choose between VCF and GEN, default is VCF", choices=set(["VCF", "GEN"]))

	parser.add_argument("--thresholds", action="store", default=[0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001], dest="thresholds", help="The p-value thresholds that controls which SNPs are used from the GWAS. Specifying the p-values simply by input one after another. default is [0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001]", nargs="+", type=float)

	parser.add_argument("--GWAS_delim", action="store", default="\t", dest="GWAS_delim", help="Delimtier of the GWAS file, default is tab-delimiter ")

	parser.add_argument("--GWAS_no_header", action="store_false", default=True, dest="GWAS_header", help="Adding this parameter signals that there is no headers for the GWAS. The default is to assume that GWAS has column names")

	parser.add_argument("--log_or", action="store_true", default=False, dest="log_or", help="Adding this parameter tells the script to log the effect sizes provided in the GWAS")

	parser.add_argument("--check_ref", action="store_true", default=False, dest="check_ref", help="Adding this option tells the script to theck reference allele when determining genoypte calls. Default is not checking")

	parser.add_argument("--app_name", action="store", default="PRS", dest="app_name", help="Give your spark application a name. Default is PRS.")

	parser.add_argument("--sample_file", action="store", dest="sample_file", default="NOSAMPLE",help="path and name of the file that contain the sample labels. It is assumed that the sample labels are already in the same order as in the genotype file.")

	parser.add_argument("--sample_file_delim", action="store", default=",", dest="sample_delim", help="Delimiter of the sample file. Default is comma")

	parser.add_argument("--sample_file_ID", action="store", default=[0], type=int, nargs="+", dest="sample_file_ID", help="Specify which columns in the sample file are used as labels. Can use one integer to specify one column, or multiple integers to specify multiple columns. Default is the first column")

	parser.add_argument("--sample_file_skip", action="store",default=1, dest="sample_skip", help="Specify how many lines to skip in the sample file, i.e. which row do the labels start. Default is 1, which assumes that the sample files has column names and the labels start on the second line", type=int)

	parser.add_argument("--no_maf", action="store_false", default=True, dest="use_maf", help="By default, the pipeline calculated the allele frequency in the genotype population. Use this flag to tell the script NOT to calculate MAF in the provided propulation and compare it with MAF in the GWAS, e.g, when the GWAS does not provide information for allele frequencies. MAF is needed to check the reference alleles of ambiguous SNPs (those whose A1 and A2 are reverese complements).  Not using this will result in ambiguous SNPs be discarded.")

	parser.add_argument("--snp_log", action="store", default=None, dest="snp_log", help="Specify the path for a log file that records the SNPs that are used at each threshold. Default is no log")

	parser.add_argument("--check_dup", action="store_true", default=False, dest="checkdup", help="Add this flag if you want to check for and discard SNPs that are duplicated, which will take extra time. By default, the script will assume there is no duplicate SNPs. You can use clean.py to remove duplicated SNPs in your data ")

	results=parser.parse_args()

    # type of files, VCF or GEN
	filetype=results.filetype

    # log file
	snp_log=results.snp_log
    ## Setting parameters
	gwas_id=results.gwas_id    # column of SNP ID
	gwas_p=results.gwas_p     # column of P value
	gwas_or=results.gwas_or    # column of odds ratio
	gwas_a1=results.gwas_a1    # column of a1 in the GWAS
	gwas_a1f=results.gwas_a1f  # column index of maf in the GWAS

    # defin column number for contents in genfile
	if filetype.lower()=="vcf":
		geno_id= 2 # column number with rsID
		geno_start=9 # column number of the 1st genotype, in the raw vcf files, after separated by the delimiter of choice
		geno_a1 = 3 # column number that contains the reference allele
		GENO_delim= "\t"
	elif filetype.lower()=="gen":
		geno_id = 1
		geno_start=5
		geno_a1=3
		GENO_delim= " "
    # List of thresholds:
	thresholds=results.thresholds
    # file delimiters:
	GWAS_delim=results.GWAS_delim

    # file names:
    #home="/Volumes/mavan/Genotyping_161114/MAVAN_imputed_161121/KIDS_info03/"  #define homefolder path

    # Name of GWAS file
	gwasFiles=results.GWAS
	GWAS_has_header=results.GWAS_header

    # programme parameter
	log_or=results.log_or  # sepcify whether you want to log your odds ratios
	check_ref=results.check_ref # if you know that there are mismatch between the top strand in the genotypes and that of the GWAS, set True. Not checking the reference allele will improve the speed
	use_maf=results.use_maf   # wheather to use MAF to check reference allele

    # sample file path and name
	sampleFilePath=results.sample_file # include the full/relative path and name of the sample file
	sampleFileDelim=results.sample_delim  # sample File Delimiter
	sampleFileID=results.sample_file_ID   # which column in the sample file has the ID
	sample_skip=results.sample_skip  # how many lines to skip so that the sample names can be matched to the genotypes 1-to-1, taking into account the header of the sample file
    ##output file information

	outputPath=results.Output

    # Sepcify whether to check for duplicate SNPs
	checkDup= results.checkdup


    # get the name of the genotype files
	genoFileNamePattern=results.GENO

    # get the whole list of the file names
	genoFileNames=glob.glob(genoFileNamePattern)

    ##  start spark context
	APP_NAME=results.app_name

	spark=SparkSession.builder.appName(APP_NAME).getOrCreate()

    # if using spark < 2.0.0, use the pyspark module to make Spark context
    # conf = pyspark.SparkConf().setAppName(APP_NAME).set()#.set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")

	sc  = spark.sparkContext

    #sc = spark.sparkContext
	sc.setLogLevel("WARN")
	log4jLogger = sc._jvm.org.apache.log4j
	LOGGER = log4jLogger.LogManager.getLogger(__name__)
	LOGGER.warn("Start Reading Files")
	LOGGER.warn("Using these genoytpe files: ")

	for filename in genoFileNames[:min(24, len(genoFileNames))]:
		LOGGER.warn(filename)
	if len(genoFileNames)>23:
		LOGGER.warn("and more...")

	LOGGER.warn("total of {} files".format(str(len(genoFileNames))))
	# 1. Load files

    # read the raw data
	genodata=sc.textFile(genoFileNamePattern)
	#LOGGER.warn("Using the GWAS file: {}".format(ntpath.basename(gwasFiles)))
	LOGGER.warn("Using the GWAS file: {}".format(gwasFiles))
	gwastable=spark.read.option("header",GWAS_has_header).option("delimiter",GWAS_delim).csv(gwasFiles).cache()
	print("Showing top 5 rows of GWAS file")
	gwastable.show(5)

	LOGGER.warn("System recognizes the following information in the GWAS :")
	LOGGER.warn("SNP ID : Column {}".format(gwas_id))
	LOGGER.warn("P-values : Column {}".format(gwas_p))
	LOGGER.warn("Effect size : Column {}".format(gwas_or))
	LOGGER.warn("Allele A1 : Column {}".format(gwas_a1))
	LOGGER.warn("Allele A2 : Column {}".format(gwas_a1+1))
	if use_maf:
		LOGGER.warn("Allele Frequencies : Column {}".format(gwas_a1f))

    # 1.1 Filter GWAS and prepare odds ratio

    # filter the genotype to contain only the SNPs less than the maximum p value threshold in the GWAS
	maxThreshold=max(thresholds)  # maximum p value
	gwasOddsMapMax=filterGWASByP_DF(GWASdf=gwastable, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=maxThreshold, logOdds=log_or)
	gwasOddsMapMaxCA=sc.broadcast(gwasOddsMapMax).value  # Broadcast the map

    # ### 2. Initial processing
    # at this step, the genotypes are already filtered to keep only the ones in 'gwasOddsMapMax'
	bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
	tic=time()
	if filetype.lower()=="vcf":
	    LOGGER.warn("Genotype data format : VCF ")

        # [chrom, bp, snpid, A1, A2, *genotype]
	    genointermediate=genodata.filter(lambda line: ("#" not in line)).map(lambda line: line.split(GENO_delim)).filter(lambda line: line[geno_id] in gwasOddsMapMaxCA).map(lambda line: line[0:5]+[chunk.split(":")[3] for chunk in line[geno_start::]]).map(lambda line: line[0:5]+[triplet.split(",") for triplet in line[5::]])

        ## (snpid, [genotypes])
        genotable=genointermediate.map(lambda line: (line[geno_id], list(itertools.chain.from_iterable(line[5::])))).mapValues(lambda geno: [float(x) for x in geno])
        if check_ref:
	        if use_maf:
	            LOGGER.warn("Correcting strand alignment, using MAF")
	            genoA1f=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1]), [float(x) for x in list(itertools.chain.from_iterable(line[5::]))])).map(lambda line: (line[0], line[1][0], line[1][1], getA1f(line[2]))).toDF(["Snpid_geno", "GenoA1", "GenoA2", "GenoA1f"])

                # 'GwasA1F' means the allele of the A1 frequency in the GWAS
	            gwasA1f=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a1+1], line[gwas_a1f])).toDF(["Snpid_gwas", "GwasA1", "GwasA2", "GwasA1F"])

                # checktable = [ geno_snpid, genoA1, genoA2, genoA1f, gwas_snpid, gwasA1, gwasA2, gwasA1f]
	            checktable=genoA1f.join(gwasA1f, genoA1f["Snpid_geno"]==gwasA1f["Snpid_gwas"], "inner").cache()
	            if checkDup:
	                flagList = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collect()  #  (snpid, flag)
	                flagMap = rmDup(flagList)
	            else:
	                flagMap = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collectAsMap()

	        else:
	            LOGGER.warn("Correcting strand alignment, without using MAF")
	            genoalleles=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1]), [float(x) for x in list(itertools.chain.from_iterable(line[5::]))])).map(lambda line: (line[0], line[1][0], line[1][1])).toDF(["Snpid_geno", "GenoA1", "GenoA2"])

	            gwasalleles=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a1+1])).toDF(["Snpid_gwas", "GwasA1", "GwasA2"])

	            checktable=genoalleles.join(gwasalleles, genoalleles["Snpid_geno"]==gwasalleles["Snpid_gwas"], "inner").cache()

	            if checkDup:
	                flagList = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collect()
	                flagMap = rmDup(flagList)
	            else:
                    # no need to check the duplicates if the data is preprocessed
	                flagMap = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collectAsMap()

	        LOGGER.warn("Generating genotype dosage while taking into account difference in strand alignment")
	        flagMap=sc.broadcast(flagMap).value
	        genotypeMax=genotable.filter(lambda line: line[0] in flagMap and flagMap[line[0]]!="discard").map(lambda line: makeGenotypeCheckRef(line, checkMap=flagMap)).cache()

        else:
	        LOGGER.warn("Generating genotype dosage without checking reference allele alignments")
	        genotypeMax=genotable.mapValues(lambda line: makeGenotype(line)).cache()
	        if checkDup:
	            genotypeCount=genotypeMax.map(lambda line: (line[0], 1)).reduceByKey(lambda a,b: a+b).filter(lambda line: line[1]==1).collectAsMap()
	            genotypeMax=genotypeMax.filter(lambda line: line[0] in genotypeCount)


	elif filetype.lower()=="gen":
	    LOGGER.warn("Genotype data format : GEN ")
	    genotable=genodata.map(lambda line: line.split(GENO_delim)).filter(lambda line: line[geno_id] in gwasOddsMapMaxCA).map(lambda line: (line[geno_id], line[geno_start::])).mapValues(lambda geno: [float(call) for call in geno])

	    if check_ref:
	        if use_maf:
	            LOGGER.warn("Correcting strand alignment, using MAF")
	            genoA1f=genodata.map(lambda line: line.split(GENO_delim)).map(lambda line: (line[geno_id], line[geno_a1], line[geno_a1+1], getA1f([float(x) for x in line[geno_start::]]))).toDF(["Snpid_geno", "GenoA1", "GenoA2", "GenoA1f"])

	            gwasA1f=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a1+1], line[gwas_a1f])).toDF(["Snpid_gwas", "GwasA1", "GwasA2", "GwasA1f"])
	            checktable=genoA1f.join(gwasA1f, genoA1f["Snpid_geno"]==gwasA1f["Snpid_gwas"], "inner").cache()
	            if checkDup:
	                LOGGER.warn("Searching and removing duplicated SNPs")
	                flagList = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collect()
	                flagMap = rmDup(flagList)
	            else:
	                flagMap = checktable.rdd.map(lambda line: checkAlignmentDF(line, bpMap)).collectAsMap()
	        else:
	            LOGGER.warn("Correcting strand alignment, without using MAF")
	            genoalleles=genodata.map(lambda line: line.split(GENO_delim)).map(lambda line: (line[geno_id], line[geno_a1], line[geno_a1+1])).toDF(["Snpid_geno", "GenoA1", "GenoA2"])
	            gwasalleles=gwastable.rdd.map(lambda line:(line[gwas_id], line[gwas_a1], line[gwas_a1+1])).toDF(["Snpid_gwas", "GwasA1", "GwasA2"])
	            checktable=genoalleles.join(gwasalleles, genoalleles["Snpid_geno"]==gwasalleles["Snpid_gwas"], "inner").cache()

	            if checkDup:
	                LOGGER.warn("Searching and removing duplicated SNPs")
	                flagList = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collect()
	                flagMap = rmDup(flagList)
	            else:
	                flagMap = checktable.rdd.map(lambda line: checkAlignmentDFnoMAF(line, bpMap)).collectAsMap()

	        LOGGER.warn("Generating genotype dosage while taking into account difference in strand alignment")
	        flagMap=sc.broadcast(flagMap).value
	        genotypeMax=genotable.filter(lambda line: line[0] in flagMap).map(lambda line: makeGenotypeCheckRef(line, checkMap=flagMap)).cache()

	    else:
	        LOGGER.warn("Generating genotype dosage without checking strand alignments")
	        genotypeMax=genotable.mapValues(lambda line: makeGenotype(line)).cache()
	        if checkDup:
	            genotypeCount=genotypeMax.map(lambda line: (line[0], 1)).reduceByKey(lambda a,b: a+b).filter(lambda line: line[1]==1).collectAsMap()
	            genotypeMax=genotypeMax.filter(lambda line: line[0] in genotypeCount)

	LOGGER.warn("Dosage generated in {:f} seconds".format(time()-tic) )
	samplesize=int(len(genotypeMax.first()[1]))
	LOGGER.warn("Detected {} samples" .format(str(samplesize)))

	#genoa1f.map(lambda line:"\t".join([line[0], "\t".join(line[1]), str(line[2])])).saveAsTextFile("../MOMS_info03_maf")

	# Calculate PRS at the sepcified thresholds

	def calcPRSFromGeno(genotypeRDD, oddsMap):
	    totalcount=genotypeRDD.count()
	    multiplied=genotypeRDD.map(lambda line:[call * oddsMap[line[0]] for call in line[1]])
	    PRS=multiplied.reduce(lambda a,b: map(add, a, b))
	    normalizedPRS=[x/totalcount for x in PRS]
	    return (totalcount,PRS)

	def calcAll(genotypeRDD, gwasRDD, thresholdlist, logsnp):
	    prsMap={}
        thresholdNoMaxSorted=sorted(thresholdlist, reverse=True)

	    thresholdmax=max(thresholdlist)
	    idlog={}
	    start=time()
	    for threshold in thresholdNoMaxSorted:
	        tic=time()
	        gwasFilteredBC=sc.broadcast(filterGWASByP_DF(GWASdf=gwasRDD, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=threshold, logOdds=log_or))
	        #gwasFiltered=spark.sql("SELECT snpid, gwas_or_float FROM gwastable WHERE gwas_p_float < {:f}".format(threshold)
	        LOGGER.warn("Filtered GWAS at threshold of {}. Time spent : {:f} seconds".format(str(threshold), time()-tic))
	        checkpoint=time()
	        filteredgenotype=genotypeRDD.filter(lambda line: line[0] in gwasFilteredBC.value)

	        if not filteredgenotype.isEmpty():
	            if logsnp:
	                idlog[threshold]=filteredgenotype.map(lambda line:line[0]).collect()
	            prsMap[threshold]=calcPRSFromGeno(filteredgenotype, gwasFilteredBC.value)

	            LOGGER.warn("Finished calculating PRS at threshold of {}. Time spent : {:f} seconds".format(str(threshold), time()-checkpoint))
	    return prsMap, idlog

	prsDict, snpids=calcAll(genotypeMax,gwastable, thresholds, logsnp=snp_log)

    # log which SNPs are used in PRS
	if snp_log:
	    logoutput=writeSNPlog(snpids, snp_log)
    # generate labels for samples
    if filetype.lower()=="vcf":
        subjNames=genodata.filter(lambda line: "#CHROM" in line).map(lambda line: line.split(GENO_delim)[9::])
	elif sampleFilePath!="NOSAMPLE":
        # get sample name from the provided sample file
        subjNames=getSampleNames(sampleFilePath,sampleFileDelim,sampleFileID, skip=sample_skip)
	    LOGGER.warn("Extracted {} sample labels".format(len(subjNames)))
	    output=writePRS(prsDict,  outputPath, samplenames=subjNames)
	else:
       LOGGER.warn("No sample file detected, generating labels for samples.")
       output=writePRS(prsDict,  outputPath, samplenames=None)


	sc.stop()
