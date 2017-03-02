import numpy as np
import csv









def get_snplog(path,pvalues):

    #pvalues = np.arange(range[0],range[1],range[2]).tolist()
    
    t = [[] for i in range(len(pvalues))]
    with open(path,'r') as snplog:
        header = next(snplog)
        for snp in snplog:
            snp = snp.split(',')
            for i in range(len(t)):
                if snp[i*2]!='':
                    t[i].append(snp[i*2])
    print t[0]
    return t



  
path_to_gwas=  'GWAS_SCZ_beta_noaf_full_noambi_NEURON_KIDS_noamb_hc09_clump.txt'
  
def get_gwas(path_to_gwas):
    g_dict = {}      
    f=open(path_to_gwas,'r')
   
    gwas = csv.reader(f,delimiter='\t')
    header = next(gwas)
    for l in gwas:
        g_dict[l[0]] = l[1:]
    return header,g_dict


#output the snps that are considered to be adding noise to prs
def DRSFiltering(rs,snplog,gwas):
    re = []
    for i in range(len(snplog)-1).reverse():
        if rs[i+1] < rs[i]:
          re = re +  list(set(snplog[i+1])-set(snplog[i]))
    for s in re:
        del gwas[s]
    
    return gwas

if __name__ == '__main__':
    path_to_gwas=  'GWAS_SCZ_beta_noaf_full_noambi_NEURON_KIDS_noamb_hc09_clump.txt'

    range = [0.001, 0.901, 0.001]
    gwas = get_gwas(path_to_gwas)
    
    print gwas['rs55777713']