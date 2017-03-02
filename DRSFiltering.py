import numpy as np

pvalues = np.arange(0.001, 0.901, 0.001).tolist()

t = [[] for i in range(len(pvalues))]
with open('susceptibility.snplog','r') as snplog:
    header = next(snplog)
    for snp in snplog:
        snp = snp.split(',')
        for i in range(len(t)):
            if snp[i*2]!='':
                t[i].append(snp[i*2])
print len(t[0])
print len(t[1])


        

        


