"""
This file holds the plotting components for the PRS_on_Spark
@author: zhuyuecai
@note: the main function is only for testing purpose
@requires: matplotlib, for testing, you need nupmy
"""

import matplotlib.pyplot as plt
from astropy.units import count



"""
@author: zhuyuecai
@param pheno_rs: a dictionary with the pheno type name as keys, and a list of ordered R square as values
@param p_values: a list of ordered p values range from 0 to 1  
@param width: number of subplot in one row, default 2
"""

def r_square_plots(pheno_rs, p_values,width = 2):
    n_plot = len(pheno_rs)
    fig = plt.figure()
    axs = []
    for i in range(n_plot):
        axs.append(fig.add_subplot(n_plot/2,width, i+1))
    counter = 0
    anno = '%1.2f'
    for k,v in pheno_rs.iteritems():
        ms = _get_max_positions_(v)

        axs[counter].plot(p_values,v)
        for m in ms:
            axs[counter].plot([p_values[m],p_values[m]],[v[m],0])
            axs[counter].annotate(anno%(p_values[m]),[p_values[m],v[m]])
            
            
        axs[counter].set_title(k)
        counter+=1

    plt.show()
        
        


    
    
    
    
    
    



    
"""
@author: zhuyuecai
@note: helper function, used in other functions
@param rs: a list of numerical values
@return: a list of the positions of the maximums in rs
"""
    

def _get_max_positions_(rs):
    m = max(rs)
    return [i for i, j in enumerate(rs) if j == m]



if __name__ == '__main__':
    import numpy as np
    prs = {'t1': [1,2,5,4,1,8,2],'t2': [6,7,5,4,1,8,2]}
    ps = np.linspace(0,1,7)
    r_square_plots(prs,ps)
    
    
    