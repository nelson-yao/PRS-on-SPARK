"""
This file holds the plotting components for the PRS_on_Spark
@author: zhuyuecai
@note: the main function is only for testing purpose
@requires: matplotlib, for testing, you need nupmy
"""

import matplotlib.pyplot as plt
from astropy.units import count
from matplotlib import cm
import matplotlib as mpl
import math
"""
@author: zhuyuecai
@param pheno_rs: a dictionary with the pheno type name as keys, and a list of ordered R square as values
@param p_values: a list of ordered p values range from 0 to 1  
@param width: number of subplot in one row, default 2
"""

def r_square_plots(pheno,rs,p_for_rs, p_values,width = 2,bar_width = 0.01):
    n_plot = len(pheno)
    fig = plt.figure()
    axs = []
    print n_plot/width+1
    for i in range(n_plot):
        axs.append(fig.add_subplot(math.ceil(float(n_plot)/width+1),width, (i%width)+1))
    axs.append(fig.add_subplot(n_plot/width+1,1,n_plot%width))
    counter = 0
    anno = '%1.2f'
    for i in range(n_plot):
        ms = _get_max_positions_(rs[i])

        #axs[counter].bar(p_values,rs[i],(max(p_values)-p_values[0])/len(p_values),color =  cm.get_cmap('cool')(p_for_rs[i]))
        axs[counter].bar(p_values,rs[i],bar_width,color =  cm.get_cmap('cool')(p_for_rs[i]))
        
        for m in ms:
            #axs[counter].scatter(p_values[m],rs[i][m]+bar_width,color = 'red',marker = "d")
            axs[counter].annotate(anno%(p_values[m]),[p_values[m],rs[i][m]])
            
            
        axs[counter].set_title(pheno[i])
        axs[counter].set_ylim([min(rs[i]),max(rs[i])])
        counter+=1
    #axs[n_plot]
    
    cmap = mpl.cm.cool
    norm = mpl.colors.Normalize(vmin=0, vmax=1.0)
    cb1 = mpl.colorbar.ColorbarBase(axs[counter], cmap=cmap,
                                norm=norm,orientation='horizontal')
    #plt.colorbar()
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
    pheno = ['t1','t2','t3']
    prs = [[1,2,5,4,1,8,2],[6,7,5,4,1,8,2],[6,7,5,4,1,8,2]]
    p_rs = [np.linspace(0,1,7),np.linspace(0,1,7),np.linspace(0,1,7)]
    ps = np.linspace(0,1,7)
    r_square_plots(pheno,prs,p_rs,ps)
    
    
    