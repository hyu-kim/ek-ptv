# %%
"""
K-mean algorithm implemented to single cell mobilities having bimodal distribution (one at zero, another at unknown)

Created on Oct 29, 2021

@author: Hyu Kim (hskimm@mit.edu)

"""
import pandas as pd
import numpy as np

def k_means(l, convergence_limit=0.05):
    """
    Filters out zero-centered values and return non-zero mean
    
    Parameters
    ----------
    l: Dataframe
       1-D list filtered with a certain column from tr_av
    convergence_limit: Numeric
                       threshold to stop while loop
    Returns
    -------
    avg_new: Float
             Mean value of the cluster
    l_kmean: Array
             All values filtered out
    """
    arr = np.array(l)
    arr = arr * 1e+9 # change unit
    c = np.zeros(len(arr))
    # convergence_limit = 0.1
    mu0 = 0 # fix this
    avg_old = np.mean(arr) # initial guess
    avg_new = 0
    while abs(avg_old - avg_new)>convergence_limit:
        for i in range(len(c)):
            if abs(arr[i]-avg_old)<abs(arr[i]-mu0):
                c[i] = 1
                # print('updated')
        avg_old = avg_new
        avg_new = sum(np.multiply(c, arr)) / sum(c)
    l_kmean = arr[c==1] * 1e-9

    return avg_new, l_kmean
# # %%
# path_load = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/2020-09-25 Pt mobility/single_cell/'
# s = path_load + 'Ch12_+1R1_R3_single.csv'
# arr = pd.read_csv(s, header=0)
# arr = arr['mobility']
# arr = np.array(arr)
# arr = arr * 1e+9 # change unit
# # %%
# c = np.zeros(len(arr))
# convergence_limit = 0.1
# mu0 = 0 # fix this
# avg_old = np.mean(arr) # initial guess
# mu1_new = 0
# while abs(avg_old - mu1_new)>convergence_limit:
#     for i in range(len(c)):
#         if abs(arr[i]-avg_old)<abs(arr[i]-mu0):
#             c[i] = 1
#             # print('updated')
#     avg_old = mu1_new
#     mu1_new = sum(np.multiply(c, arr)) / sum(c)
