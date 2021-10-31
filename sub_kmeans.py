# %%
"""
K-mean algorithm implemented to single cell mobilities having bimodal distribution (one at zero, another at unknown)

Created on Oct 29, 2021

@author: Hyu Kim (hskimm@mit.edu)

"""
import pandas as pd
import numpy as np
# %%
def k_means(tr_av_mobility, convergence_limit=0.05):
    """
    Filters out zero-centered values and return non-zero mean
    
    Parameters
    ----------
    tr_av_mobility: Dataframe
           trace matrix filtered with mobility value from tr_av
    convergence_limit: Numeric
                       threshold to stop while loop
    Returns
    -------
    mobility_av: Float
                 Mean value of mobility cluster
    tr_av_kmean: Array
                 trace matrix with mobilites filtered out
    """
    arr = np.array(tr_av_mobility)
    arr = arr * 1e+9 # change unit
    c = np.zeros(len(arr))
    # convergence_limit = 0.1
    mu0 = 0 # fix this
    mu1_old = np.mean(arr) # initial guess
    mu1_new = 0
    while abs(mu1_old - mu1_new)>convergence_limit:
        for i in range(len(c)):
            if abs(arr[i]-mu1_old)<abs(arr[i]-mu0):
                c[i] = 1
                # print('updated')
        mu1_old = mu1_new
        mu1_new = sum(np.multiply(c, arr)) / sum(c)
    tr_av_kmean = arr[c==1] * 1e-9

    return mu1_new, tr_av_kmean
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
# mu1_old = np.mean(arr) # initial guess
# mu1_new = 0
# while abs(mu1_old - mu1_new)>convergence_limit:
#     for i in range(len(c)):
#         if abs(arr[i]-mu1_old)<abs(arr[i]-mu0):
#             c[i] = 1
#             # print('updated')
#     mu1_old = mu1_new
#     mu1_new = sum(np.multiply(c, arr)) / sum(c)
