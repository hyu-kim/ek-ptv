# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Nov 30 2021

For LPS-DEP project. Updates --
Oct 30: 1) file path, 2) noise reduction for bacterial cell tracking
Nov 30: Added a cell that reads a list of files in directory to create 'info.txt'

@author: Hyu Kim (hskimm@mit.edu)
"""
from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sub,sub2, pims
import trackpy as tp
import time

# from trackpy.sub import _
# %%
"""
Reads a list of files in a directory and creates a draft 'info.txt'
format as '[treatment]_[replicate]_[channel]_[fluorescence]_[voltage]_[magnification]_[run_no].ome.tif'
"""
exp_date = '2021-11-10'
path_dir = '/Volumes/LEMI_HK/LPS-DEP/' + exp_date + '/tif'
l = os.listdir(path_dir)

# remove unnecessary elements
ind = 0
while ind<len(l):
    if l[ind][0:2]=='._':
        del l[ind]
    else:
        ind = ind+1

# create and add video info
info = pd.DataFrame(columns=['date','channel','cond','rep','voltage','fps','front','back'])
for ind in range(len(l)):
    chl = l[ind]
    new_row = {'date':exp_date, 'channel':}
# %%
path_info = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/info_2021-10-23.txt'
path_plot = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis'
info = pd.read_csv(path_info, delimiter=',', header=0)
# np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=x_cols)

i = 3
path_tif = '/Volumes/LEMI_HK/LPS-DEP/XXXX-XX-XX/adjusted'
path_tif = path_tif.replace('XXXX-XX-XX',info.values[i,0])
s = path_tif + '/' + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.ome_v2.tif' % (info.cond[i], info.rep[i], info.channel[i], info.voltage[i])
frame = pims.open(s)

t1 = time.time()
f = sub.pile(frame[1:], diam=25, topn=10) # exclude the first frame it has been subtracted to remove background

pred = tp.predict.NearestVelocityPredict()
tr = pd.concat(pred.link_df_iter(f, search_range=25))
# tr = tp.link(f, 50); # not anymore
# tr = tr[(tr['x'] < x_hi) & (tr['x'] > x_lo)]; # not recommended to use

tr = sub.filter_ephemeral(tr, thres=5)

tr_v = sub.get_v(tr)
# sub.plot_tr_v(tr_v)

t2 = time.time();
print("elapsed : %s sec" % (t2-t1))

tr_v2 = sub.filter_v(tr_v, xlim=10, ylim1=20, ylim2=-10, direction=True)
sub.plot_tr_v(tr_v2)

info = pd.read_csv(path_info, delimiter=',', header=0) # update info
tr_v3 = sub.convert_tr(tr_v2, front=info.values[i,-2], back=info.values[i,-1])

tr_av = sub.each_particle(tr_v3, vol_init=15)

mu, tr_av2 = sub2.k_means(tr_av['mobility'])

# %% Export to comma delimited text file
path_sav = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis/'
tr_sav = pd.DataFrame(data = tr_av2, columns=['mobility'])
# tr_sav = get_tr_sav(tr_av, ind, info)   #ignore in this updated version 
s = path_sav + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.ome.csv' % (info.cond[i], info.rep[i], info.channel[i], info.voltage[i])
tr_sav.to_csv(s, index = False)

# %%
# sub.plot_v2(v2, path_plot, plotinfo='%s_R%d_%s' % (info.values[i,2], info.values[i,3], info.values[i,0]))
