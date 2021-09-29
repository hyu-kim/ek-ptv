# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on 9/28/2021
@author: Hyu Kim (hskimm@mit.edu)

Objectives
1) Obtain every cell speed and create barplot (no statitic use)
2) Analyze every cell under increasing voltage
"""
from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sub, pims
import trackpy as tp
import time
# %% Track single cells and obtain zeta
path_info = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/2020-09-25 Pt mobility/code/info.txt'
path_plot = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/2020-09-25 Pt mobility/plot'
info = pd.read_csv(path_info, delimiter=',', header=0)

ind = 1;
path_tif = '/Volumes/LEMI_HK/LLNL BioSFA/EK/XXXX-XX-XX/tif_v2'
path_tif = path_tif.replace('XXXX-XX-XX',info.values[ind,0])
s = path_tif + '/' + '%s_R%d_Ch%02d_TxRed_10-60V_1Vps_10X_001.ome_v2.tif' % (info.values[ind,2], info.values[ind,3], info.values[ind,1])
frame = pims.open(s)
# t1 = time.time();
f = sub.pile(frame, diam=35, topn=25);
pred = tp.predict.NearestVelocityPredict();
tr = pd.concat(pred.link_df_iter(f, search_range=40));
tr = sub.filter_ephemeral(tr);
tr_v = sub.scatter_v(tr);
# print("elapsed : %s sec" % (time.time()-t1));
tr_v2 = sub.filter_v(tr_v, xlim=4, ylim1=10, ylim2=-40, direction=True);

# v = sub.plot_v_quantile(tr_v2, 'mean');
# info = pd.read_csv(path_info, delimiter=',', header=0) # updated info after seeing v

tr_v3 = convert_tr(tr_v2, front=info.values[ind,-2], back=info.values[ind,-1])
tr_av = each_particle(tr_v3)

##%% Export to comma delimited text file
path_sav = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/2020-09-25 Pt mobility/single_cell/'
tr_sav = get_tr_sav(tr_av, ind, info)
s = path_sav + 'Ch%02d_%s_R%d_single.csv' % (info['channel'][ind], info['cond'][ind], info['rep'][ind])
tr_sav.to_csv(s, index = False)