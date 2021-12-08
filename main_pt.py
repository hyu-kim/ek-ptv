# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Dec 08 2021

For LLNL EK project. Updates --
Dec 08: Duplicated from 'main_bact.py'. Updated path info.

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

# %%
exp_date = '2021-11-29'
path_info = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' Pt mobility 2/' + 'info_' + exp_date + '.txt'
path_plot = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' Pt mobility 2/analysis'
info = pd.read_csv(path_info, delimiter=',', header=0)

path_tif = '/Volumes/LEMI_HK/LLNL-BioSFA/EK/XXXX-XX-XX/adjusted'
for i in range(len(info)):
    path_tif = path_tif.replace('XXXX-XX-XX',info.date[i])
    s = path_tif + '/' + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.ome.tif' % (info.cond[i], info.rep[i], info.channel[i], info.voltage[i])
    frame = pims.open(s)

    t1 = time.time()
    mid = (info.front[i] + info.back[i])//2
    b, cnt = sub.binarize(frame[mid])
    frame2, _ = sub.binarize_batch(frame) # for validating tp.annotate

    f = pile(frame[1:], topn=cnt//2) # exclude the first frame it has been subtracted to remove background
    # f = pile(frame[1:], topn=5) # use this when cell number is too low (i.e. wrong cnt)
    # tp.annotate(f[100], frame2[100]) # run this to check if cells were properly detected

    pred = tp.predict.NearestVelocityPredict()
    tr = pd.concat(pred.link_df_iter(f, search_range=25))
    # tr = tp.link(f, 50); # not anymore
    # tr = tr[(tr['x'] < x_hi) & (tr['x'] > x_lo)]; # not recommended to use

    tr = sub.filter_ephemeral(tr, thres=5)

    tr_v = sub.get_v(tr)
    # sub.plot_tr_v(tr_v)

    t2 = time.time()
    print("elapsed : %s sec" % (t2-t1))

    tr_v2 = sub.filter_v(tr_v, xlim=10, ylim1=info.voltage[i], ylim2=-10, direction=True)
    sub.plot_tr_v(tr_v2)

    info = pd.read_csv(path_info, delimiter=',', header=0) # update info
    tr_v3 = sub.convert_tr(tr_v2, front=info.front[i], back=info.back[i], rate_time=1/info.fps[i])

    tr_av = sub.each_particle(tr_v3, vol_init=info.voltage[i])

    mu, tr_av_vel = sub2.k_means(tr_av['velocity'])

    # %% Export to comma delimited text file
    path_sav = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis/' + exp_date + '/'
    tr_sav = pd.DataFrame(data = tr_av_vel, columns=['velocity'])
    # tr_sav = get_tr_sav(tr_av, ind, info)   #ignore in this updated version 
    s = path_sav + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.ome.csv' % (info.cond[i], info.rep[i], info.channel[i], info.voltage[i])
    tr_sav.to_csv(s, index = False)