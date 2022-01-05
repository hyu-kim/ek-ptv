# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last run on Jan 05 2021

For LPS-DEP project. Updates --
Oct 30: 1) file path, 2) noise reduction for bacterial cell tracking
Nov 30: Added a cell that reads a list of files in directory to create 'info.txt'
Dec 04: Use v_y instead of mobility for export
Dec 11: Export trace dataframe for record. Increase the upper cutoff filtering velocity range
Dec 17: 1) updated module sub ("get_v" and "convert_tr")
        2) Hold back using kmeans clustering -- as the mobility computed too high than expected

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
exp_date = '2022-01-03_BL21'
path_info = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/info_' + exp_date + '.txt'
path_plot = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis'
info = pd.read_csv(path_info, delimiter=',', header=0, dtype={'od': 'str'})

path_tif = '/Volumes/LEMI_HK/LPS-DEP/XXXX-XX-XX/adjusted'
for i in range(len(info)):
    path_tif = path_tif.replace('XXXX-XX-XX',info.date[i])
    s = path_tif + '/' + '%s_OD%s_Ch%02d_GFP_%02dV_10X_001.ome.tif' % (info.cond[i], info.od[i], info.channel[i], info.voltage[i])
    # s = path_tif + '/' + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.ome_test.tif' % (info.cond[i], info.od[i], info.channel[i], info.voltage[i])
    frame = pims.open(s)

    t1 = time.time()
    mid = (info.front[i] + info.back[i])//2
    b, cnt = sub.binarize(frame[mid])
    cnt = cnt * ((cnt>=10) & (cnt<=200)) + 10 * ((cnt<10) | (cnt>200))
    frame2, _ = sub.binarize_batch(frame) # for validating tp.annotate

    f = pile(frame[1:], diam=11, topn=cnt//2) # exclude the first frame it has been subtracted to remove background
    # f = pile(frame[1:], topn=5) # use this when cell number is too low (i.e. wrong cnt)
    # tp.annotate(f[100], frame2[100]) # run this to check if cells were properly detected

    pred = tp.predict.NearestVelocityPredict()
    tr = pd.concat(pred.link_df_iter(f, search_range=25))
    
    tr = sub.filter_ephemeral(tr, thres=5)

    tr_v = sub.get_v(tr)
    # sub.plot_tr_v(tr_v)

    t2 = time.time()
    print("elapsed : %s sec" % (t2-t1))

    tr_v2 = sub.filter_v(tr_v, xlim=5, ylim1=2*info.voltage[i], ylim2=-3, direction=True)
    sub.plot_tr_v(tr_v2)

    tr_v3 = sub.convert_tr(tr_v2, front=info.front[i], back=info.back[i], rate_time=1/info.fps[i], rate_space=1.95)

    tr_av = sub.each_particle(tr_v3, vol_init=info.voltage[i])

#     mu, tr_av_vel = sub2.k_means(tr_av['velocity']) # hold back until QW issue resolves

    # %% Export tr_av and tr_av_vel to comma delimited text file
    path_sav_vy = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis/' + exp_date + '/vy/'
    path_sav_tr = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/analysis/' + exp_date + '/tr/'
    # tr_sav = pd.DataFrame(data = tr_av_vel, columns=['velocity']) # hold back until QW issue resolves
    tr_sav = pd.DataFrame(data = tr_av['velocity'].values, columns=['velocity'])
    s = path_sav_vy + '%s_OD%s_Ch%02d_GFP_%02dV_10X_001.ome.csv' % (info.cond[i], info.od[i], info.channel[i], info.voltage[i])
    s2 = path_sav_tr + '%s_OD%s_Ch%02d_GFP_%02dV_10X_001.ome.csv' % (info.cond[i], info.od[i], info.channel[i], info.voltage[i])
    tr_sav.to_csv(s, index = False)
    tr_av.to_csv(s2, index = False)