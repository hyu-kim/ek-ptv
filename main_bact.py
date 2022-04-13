# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last run on Feb 18 2022

For LLNL-EK project (git branch = "amp"). Updates --
Oct 30: 1) file path, 2) noise reduction for bacterial cell tracking
Nov 30: Added a cell that reads a list of files in directory to create 'info.txt'
Dec 04: Use v_y instead of mobility for export
Dec 11: Export trace dataframe for record. Increase the upper cutoff filtering velocity range
Dec 17: 1) updated module sub ("get_v" and "convert_tr")
        2) Hold back using kmeans clustering -- as the mobility computed too high than expected
Feb 18: Bring code from "main_pt.py"

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
exp_date = '2022-04-06'
path_info = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' bact ek/' + 'info_' + exp_date + '.txt'
path_plot = '/Users/hk/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' bact ek/'
path_sav_vy = path_plot + 'vy/'
path_sav_tr = path_plot + 'tr/'
info = pd.read_csv(path_info, delimiter=',', header=0)

path_tif = '/Volumes/LEMI_HK/LLNL BioSFA/EK/XXXX-XX-XX/tif_v2'
for i in range(30,60):
    print('iteration', i)
    
    path_tif = path_tif.replace('XXXX-XX-XX',info.date[i])
    s = path_tif + '/' + '%s_R%d_Ch%02d_%s_%02dV_10X_001.ome.tif' % (info.cond[i], info.rep[i], info.channel[i], info.light[i], info.voltage[i])
    frame = pims.open(s)

    t1 = time.time()
    mid = (info.front[i] + info.back[i])//2
    b, cnt = sub.binarize(frame[mid])
    cnt = cnt * ((cnt>=10) & (cnt<=100)) + 10 * (cnt<10) + 100 * (cnt>100)
    frame2, _ = sub.binarize_batch(frame) # for validating tp.annotate

    f = pile(frame[1:], topn=cnt) # exclude the first frame it has been subtracted to remove background
    # f = pile(frame[1:], topn=5) # use this when cell number is too low (i.e. wrong cnt)
    # tp.annotate(f[100], frame2[100]) # run this to check if cells were properly detected

    pred = tp.predict.NearestVelocityPredict()
    tr = pd.concat(pred.link_df_iter(f, search_range=40))

    tr = sub.filter_ephemeral(tr, thres=5)

    tr_v = sub.get_v(tr)
    # sub.plot_tr_v(tr_v)

    t2 = time.time()
    print("elapsed : %s sec" % (t2-t1))

    tr_v2 = sub.filter_v(tr_v, xlim=10, ylim1=info.voltage[i], ylim2=-10, direction=True)
    sub.plot_tr_v(tr_v2)

    info = pd.read_csv(path_info, delimiter=',', header=0) # update info
    tr_v3 = sub.convert_tr(tr_v2, front=info.front[i], back=info.back[i], rate_time=1/info.fps[i], rate_space=1.29) # mag 10x, binned 2 by 2

    tr_av = sub.each_particle(tr_v3, vol_init=info.voltage[i])

    mu, tr_av_vel = sub2.k_means(tr_av['velocity'])
    # tr_av_vel = tr_av['velocity']

    # %% Export tr_av and tr_av_vel to comma delimited text file
    s = path_sav_vy + '%s_R%d_Ch%02d_%s_%02dV_10X_001.ome.csv' % (info.cond[i], info.rep[i], info.channel[i], info.light[i], info.voltage[i])
    s2 = path_sav_tr + '%s_R%d_Ch%02d_%s_%02dV_10X_001.ome.csv' % (info.cond[i], info.rep[i], info.channel[i], info.light[i], info.voltage[i])

    tr_sav = pd.DataFrame(data = tr_av_vel, columns=['velocity'])
    tr_sav.to_csv(s, index = False)
    tr_av.to_csv(s2, index = False)