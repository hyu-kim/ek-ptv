# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last modified on Jan 26 2022

For LLNL EK project. Updates --
Dec 08: Duplicated from 'main_bact.py'. Updated path info.
Jan 17: Info now includes pH and light filter
Jan 26: Test run for GSP2022 data

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
import math

# %%
exp_date = '2022-01-15'
path_info = '/Users/hyungseokkim/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' ep_ph/' + 'info_' + exp_date + '.txt'
path_plot = '/Users/hyungseokkim/Desktop/LEMI/SFA/Electrokinetics/' + exp_date + ' ep_ph/'
path_sav_vy = path_plot + 'vy/'
path_sav_tr = path_plot + 'tr/'
info = pd.read_csv(path_info, delimiter=',', header=0)

path_tif = '/Volumes/LEMI_HK/LLNL BioSFA/EK/XXXX-XX-XX/adjusted'
for i in range(34,41):
    path_tif = path_tif.replace('XXXX-XX-XX',info.date[i])
    s = path_tif + '/' + '%s_R%d_Ch%02d_pH%dp%02d_%s_%02dV_10X_001.ome.tif' % (info.cond[i], info.rep[i], info.channel[i], math.floor(info.ph[i]), round(100*(info.ph[i]-math.floor(info.ph[i]))), info.light[i], info.voltage[i])
    frame = pims.open(s)

    t1 = time.time()
    mid = (info.front[i] + info.back[i])//2
    b, cnt = sub.binarize(frame[mid])
    cnt = cnt * ((cnt>=10) & (cnt<=200)) + 10 * ((cnt<10) | (cnt>200))
    frame2, _ = sub.binarize_batch(frame) # for validating tp.annotate

    f = pile(frame[1:], topn=cnt*2//3) # exclude the first frame it has been subtracted to remove background
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
    tr_v3 = sub.convert_tr(tr_v2, front=info.front[i], back=info.back[i], rate_time=1/info.fps[i], rate_space=1.29) # mag 10x, binned 2 by 2

    tr_av = sub.each_particle(tr_v3, vol_init=info.voltage[i])

    mu, tr_av_vel = sub2.k_means(tr_av['velocity'])
    # tr_av_vel = tr_av['velocity']

    # %% Export tr_av and tr_av_vel to comma delimited text file
    s = path_sav_vy + '%s_R%d_Ch%02d_pH%dp%02d_%s_%02dV_10X_001.ome.txt' % (info.cond[i], info.rep[i], info.channel[i], math.floor(info.ph[i]), round(100*(info.ph[i]-math.floor(info.ph[i]))), info.light[i], info.voltage[i])
    s2 = path_sav_tr + '%s_R%d_Ch%02d_pH%dp%02d_%s_%02dV_10X_001.ome.txt' % (info.cond[i], info.rep[i], info.channel[i], math.floor(info.ph[i]), round(100*(info.ph[i]-math.floor(info.ph[i]))), info.light[i], info.voltage[i])

    tr_sav = pd.DataFrame(data = tr_av_vel, columns=['velocity'])
    tr_sav.to_csv(s, index = False)
    tr_av.to_csv(s2, index = False)

# %% FOR CORRECTING OUTPUT FILES IN TR / VY FOLDERS
# multiply by 1.95 times to correct pixel-micron ratio
l = os.listdir(path_plot + 'tr')
ind = 0
while ind<len(l):
    if l[ind][0]=='.':
        del l[ind]
    else:
        ind = ind+1
l.sort()
for s in l:
    mtx = pd.read_csv(path_plot+'tr/'+s, delimiter=',', header=0)
    mtx.velocity = mtx.velocity * 1.95
    mtx.mobility = mtx.mobility * 1.95
    mtx.zeta = mtx.zeta * 1.95
    mtx.to_csv(path_plot+'tr/'+s, index=False)

# %% re-save tr_av_vel withouth k-means filtering
for i in range(len(info)):
    s2 = path_sav_tr + '%s_R%d_Ch%02d_pH%dp%02d_%s_%02dV_10X_001.ome.txt' % (info.cond[i], info.rep[i], info.channel[i], math.floor(info.ph[i]), 100*(info.ph[i]-math.floor(info.ph[i])), info.light[i], info.voltage[i])
    tr_av = pd.read_csv(s2, delimiter=',', header=0)