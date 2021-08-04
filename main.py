#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 02:46:22 2021

@author: Hyu Kim (hskimm@mit.edu)
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

os.chdir('/Users/hk/Desktop/Research/SFA/Electrokinetics/2020-09-25 Pt mobility/code')
info = pd.read_csv('info.txt', delimiter=',', header=0)
# np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=x_cols)

path_tif = '/Volumes/LEMI_HK/LLNL BioSFA/EK/XXXX-XX-XX/tif_v2'
exp_date = ['2020-09-17','2021-04-25']
i = 0;
path_tif = path_tif.replace('XXXX-XX-XX',info.values[i,0])
os.chdir(path_tif)

# s = '+1R1_R1_Ch05_TxRed_10-60V_1Vps_10X_001';
# frame = pims.open('%s.ome_v2.tif' % s); #v2, adjusted with imageJ
frame = pims.open('%s_R%d_Ch%02d_TxRed_10-60V_1Vps_10X_001.ome_v2.tif' % (info.values[i,2], info.values[i,3], info.values[i,1]));


diam = 35;
t1 = time.time();
f0 = tp.locate(frame[i], diam, invert=False, topn=25);
f1 = tp.locate(frame[i+1], diam, invert=False, topn=25);
f = sub.pile(frame, diam, topn=25);
# tp.quiet()
t2 = time.time();
# tp.annotate(f, frame[i])
print("elapsed : %s sec" % (t2-t1));

pred = tp.predict.NearestVelocityPredict();
# tr = pd.concat(pred.link_df_iter((f0, f1), search_range=40))
# tr = pd.concat(tp.link_df_iter((f0, f1), search_range=est_vel(i)))
tr = pd.concat(pred.link_df_iter(f, search_range=40));
# tr = tp.link(f, 50);
sub.trshow(tr[5000:7000]);
# tr2 = pd.concat(pred.link_df_iter((f[100:125], f[125:150]), search_range=40));
# tr3 = pd.concat(pred.link_df_iter((f[100:125], f[125:150], f[150:175]), search_range=40));

# x_lo = 100;
# x_hi = 450;
# tr = tr[(tr['x'] < x_hi) & (tr['x'] > x_lo)]; # not recommended to use
tr = sub.filter_ephemeral(tr);
tr_v = sub.scatter_v(tr);
tr_v = sub.filter_v(tr_v, xlim=3);
v = sub.plot_v_quantile(tr_v, 'q1');
v2 = sub.convert_vy(v, front=info.values[i,-2], back=info.values[i,-1]);
