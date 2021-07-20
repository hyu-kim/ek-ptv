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
from pandas import DataFrame, Series  # for convenience
from sub import trans_contrast, trshow, pile, filter_ephemeral, filter_v, scatter_v
import pims
import trackpy as tp
import time

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(5, 10))
mpl.rc('image', cmap='gray')

# os.chdir('/Volumes/LEMI_HK/LLNL BioSFA/EK/2021-04-25/tiff')
os.chdir('/Users/hk/Desktop/Research/SFA/Electrokinetics/trackpy')

s = '+1R1_R1_Ch05_TxRed_10-60V_1Vps_10X_001.ome';
frame = pims.open('%s_v2.tif' % s); #v2, adjusted with imageJ
# frame2 = trans_contrast(frame, q1=0.7, q2=1-1e-3)

diam = 35;
i = 350;

t1 = time.time();
f0 = tp.locate(frame[i], diam, invert=False, topn=25);
f1 = tp.locate(frame[i+1], diam, invert=False, topn=25);
f = pile(frame, diam, topn=25);
# tp.quiet()
t2 = time.time();
# tp.annotate(f, frame[i])
print("elapsed : %s sec" % (t2-t1));

pred = tp.predict.NearestVelocityPredict();
# tr = pd.concat(pred.link_df_iter((f0, f1), search_range=40))
# tr = pd.concat(tp.link_df_iter((f0, f1), search_range=est_vel(i)))
tr = pd.concat(pred.link_df_iter(f, search_range=40));
# tr = tp.link(f, 50);
trshow(tr[5000:7000]);
# tr2 = pd.concat(pred.link_df_iter((f[100:125], f[125:150]), search_range=40));
# tr3 = pd.concat(pred.link_df_iter((f[100:125], f[125:150], f[150:175]), search_range=40));

# x_lo = 100;
# x_hi = 450;
# tr = tr[(tr['x'] < x_hi) & (tr['x'] > x_lo)]; # not recommended to use
tr = filter_ephemeral(tr);
tr_v = scatter_v(tr);
tr_v2 = filter_v(tr_v);
