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

import pims
import trackpy as tp
import time

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

os.chdir('/Volumes/LEMI_HK/LLNL BioSFA/EK/2021-04-25/tiff')
# os.chdir('./+1R1_R1-1');

s = '+1R1_R1_Ch05_TxRed_10-60V_1Vps_10X_001.ome';
i = 410;
frame = pims.open('%s.tif' % s);
# i = 0;
# frame = pims.open('+1R1_R1_Ch05_TxRed_10-60V_1Vps_10X_001.ome0401.jpg');

# frame_tmp = np.reshape(frame, (1,-1))
# h = plt.hist(frame_tmp[0], bins='auto')

diam = 31;

t1 = time.time();
# f = tp.locate(frame[i], diam, separation=1, percentile=80, minmass=100, threshold=15, invert=False);
f = tp.locate(frame[i], diam, invert=False);
# f = tp.batch(frame, diam, minmass=100);
t2 = time.time();
tp.annotate(f, frame[i])
print("elapsed : %s sec" % (t2-t1));