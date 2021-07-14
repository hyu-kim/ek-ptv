#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 23:35:06 2021

@author: Hyu Kim (hskimm@mit.edu)

"""

import numpy as np
import matplotlib.pyplot as plt
import pims
import pandas as pd
import copy
import trackpy as tp
@pims.pipeline

def trans_contrast(frame, q1=0.3, q2=0.995):
    """
    Transforms a sequence of image to have a higher constrast
    
    Parameters
    ----------
    frame : PIMS frame object
        image sequence.
    q1 : float, between 0 and 1
        first quantile for lower bound. The default is 0.3.
    q2 : float, between 0 and 1
        second quantile for upper bound. The default is 0.995.

    Returns
    -------
    frame2 : image sequende with a higher contrast

    """
    frame2 = copy.deepcopy(frame);
    m = round(len(frame)/2);
    x = np.reshape(frame[m],(1,-1))[0];
    a,b,c = plt.hist(x,bins=255,density=True,range=(0,255),histtype='step',cumulative=True);
    thre1 = np.where(a<q1)[0][-1];
    sz2 = len(np.where(a>=q2)[0]);
    thre2 = np.where(a>=q2)[0][round(sz2*0.9)];
    
    for i in range(len(frame)):
        frame2[i][frame[i]<thre1] = thre1;
        frame2[i][frame[i]>thre2] = thre2;
        temp = (frame2[i] - thre1) / (thre2 - thre1) * 255;
        temp = np.around(temp);
        frame2[i] = temp.astype(np.uint16);
    return frame2

def pile(frame, diam, topn):
    """
    same as the function batch. Designed to work faster
    
    Parameters
    ----------
    frames : PIMS frame object. Sequence of images.
    diam : Integer, odd. See function locate
    topn : Integer. See function locate

    Returns
    -------
    DataFrame([x, y, mass, size, ecc, signal])

    """
    f = tp.locate(frame[0], diam, invert=False, topn=topn);
    f_tup = (f,);
    for i in range(len(frame)-1):
        f_tmp = tp.locate(frame[i+1], diam, invert=False, topn=topn);
        # f = pd.concat([f, f_tmp]);
        f_tup = f_tup + (f_tmp,);
    return f

def trshow(tr, first_style='bo', last_style='gs', style='b.'):
    frames = list(tr.groupby('frame'))
    nframes = len(frames)
    for i, (fnum, pts) in enumerate(frames):
        if i == 0:
            sty = first_style
        elif i == nframes - 1:
            sty = last_style
        else:
            sty = style
        plt.plot(pts.x, pts.y, sty, markersize=2)
    tp.plot_traj(tr, colorby='frame')
    # plt.axis('equal'); 
    # ylim(ymin=-1.0, ymax=3.5)
    plt.xlabel('x')
    plt.ylabel('y')

def est_vel(i, rate=0.15):
    """
    Estimates vetical velocity based on frame number.
    Assumes 0.14 s/frame

    Parameters
    ----------
    i : integer
        frame number.
    rate : float
           Rate of velocity increase [px/frame].

    Returns
    -------
    v : float
        estimating velocity [px/frame].
    """
    v = (i>70)*rate*(i-71);
    return v

@tp.predict.predictor
def predict(i, particle):
    velocity = np.array((0, est_vel(i)))
    return particle.pos + velocity * (i - particle.t)