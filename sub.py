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
    frame : image sequende with a higher contrast

    """
    m = round(len(frame)/2);
    x = np.reshape(frame[m],(1,-1))[0];
    a,b,c = plt.hist(x,bins=255,density=True,range=(0,255),histtype='step',cumulative=True);
    thre1 = np.where(a<q1)[0][-1];
    sz2 = len(np.where(a>=q2)[0]);
    thre2 = np.where(a>=q2)[0][round(sz2*0.9)]
    
    for i in range(len(frame)):
        frame[i][frame[i]<thre1] = thre1;
        frame[i][frame[i]>thre2] = thre2;
        temp = (frame[i] - thre1) / (thre2 - thre1) * 255;
        temp = np.around(temp);
        frame[i] = temp.astype(np.uint16);
    return frame

def pile():
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
    return 0

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
        plot(pts.x, pts.y, sty)
    trackpy.plot_traj(tr, colorby='frame', ax=gca())
    axis('equal'); ylim(ymin=-1.0, ymax=3.5)
    xlabel('x')
    ylabel('y')    