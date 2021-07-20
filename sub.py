#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 23:35:06 2021

@author: Hyu Kim (hskimm@mit.edu)

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
    tupule of frames

    """
    f = tp.locate(frame[0], diam, invert=False, topn=topn);
    f_tup = (f,);
    for i in range(len(frame)-1):
        f_tmp = tp.locate(frame[i+1], diam, invert=False, topn=topn);
        # f = pd.concat([f, f_tmp]);
        f_tup = f_tup + (f_tmp,);
    return f_tup

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
    tp.plot_traj(tr, colorby='frame');
    plt.grid();
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

def filter_ephemeral(tr, thres=9):
    """
    Removes traces appearing at no more than specific number of frames

    Parameters
    ----------
    tr : Dataframe
        trace matrix size of (topn * no. frames, 10).
    thres : Integer
        maximum number of frames to filter.

    Returns
    -------
    tr2 : Dataframe
        Filtered trace dataframe
    """
    tr_particle = tr.loc[:,'particle'];
    tr2 = pd.concat([tr.loc[(tr.loc[:,'particle']==i) & (len(tr.loc[tr_particle==i])>thres)] \
                     for i in range(max(tr_particle))]);
    return tr2

def scatter_v(tr):
    """
    Plots scatter plot of trace dataframe. Generates a trace dataframe added with velocity

    Parameters
    ----------
    tr : dataframe
        trace.

    Returns
    -------
    tr_v : dataframe
        trace added with velocity

    """
    tr_v = tr.copy();
    tr_v['v_x'] = tr['ep'];
    tr_v['v_y'] = tr['ep'];
    fr = 0;
    i_prev = min(tr['particle'])-1;
    for i in tr['particle']:
        if i==i_prev:
            ind1 = (tr['particle']==i)&(tr['frame']==fr);
            fr = fr+1;
            ind2 = (tr['particle']==i)&(tr['frame']==fr);
            tr_v.loc[ind2, 'v_x'] = tr_v['x'][ind2].values[0] - tr_v['x'][ind1].values[0];
            tr_v.loc[ind2, 'v_y'] = tr_v['y'][ind2].values[0] - tr_v['y'][ind1].values[0];
        else:
            fr = min(tr['frame'][tr['particle']==i]);
        i_prev = i;
    mpl.rc('figure',  figsize=(10, 10));
    plt.plot(tr_v['v_x'], tr_v['v_y'], 'b.', markersize=2);
    plt.grid();
    return tr_v

def filter_v(tr_v, xlim=5, ylim1=5, ylim2=-35, direction=True):
    """
    Excludes velocity values out of bound

    Parameters
    ----------
    tr_v : Dataframe
        input trace+velocity matrix.
    xlim : Float
        upper and lower bounds of v_x
    ylim1 : Float
        upper bound of v_y
    ylim2 : Float
        lower bound of v_y
    direction : Bool
        True if velocity_y < 0, False else

    Returns
    -------
    tr_v2 : dataframe, trace+velocity filtered

    """
    if not direction:
        tr_v['v_y'] = -tr_v['v_y'];
    tr_v2 = tr_v[(abs(tr_v['v_x'])<xlim) & (tr_v['v_y']<ylim1) & (tr_v['v_y']>ylim2)];
    mpl.rc('figure',  figsize=(10, 10));
    plt.plot(tr_v2['v_x'], tr_v2['v_y'], 'b.', markersize=2);
    plt.grid();
    return tr_v2

def distrib_v(tr_v, fr):
    """
    Plots distribution of velocity x,y of a frame (sequence)

    Parameters
    ----------
    tr_v : Dataframe
        trace with velocity matrix.
    fr : Numeric or list
        Frame.

    Returns
    -------
    None.

    """
    tr_vf = tr_v[tr_v['frame']==fr];
    mpl.rc('figure',  figsize=(10, 10));
    plt.plot(tr_vf['v_x'], tr_vf['v_y'], 'b.', markersize=2);
    plt.grid();
    
def plot_v_quantile(tr_v, s):
    """
    Plots a statistic of v_y (max / quantile) by frame number.

    Parameters
    ----------
    tr_v : Dataframe
        trace matrix.
    s : Str
        'max', 'mean', 'q1', 'q2'

    Returns
    -------
    None.

    """
    if s=='max':
        fn = max;
    elif s=='q1':
        fn = lambda x: np.quantile(x, 0.5);
        
    for i in range(max(tr_v['frame'])):
        tr_vf = tr_v['v_y'][tr_v['frame']==i];
        print("%s" % fn(tr_vf));
        