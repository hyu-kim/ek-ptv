#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 23:35:06 2021

@author: Hyu Kim (hskimm@mit.edu)

"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pims
import pandas as pd
import copy
import trackpy as tp

#%%
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
    tr_v = tr_v[~np.isnan(tr_v['v_x'])];
    mpl.rc('figure',  figsize=(10, 10));
    plt.figure();
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
    plt.figure();
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
    plt.figure();
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
    v : Array
        list of statistics of v_y by each frame

    """
    if s=='max':
        fn = min;
    elif s=='q1':
        fn = lambda x: np.quantile(x, 0.25);
    elif s=='mean':
        fn = np.mean;
    elif s=='q2':
        fn = lambda x: np.quantile(x, 0.5);
    else:
        raise ValueError("try different statistic")
    
    n_fr = max(tr_v['frame']);
    v = np.zeros(n_fr);
    for i in range(n_fr-1):
        tr_vf = tr_v['v_y'][tr_v['frame']==i+1];
        v[i] = fn(tr_vf);
        
    mpl.rc('figure',  figsize=(10, 10));
    plt.figure();
    plt.plot(range(n_fr), v, '-k', markersize=2);
    plt.grid();
    return v

def convert_vy(v, front, back, rate_time=0.138, rate_space=1.288):
    """
    Trims starting / ending frames. Converts units in time / space

    Parameters
    ----------
    v : np array
        list of v_y statistics per frame [px/frame]
    front : numeric
        initial frame number to include.
    back : numeric
        last frame number to include.
    rate_time : float
        conversion factor [s/frame]
    rate_space : float
        conversion factor [µm/px]

    Returns
    -------
    v2 : array
        converted v_y.

    """
    v = v[front:back+1] * rate_space / rate_time;
    t = np.array(range(len(v))) * rate_time;
    v2 = np.hstack((np.transpose([t]), np.transpose([-v])));
    return v2

def convert_tr(tr_v, front, back, rate_time=0.138, rate_space=1.288):
    """
    Trims the first and last frames and convert units in time / space
    
    Parameters
    ----------
    tr_v : dataframe
           trace matrix with velocity
    front : numeric
            starting frame number
    back : numeric
           ending frame number
    rate_time : float
        conversion factor [s/frame]
    rate_space : float
        conversion factor [µm/px]
    """
    tr_v2 = tr_v.copy()
    tr_v2 = tr_v2[(tr_v2['frame']>=front)&(tr_v2['frame']<back)]
    tr_v2['frame'] = (tr_v2['frame'] - front) * rate_time
    tr_v2['x'] = tr_v2['x'] * rate_space
    tr_v2['y'] = tr_v2['y'] * rate_space
    tr_v2['v_x'] = tr_v2['v_x'] * rate_space / rate_time
    tr_v2['v_y'] = tr_v2['v_y'] * rate_space / rate_time
    tr_v2 = tr_v2.rename({'frame':'time'}, axis=1)
    return tr_v2

def calc_param(v2):
    """
    Calculates mobility and zeta potential from time-velocity data via linear regression
    
    """
    # constants
    eps_r = 80; # relative permitivity []
    eps_0 = 8.854e-12; # vacuum permitivity [F/m]
    eta = 8.66e-4; # viscosity [Pa-s]
    l = 10e-3; # channel length [m]
    
    coef = np.polyfit(v2[:,0], v2[:,1], 1) # [0,] µm / V-s
    mu = coef[0] * l * 1e-6;
    zeta = mu * eta / (eps_r * eps_0)
    return mu, zeta

def plot_v2(v2, path = None, plotinfo = None):
    mpl.rc('figure',  figsize=(10, 10));
    plt.rc('font', size=16)
    plt.figure();
    plt.plot(v2[:,0], v2[:,1], '.k', markersize=5);
    # plt.grid();
    if path is not None:
        path = path + '/' + plotinfo
        plt.savefig(path)

def each_particle(tr_v):
    """
    Obtain average velocity, mobility etc for each particle

    Parameters
    ----------
    tr_v : Dataframe
           trace matrix with velocity in conversed units at each frame and particle
    Returns
    -------
    tr_zeta : Dataframe
              mean velocity and zeta potential for each particle (averaged by sequence frames)
    """
    # constants
    eps_r = 80; # relative permitivity []
    eps_0 = 8.854e-12; # vacuum permitivity [F/m]
    eta = 8.66e-4; # viscosity [Pa-s]
    l = 10e-3; # channel length [m]
    n_prt = max(tr_v['particle'])
    colname = ['particle', 'time', 'duration', 'velocity', 'voltage', 'mobility', 'zeta']
    tr_avg = pd.DataFrame(columns = colname)
    
    for i in range(min(tr_v['particle']), max(tr_v['particle'])+1):
        ind = tr_v['particle']==i
        if len(tr_v['time'][ind]) > 0:
            par = i # particle no
            t0 = min(tr_v['time'][ind])
            t1 = max(tr_v['time'][ind])
            t_a = (t0 + t1) / 2 # time sec
            dur = t1 - t0 # duration sec
            vel = np.mean(tr_v3['v_y'][ind]) # velocity µm/s
            vol = t_a # voltage, ramp rate 1 Vps
            mob = vel / (vol/l) * 1e-6 # mobility, [m2/s-V]
            zet = mob * eta / (eps_r * eps_0) # zeta potential, V
            tup = pd.DataFrame([[par, t_a, dur, vel, vol, mob, zet]], columns = colname)
            tr_avg = pd.concat([tr_avg, tup], axis=0)

    return tr_avg