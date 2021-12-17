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
from skimage import measure, filters
# from skimage.filters import threshold_otsu



def binarize(frame, sigma=2):
    """
    Otsu's thresholding after gaussian filter

    Parameters
    ----------
    frame : array / PIMS frame object
            An image to binarize (bw)
    
    Returns
    -------
    binary : array
             Binarized image
    n : integer
        Number of counts
    """
    frame2 = filters.gaussian(frame, sigma=sigma)
    thresh = filters.threshold_otsu(frame2)
    binary = frame2 > thresh
    labels = measure.label(binary)
    n = np.amax(labels)
    
    return binary, n

def binarize_batch(frames, sigma=2, front=0, back=200):
    """
    Binarizes multiple images for smoother cell tracking
    Not likely to be used (incompatible with function "tp.locate")
    """
    frames_b = []
    for i in range(len(frames)):
        b, n = binarize(frames[i], sigma=sigma)
        b = b*255
        frames_b = frames_b + [b]
    b, n = binarize(frames[(front+back)//2], sigma=sigma)

    return frames_b, n

@pims.pipeline
def pile(frame, diam=15, topn=None, minmass=None):
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
    f = tp.locate(frame[0], diam, invert=False, topn=topn, minmass=minmass)
    f_tup = (f,)
    for i in range(len(frame)-1):
        f_tmp = tp.locate(frame[i+1], diam, invert=False, topn=topn, minmass=minmass)
        # f = pd.concat([f, f_tmp]);
        f_tup = f_tup + (f_tmp,)
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
        number of frames to be set as a threshold.

    Returns
    -------
    tr2 : Dataframe
        Filtered trace dataframe
    """
    tr_particle = tr.loc[:,'particle'];
    tr2 = pd.concat([tr.loc[(tr.loc[:,'particle']==i) & (len(tr.loc[tr_particle==i])>thres)] \
                     for i in range(max(tr_particle))]);
    return tr2

def get_v(tr):
    """
    Generates a trace dataframe added with velocity

    Parameters
    ----------
    tr : dataframe
        trace.

    Returns
    -------
    tr_v : dataframe
        trace added with velocity

    """
    tr_v = tr.copy()
    tr_v['v_x'] = float('nan')
    tr_v['v_y'] = float('nan')
    tr_v.frame = tr_v.frame - min(tr_v.frame)
    fr = 0
    i_prev = min(tr['particle'])-1

    # convert tr['particle'] into a list containing non-repeating values (to avoid case issue in for loop below)
    p_ext = []
    [p_ext.append(x) for x in tr['particle'] if x not in p_ext]
    p_ext.sort()

    for i in p_ext:
        fr = min(tr_v['frame'][tr_v['particle']==i])
        for j in range(sum(tr_v['particle']==i)-1):
            ind1 = (tr_v['particle']==i)&(tr_v['frame']==fr)
            fr = fr + 1
            ind2 = (tr_v['particle']==i)&(tr_v['frame']==fr)
            tr_v.loc[ind2, 'v_x'] = tr_v['x'][ind2].values[0] - tr_v['x'][ind1].values[0]
            tr_v.loc[ind2, 'v_y'] = tr_v['y'][ind2].values[0] - tr_v['y'][ind1].values[0]

    tr_v = tr_v[~np.isnan(tr_v['v_x'])]
    
    return tr_v

def plot_tr_v(tr_v):
    """
    Plots scatters of velocity from tr_v
    """
    mpl.rc('figure',  figsize=(10, 10));
    plt.figure();
    plt.plot(tr_v['v_x'], tr_v['v_y'], 'b.', markersize=2);
    plt.grid();

def filter_v(tr_v, xlim=5, ylim1=5, ylim2=-35, direction=True):
    """
    Excludes velocity values out of bound, unit in px/frame

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
        tr_v['v_y'] = -tr_v['v_y']
    tr_v2 = tr_v[(abs(tr_v['v_x'])<xlim) & (tr_v['v_y']<ylim1) & (tr_v['v_y']>ylim2)]
    # mpl.rc('figure',  figsize=(10, 10));
    # plt.figure();
    # plt.plot(tr_v2['v_x'], tr_v2['v_y'], 'b.', markersize=2);
    # plt.grid();
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
    # mpl.rc('figure',  figsize=(10, 10));
    # plt.figure();
    # plt.plot(tr_vf['v_x'], tr_vf['v_y'], 'b.', markersize=2);
    # plt.grid();
    
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
        
    # mpl.rc('figure',  figsize=(10, 10));
    # plt.figure();
    # plt.plot(range(n_fr), v, '-k', markersize=2);
    # plt.grid();
    return v

def convert_tr(tr_v, front, back, rate_time=0.1, rate_space=1): # default set to 10 fps and 20X with 3x3 binning
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
    """
    Exports figure (velocity plot) to a path
    """
    # mpl.rc('figure',  figsize=(10, 10));
    # plt.rc('font', size=16)
    # plt.figure();
    # plt.plot(v2[:,0], v2[:,1], '.k', markersize=5);
    # plt.grid();
    if path is not None:
        path = path + '/' + plotinfo
        plt.savefig(path)

def each_particle(tr_v, vol_init=10, ramp_rate=0):
    """
    Obtain average velocity, mobility etc for each particle

    Parameters
    ----------
    tr_v : Dataframe
           trace matrix with velocity in conversed units at each frame and particle
    ramp_rate : float
                increase rate of voltage by time (V/s). Default by zero
    Returns
    -------
    tr_avg : Dataframe
             calculation for each particle appearing in a sequence of frames
    """
    # constants
    eps_r = 80; # relative permitivity []
    eps_0 = 8.854e-12; # vacuum permitivity [F/m]
    eta = 8.66e-4; # viscosity [Pa-s]
    l = 10e-3; # channel length [m]
    n_prt = max(tr_v['particle'])
    colname = ['particle', 'time', 'duration', 'velocity', 'voltage', 'mobility', 'zeta']
    tr_avg = pd.DataFrame(columns = colname)
    
    for i in range(n_prt):
        ind = tr_v['particle']==i
        if len(tr_v['time'][ind]) > 4: # collects only when observed more than 4 frames
            ind2 = tr_v['time'][ind]!=min(tr_v['time'][ind])
            par = i # particle no
            t0 = min(tr_v.time[ind][ind2])
            t1 = max(tr_v['time'][ind][ind2])
            t_a = (t0 + t1) / 2 # time [sec]
            dur = t1 - t0 # duration [sec]
            vel = np.mean(tr_v['v_y'][ind][ind2]) # velocity [µm/s]
            vol = t_a * ramp_rate + vol_init # voltage, mutlplied by ramp rate [V]
            mob = vel / (vol/l) * 1e-6 # mobility, [m2/s-V]
            zet = mob * eta / (eps_r * eps_0) # zeta potential, [V]
            tup = pd.DataFrame([[par, t_a, dur, vel, vol, mob, zet]], columns = colname)
            tr_avg = pd.concat([tr_avg, tup], axis=0)
            # print('particle no.', i)

    return tr_avg

def get_tr_sav(tr_av, ind, info):
    """
    Exports mobility and zeta potential from tr_av with treatment information
    Update: no longer used since K-means algorithm

    Parameters
    ----------
    tr_av : Dataframe
            trace matrix with zeta and mobility
    ind : numeric
          index required to save
    info : dataframe
           matrix that includes experimetal condition
    """
    colname = ['mobility', 'zeta']
    tr_sav = pd.DataFrame(columns = colname)
    for i in tr_av['particle']:
        tup = pd.DataFrame([[
            tr_av.values[tr_av['particle']==i][0][5],
            tr_av.values[tr_av['particle']==i][0][6]
        ]], columns = colname)
        tr_sav = pd.concat([tr_sav, tup], axis=0)

    return tr_sav

# def trans_contrast(frame, q1=0.3, q2=0.995):
#     """
#     Transforms a sequence of image to have a higher constrast
    
#     Parameters
#     ----------
#     frame : PIMS frame object
#         image sequence.
#     q1 : float, between 0 and 1
#         first quantile for lower bound. The default is 0.3.
#     q2 : float, between 0 and 1
#         second quantile for upper bound. The default is 0.995.

#     Returns
#     -------
#     frame2 : image sequende with a higher contrast

#     """
#     frame2 = copy.deepcopy(frame);
#     m = round(len(frame)/2);
#     x = np.reshape(frame[m],(1,-1))[0];
#     a,b,c = plt.hist(x,bins=255,density=True,range=(0,255),histtype='step',cumulative=True);
#     thre1 = np.where(a<q1)[0][-1];
#     sz2 = len(np.where(a>=q2)[0]);
#     thre2 = np.where(a>=q2)[0][round(sz2*0.9)];
    
#     for i in range(len(frame)):
#         frame2[i][frame[i]<thre1] = thre1;
#         frame2[i][frame[i]>thre2] = thre2;
#         temp = (frame2[i] - thre1) / (thre2 - thre1) * 255;
#         temp = np.around(temp);
#         frame2[i] = temp.astype(np.uint16);
#     return frame2

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
