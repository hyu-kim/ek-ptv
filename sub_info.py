# %%
"""
Last modified on Jan 17 2022
Jan 17: filename format now includes pH information
"""
import os
import pandas as pd
# %%
def str2df(s, date):
    """
    Extracts info from a filename .tif of a string, delmited by _

    Parameters
    ----------
    s : string. filename
    date : string. date
    
    Returns
    -------
    df : dataframe
         infos extracted from filename
    """
    i = 0
    j1 = 0 # read begins
    j2 = 0 # read ends
    l = [] # list of reads
    while s[i]!='.':
        if s[i]=='_':
            l = l + [s[j1:j2+1]]
            j1 = i+1
        else:
            j2 = i
        i = i+1
    df = pd.DataFrame({
        'date': [date], 'channel': [int(l[2][2:])], 'cond': [l[0]], 
        'rep': [int(l[1][1])], 'ph': [int(l[3][2]) + 0.01*float(l[3][4:6])], 
        'light': [l[4]], 'voltage': [l[5][:-1]],
        'fps':[10], 'front':[0], 'back':[0]
    })
    return df

def create_info(exp_date='2021-11-10', path='/Volumes/LEMI_HK/LPS-DEP/', path_out='/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/'):
    """
    Reads a list of files in a directory and creates a draft 'info.txt'
    format as '[treatment]_[replicate]_[channel]_[fluorescence]_[voltage]_[magnification]_[run_no].ome.tif'

    Parameters
    ----------
    exp_date : string
               date of experiment to analyze
    """
    path_dir = path + exp_date + '/tif'
    l = os.listdir(path_dir)

    # remove unnecessary elements
    ind = 0
    while ind<len(l):
        if l[ind][0:2]=='._':
            del l[ind]
        else:
            ind = ind+1
    l.sort()

    # create and add video info
    info = pd.DataFrame(columns=['date','channel','cond','rep','ph','light','voltage','fps','front','back'])
    for s in l:
        new_row = str2df(s, exp_date)
        info = info.append(new_row)

    # export to txt file
    info_dir = path_out + 'info_' + exp_date + '.txt'
    info.to_csv(info_dir, index = False)

def write_fps(exp_date='2021-11-10', path='/Volumes/LEMI_HK/LPS-DEP/'):
    """
    Reads a list of time record in subfolder 'time' then writes average fps into 'info_.txt'
    formatted '[treatment]_[replicate]_[channel]_[fluorescence]_[voltage]_[magnification]_[run_no].ome.tif'

    Parameters
    ----------
    exp_date : string
               date of experiment to analyze
    Returns
    -------
    """
    path_info = '/Users/hk/Desktop/LEMI/DEP-LPS/Linear EK/info_' + exp_date + '.txt'
    path_time = path + 'XXXX-XX-XX/time'
    path_time = path_time.replace('XXXX-XX-XX',exp_date)
    info = pd.read_csv(path_info, delimiter=',', header=0)
    for i in range(len(info)):
        s = path_time + '/' + '%s_R%d_Ch%02d_GFP_%02dV_20X_001.txt' % (info.cond[i], info.rep[i], info.channel[i], info.voltage[i])
        info_nd2 = pd.read_csv(s, header=59, sep='\t')
