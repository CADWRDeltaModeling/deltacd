import pyhecdss
import pandas as pd
from numpy import arange, pi, array, zeros

import os
import sys
import string
import math
import numpy
from os import listdir
from math import cos, sin, tan, atan, sqrt, pi, pow


def parse_parts(path_parts):
    pi = [path_parts.find('%s=' % x) for x in 'ABCDEF']
    pi.append(len(path_parts))
    parts = [str.strip(path_parts[pi[i]+2:pi[i+1]]) for i in range(len(pi)-1)]
    return dict(zip(list('ABCDEF'), parts))


def hecdt2time(dstr):
    '''Takes a string in HEC DSS format and returns a time stamp
    E.g. 01OCT1921 2400 is parsed as a date and then 2400 is parsed as hour min
    '''
    d1, d2 = dstr.split(' ')
    return pd.to_datetime(d1)+pd.to_timedelta(int(d2[0:2]), 'h')+pd.to_timedelta(int(d2[2:4]), 'm')


def save_from_dssts_format(fname):
    df = pd.read_csv(fname,header=None)
    dssfname = df.iloc[0,0].strip()
    ai = df[df.iloc[:, 0].str.match('A=')].index
    ai=ai.append(pd.Index([len(df)-2])) # add end of the last data frame
    with pyhecdss.DSSFile(dssfname,create_new=True) as dh:
        for i in range(len(ai)-1):
            dfi = pd.to_numeric(df.iloc[ai[i]+4:ai[i+1]-1, 0])#.reset_index().drop('index', axis=1)
            path_parts = df.iloc[ai[i], 0]
            cunits = df.iloc[ai[i]+1, 0].strip()
            parts = parse_parts(path_parts)
            sdate = hecdt2time(df.iloc[ai[i]+3,0])
            freq=pyhecdss.DSSFile.get_freq_from_epart(parts['E'])
            dfi.index = pd.period_range(start=sdate-freq,freq=freq,periods=len(dfi))
            path='/'+'/'.join(parts.values())+'/'
            dh.write_rts(path, dfi, cunits, 'PER-AVER')

if __name__ == "__main__":
    pyhecdss.set_message_level(0)
    f1 = open(sys.argv[1])
    save_from_dssts_format(f1)
