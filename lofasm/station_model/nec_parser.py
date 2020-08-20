#!/usr/bin/env python3.7
import os
import re
import numpy as np
import pandas as pd
import sys

FREQ = int(sys.argv[1])
FILE = "../agmp/alpha_{0}_{1}_{freq}.out"
OILE = "./alpha_tpmp_{freq}.csv"
PW = re.compile ("PLANE WAVE - THETA:\s*(?P<theta>\d+\.\d*)\sdeg, PHI:\s*(?P<phi>\d+\.\d*)\sdeg")
MP_PATTERN = "{0}\s*{1}\s*-?\d+\.\d*\s*-?\d+\.\d*\s*-?\d+\.\d*\s*-?\d+\.\d*\s*-?\d+\.\d*E[+-]?\d*\s*-?\d+\.\d*E[+-]?\d*\s*(?P<mg>-?\d+\.\d*E[+-]?\d*)\s*(?P<ph>-?\d+\.\d*)"

########################
def action (ff, pw, mp):
    with open(ff, 'r') as f:
        FREAD = f.read()
        tp = np.array (pw.findall (FREAD), dtype=np.float32).T
        ii = np.array (mp.findall (FREAD), dtype=np.float32).T
    return tp,ii

dc = dict()
for xy,seg,tag in zip(['x', 'y'], [106,213], [21,42]):
    for pt in ['p', 't']:
        ff = FILE.format (xy, pt, freq=FREQ)
        print (ff)
        mp = re.compile (MP_PATTERN.format(seg, tag))
        tp,ii = action (ff, PW, mp)
        lab = xy + pt
        dc['{0}_theta'.format(lab)] = tp[0]
        dc['{0}_phi'.format(lab)] = tp[1]
        dc['{0}_mag'.format(lab)] = ii[0]
        dc['{0}_phase'.format(lab)] = ii[1]
########################
df = pd.DataFrame (dc)
df.to_csv (OILE.format(freq = FREQ), index=False)
