import numpy as np
import lofasm_nec as lfn
import sys
FREQ = int(sys.argv[1])
#############
FREQS = np.array ([FREQ], dtype=np.uint8)
LF    =  lfn.LoFASM_NEC (freqs=FREQS, root_dir="tt/")
THETA = np.arange (91)
PHI   = np.arange (360)
XT    = []
YT    = []
XP    = []
YP    = []
#############
for t in THETA:
    for p in PHI:
        rr = LF.response (t, p)
        XT.append (rr['xt'])
        YT.append (rr['yt'])
        XP.append (rr['xp'])
        YP.append (rr['yp'])
#############
XT = np.array (XT)
YT = np.array (YT)
XP = np.array (XP)
YP = np.array (YP)
#############
ddd = {'xt':XT, 'yt':YT, 'xp':XP, 'yp':YP}
with open ("./lofasm_bp/lofasm_bp_{0:d}.pkl".format(FREQ), 'wb') as f:
    import pickle as pkl
    pkl.dump (ddd, f)
