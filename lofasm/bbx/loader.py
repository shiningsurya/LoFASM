"""
LoFASM bbx loaders
"""
import sys
import os

import gzip
import numpy as np
import struct
from copy import deepcopy
from random import choice
from string import ascii_letters, digits
from shutil import copyfileobj

import bbx

AUTOPOLS  = ['AA', 'BB', 'CC', 'DD']
CROSSPOLS = ['AB', 'AC', 'AD', 'BC', 'BD', 'CD']
ALLPOLS   = AUTOPOLS + CROSSPOLS


class DimSet (object):
    """
    Class object to handle dimension properties

    Start, stop, step, number
    """
    def __init__ (self, start, stop=None, step=None, n=None, span=None,label=None):
        """
        Arguments:
            start (float):
            stop  (float)
        """
        self.start   = start
        self.stop    = stop
        self.step    = step
        self.n       = n
        self.label   = label
        self.range   = span
        if self.stop is None and self.range is not None:
            self.stop = self.start + self.range
        #
        if self.step is None and self.n is None and self.stop is None:
            raise ValueError ("Need to pass atleast one of step, number, and stop.")
        elif self.step is not None and self.n is not None:
            self.range  = self.n * self.step
            self.stop   = self.start + self.range
        elif self.step is not None and self.stop is not None:
            self.range  = self.stop - self.start
            self.n      = int (self.range / self.step)
        elif self.stop is not None and self.n is not None:
            self.range  = self.stop - self.start
            self.step   = self.range / self.n

    def get_array (self, use_linspace=False):
        """Generates an array"""
        if use_linspace:
            return np.linspace (self.start, self.stop, self.n)
        else:
            return np.arange (self.start, self.stop, self.step)

    def __repr__ (self):
        return "{0} start={1} stop={2} n={3:d} step={4}".format (self.label,self.start, self.stop, self.n, self.step)

    def find_closest (self, x, return_float=False):
        """
        Finds the closest index to x
        """
        dx = ( x - self.start ) / self.step
        if return_float:
            return x
        else:
            return int (dx)

    def __eq__ (self, o):
        a   = self.start == o.start
        b   = self.stop  == o.stop
        c   = self.step  == o.step
        return a and b and c

class LofasmDataloader (object):
    """
    Taken from calibrate/lofasmcal.py
    Edited accordingly

    Since this functionality is very generic, 
    moving it up 
    """
    def __init__(self, files, verbose=False):
        """Reads lofasm files and outputs data ready for calibration.
            Arguments:
                files (str): wildcard string
        """
        self.filelist = sorted (glob.glob (files))
        self.nfiles   = len (self.filelist)
        self.v        = verbose

        #Remove non lofasm files from filelist
        for f in reversed(self.filelist):
            if not bbx.is_lofasm_bbx(f):
                self.filelist.remove(f)
        # read/populate times and stuff
        self.times = []
        self.freqs = []
        self.corr  = []
        self.sta   = []
        self.mjd   = []
        self.__read_header__ ()
        self.index  = np.argsort (self.mjd)

    def __read_header__ (self,):
        for i, f in enumerate (self.filelist):
            bx     = bbx.LofasmFile (f)
            head   = bx.header
            hmet   = head['metadata']
            #
            t_start  = head['start_time']
            t_bins   = hmet['dim1_len']
            t_len    = float (head['dim1_span'])
            t_dt     = self.__resolve_time (t_start)
            t_dim    = DimSet (t_dt, n=t_bins, span=t_len)
            #
            f_bins   = int (hmet['dim2_len'])
            f_bw     = float (head['dim2_span']/1e6)
            f_start  = float (head['dim2_start']/1e6)
            f_dim    = DimSet (f_start, n=f_bins, span=f_bw, label='Freq')
            # 
            self.times.append (t_dim)
            self.freqs.append (f_dim)
            self.corr.append  (head['channel'])
            self.sta.append   (head['station'])
            self.mjd.append   (head['mjd'])

    def __resolve_time (self, tt):
        """
        Convert header start_time to datetime object
        ================TODO=====================================
        Due to a bug some LoFASM data has messed up start times.
        What follows is a hotfix so that this library can support
        both time formats. In the future, the bug should be fixed
        and this library should only have to support 1 time format.
        
        The following loop tries each format specified in fmts.
        If the time is parsed correctly (and does not throw an
        exception) then the loop will be broken.
        ==========================================================
        """
        fmts = ["%Y-%m-%dT%H:%M:%S",  # correct time format for bbx files
                "%Y%m%d_%H%M%S"]  # additional format found in LoFASM I files
        startt_repr = tt[:-8] if 'T' in tt else tt 
        time_obj = None
        for fmt in fmts:
            try:
                time_obj = datetime.datetime.strptime(startt_repr, fmt)
                break
            except ValueError:
                pass
        else:
            print ( "Cannot parse start time header field {}".format(tt) )
        return time_obj

    def read_freq (self, freq):
        """
        freq in MHz
        """
        total_n  = sum([ds.n for ds in self.times])
        ret = np.empty (total_n, dtype=np.float32)
        II,JJ = 0,0
        for i,f in enumerate (self.filelist):
            ifreq = self.freqs[i].find_closest (freq)
            JJ    = II + self.times[i].n
            bf    = bbx.LofasmFile (f)
            bf.read_data ()
            ret[II:JJ] = bf.data[ifreq, :]
            bf.close ()
            II    = JJ
        return ret
    
    def read_ifreq (self, ifreq):
        """
        ifreq is channel index
        can be array or list
        """
        if isinstance (ifreq, int):
            ifreq = [ifreq]
        nf       = len(ifreq)
        total_n  = sum([ds.n for ds in self.times])
        ret = np.empty ((nf, total_n), dtype=np.float32)
        II,JJ = 0,0
        for i,f in enumerate (self.filelist):
            JJ    = II + self.times[i].n
            bf    = bbx.LofasmFile (f)
            bf.read_data ()
            for iidx, ichan in enumerate(ifreq):
                ret[iidx, II:JJ] = bf.data[ichan, :]
            bf.close ()
            II    = JJ
        return ret

    def read_files(self, files, freq, verbose=True):
        """Creates data attribute using the minimum value in frequency and time
           for each file.

           PRELIMINARY
           DONOT USE
        """
        self.freqmhz = freq

        re = 'Reading data...'
        for i in range(len(self.filelist)):

            # Find the frequency bin corresponding to the given frequency
            bw = (float(head['dim2_span'])/1000000.0)/head['metadata']['dim2_len']
            freqbin = int((self.freqmhz-(float(head['dim2_start'])/1000000.0))/bw)
            self.freq_bin = np.append(self.freq_bin, freqbin)

            f.read_data()
            freqbin_range = 10 #bins around specified 'freq' to read
            freqbins = f.data[freqbin-freqbin_range/2:freqbin+freqbin_range/2,:]
            min_vals_per_timebin = [np.min(x) for x in np.rot90(freqbins, k=-1)]
            min_file_val = np.min(min_vals_per_timebin)
            f.close()

            # Compute datetime of the selected timebin
            index = min_vals_per_timebin.index(min_file_val)
            seconds_into_file = (float(index)/timebins)*timelength
            seconds_into_file = datetime.timedelta(seconds=seconds_into_file)
            time_from_bin = time_obj + seconds_into_file

            self.data = np.append(self.data, min_file_val)
            self.times_array = np.append(self.times_array, time_from_bin)

            if verbose == True:
                p = (str(i*100/len(self.filelist)) + '%')
                if i+1 not in range(len(self.filelist)):
                    p = 'Done'
                sys.stdout.write("\r%s%s" % (re,p))
                sys.stdout.flush()

if __name__ == "__main__":
    print ()
    a = DimSet (100., span=20., step=1.0, label='Test')
    print (a)
    import datetime as dt
    start = dt.datetime.now()
    dt    = dt.timedelta (minutes=2)
    b = DimSet (start, n=10, step=dt)
    print (b)
