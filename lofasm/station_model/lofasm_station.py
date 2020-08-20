import os
import numpy as np
import pickle as pkl
import astropy.time as at
import astropy.coordinates as asc
import astropy.units as au
import scipy.interpolate as sin
"""
EarthLocations taken from shiningsurya@github.com/pulsar_search
"""
LF1_LAT   = asc.Latitude ( "26:33:19.676", unit=au.degree)
LF1_LONG  = asc.Longitude( "-97:26:31.174",unit=au.degree)
LF1_ELVF  = ( 10.36 * au.imperial.foot ).to (au.meter)
LF2_LAT   = asc.Latitude ( "34:4:43.497" , unit=au.degree)
LF2_LONG  = asc.Longitude( "-107:37:5.819",unit=au.degree)
LF2_ELVF  = ( 6967.08 * au.imperial.foot ).to (au.meter)
LF3_LAT   = asc.Latitude ( "38:25:59.0" ,  unit=au.degree)
LF3_LONG  = asc.Longitude( "-79:50:23.0",  unit=au.degree)
LF3_ELVF  = ( 2464.89 * au.imperial.foot ).to (au.meter)
LF4_LAT   = asc.Latitude ( "34:12:3.0" ,   unit=au.degree)
LF4_LONG  = asc.Longitude( "-118:10:18.0", unit=au.degree)
LF4_ELVF  = ( 1167.89 * au.imperial.foot ).to (au.meter)
# Earth locations
LF1       = asc.EarthLocation (LF1_LONG, LF1_LAT, LF1_ELVF)
LF2       = asc.EarthLocation (LF2_LONG, LF2_LAT, LF2_ELVF)
LF3       = asc.EarthLocation (LF3_LONG, LF3_LAT, LF3_ELVF)
LF4       = asc.EarthLocation (LF4_LONG, LF4_LAT, LF4_ELVF)

class LoFASM_Station (object):
    """Models LoFASM station coordinate systems and stuff"""
    def __init__ (self, 
            name,
            el,
            nalt=910,
            naz=3590,
        ):
        #
        self.name       = name
        self.location   = el
        ##
        self.gal        = asc.Galactic
        self.alt        = np.linspace (0., 90.0, nalt) * au.degree
        self.az         = np.linspace (0., 360.0, naz) * au.degree
        self.tt, self.pp= np.meshgrid (self.alt, self.az)

    def __quadrature_sum (self, x, y):
        return np.sqrt ( np.power(x, 2) + np.power(y,2) )

    def beam_pattern (self, freq, root_dir="./"):
        """
        Reads the pickle file and interpolates the beam pattern
        """
        FF = "lofasm_bp_{0:d}.pkl"
        with open (os.path.join (root_dir, FF.format(freq)), 'rb') as f:
            # bp = pkl.load (f, encoding='ascii')
            bp = pkl.load (f,encoding='latin1')
        xt = bp['xt'].reshape ((91,-1))
        yt = bp['yt'].reshape ((91,-1))
        yp = bp['yp'].reshape ((91,-1))
        xp = bp['xp'].reshape ((91,-1))
        AXYTP = 0.5 * ( self.__quadrature_sum (xt, yt) + self.__quadrature_sum (xp, yp) )
        thetas  = np.arange (91.0)
        phis    = np.arange (360.0)
        ####
        interp  = sin.RectBivariateSpline (thetas, phis, AXYTP)
        ZZ = interp (self.alt, self.az)
        return ZZ

    def __call__ (self, t):
        """Selects the t
            t should be utc and in datetime or isot format
        Arguments:
            t (datetime, str) : should be parsable by astropy.Time
        """
        ttime = at.Time (t, scale='utc', location=self.location)
        aa = asc.AltAz (alt=self.tt, az=self.pp, obstime=ttime, location=self.location)
        gc = aa.transform_to (self.gal)
        return gc
