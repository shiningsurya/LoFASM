import os
import numpy as np
import pandas as pd
import glob

SPEED_LIGHT = 299792458.0 #Units = m/s


class LoFASM_NEC (object):
    """LoFASM NEC modeling"""
    def __init__ (self, 
            Nring=2, 
            Nperring=6, 
            freqs=np.arange(10, 89, dtype=np.uint8), 
            root_dir="./", 
            initial_r=441, 
            r_fact=np.sqrt(3),
            file_format = "alpha_tpmp_{0:d}.csv",
            ntheta = 91,
            nphi   = 360,
        ):
        """
        Arguments:
            freqs (array, optional):  NumPy array of LoFASM 
            root_dir (str, optional): Directory holding the csv files
        """
        self.ir       = initial_r
        self.rfact    = r_fact
        self.nrings   = Nring
        self.nperring = Nperring
        ## theta, phi
        self.ntheta   = ntheta
        self.nphi     = nphi
        ## freqs, wave vectors
        self.freqs    = freqs
        omegas        = 2 * np.pi * self.freqs * 1e6
        self.k        = omegas / SPEED_LIGHT
        ###
        # geometry
        # planar geometry so only two components
        # first dim - x, second dim - y
        self.geom = np.zeros ((self.nrings, self.nperring, 2))
        p_per_ring = np.linspace (0, 2*np.pi, self.nperring)
        for i in range (self.nrings):
            this_r  = self.ir * self.rfact ** i
            self.geom[i, :, 0] = this_r * np.cos (p_per_ring)
            self.geom[i, :, 1] = this_r * np.sin (p_per_ring)
        ###
        # read single antenna response into memory
        self.single = dict()
        for f in self.freqs:
            fff = os.path.join (root_dir, file_format.format(f))
            self.single[f] = pd.read_csv (fff)

    def resolve_index (self, theta, phi):
        """
        Extremely dependent on the csv
        faster in theta, 
        """
        return int( phi*self.ntheta + theta )

    def response (self, theta, phi):
        """
        Response w.r.t to the station center.
        """
        idx = self.resolve_index (theta, phi)
        hat = np.array ([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta)])
        ##
        xt  = np.zeros (self.freqs.shape, dtype=np.complex)
        yt  = np.zeros (self.freqs.shape, dtype=np.complex)
        xp  = np.zeros (self.freqs.shape, dtype=np.complex)
        yp  = np.zeros (self.freqs.shape, dtype=np.complex)
        # work loop
        for i,f in enumerate (self.freqs):
            ttt  = self.single[f].iloc[idx]
            for ir in range (self.nrings):
                for ip in range (self.nperring):
                    ppp  = -1j * self.k[i] * np.dot (self.geom[ir,ip], hat)
                    xt[i] += ttt.xt_mag * np.exp ( (1j*ttt.xt_phase) + ppp )
                    yt[i] += ttt.yt_mag * np.exp ( (1j*ttt.yt_phase) + ppp )
                    xp[i] += ttt.xp_mag * np.exp ( (1j*ttt.xp_phase) + ppp )
                    yp[i] += ttt.yp_mag * np.exp ( (1j*ttt.yp_phase) + ppp )
        # compute powers and return
        ret = dict()
        ret['xt'] = np.abs (xt)
        ret['yt'] = np.abs (yt)
        ret['xp'] = np.abs (xp)
        ret['yp'] = np.abs (yp)
        return ret

if __name__ == "__main__":
    test_freqs = np.array ([74], dtype=np.uint8)
    test_lf = LoFASM_NEC (freqs=test_freqs,root_dir="tt/")
    test_theta = 45
    test_phi   = 60
    test_response = test_lf.response (test_theta, test_phi)
    print ("For theta={0} phi={1}".format (test_theta, test_phi))
    print ("Response is {0}".format (test_response))


