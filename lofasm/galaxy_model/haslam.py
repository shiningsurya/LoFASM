import numpy as np
import pkg_resources
import healpy as hp

class Galaxy_Model (object):
    """Interfaces with the galaxy modeling"""
    def __init__ (self, 
            mapfile, 
            map_freq=408., 
            power_scale=-2.55
        ):
        self.mfile          = mapfile
        self.hmap           = hp.read_map (mapfile, verbose=False)
        self.map_freq       = map_freq
        self.power_scale    = power_scale

    def __call__ (self, gl, gb, freq=74.):
        """interpolates
        """
        freq_fac = ( freq / self.map_freq ) ** self.power_scale
        return hp.get_interp_val (self.hmap, gl, gb, lonlat=True) * freq_fac

class Haslam408 (Galaxy_Model):
    """Inherits Galaxy_Model"""
    def __init__ (self, power_scale=-2.55, map_freq=408.,):
        mapfile = pkg_resources.resource_filename (__name__, 'lambda_haslam408_dsds.fits')
        super (Haslam408, self).__init__ (mapfile, map_freq=map_freq, power_scale=power_scale)

if __name__ == "__main__":
    hs  = Galaxy_Model ("./lambda_haslam408_dsds.fits")
