import numpy as np
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

if __name__ == "__main__":
    hs  = Galaxy_Model ("./lambda_haslam408_dsds.fits")
