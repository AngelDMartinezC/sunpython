import logging
import sys
import warnings
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
logging.getLogger("sunpy").setLevel(logging.WARNING)
# Suppress all SunPy warnings
warnings.filterwarnings("ignore", module="sunpy")
from sunpy.coordinates import frames
from sunpy.coordinates.utils import get_heliocentric_angle
from sunpy.map import Map

smap_str = sys.argv[1]
hpx = float(sys.argv[2])
hpy = float(sys.argv[3])

smap = Map(smap_str)

# At the center of the solar disk
hpc_coord_center = SkyCoord(
        hpx*u.arcsec,
        hpy*u.arcsec,
        # frame='helioprojective',
        # observer="earth",
        # obstime="2024-12-08"
        frame=frames.Helioprojective,
        observer=smap.observer_coordinate,
        obstime=smap.date,
        )

angle = get_heliocentric_angle(hpc_coord_center)

# mu
mu = np.cos(angle.to_value(u.rad))

print(mu)
