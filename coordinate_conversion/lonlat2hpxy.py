#!/usr/bin/env python
'''
    Find x, y pixel location from Carrington Longitude and Latitude
    (Heliographic coordinates)

    Usage: lonlat2hpxy.py REF_FILE[FITS] CRLN CRLT

    Returns: (hpx, hpy)
'''


def lonlat2hpxy(map_input, crln, crlt):
    '''
    Converts from Heliographic Carrington to helioprojective cartesian

    Parameters
    ----------
        map_input : sunpy.map.GenericMap
        crln : float or astropy.units.degree
            Heliographic Carrington Longitude
        crlt : float or astropy.units.degree
            Heliographic Carrington Latitude

    Returns
    -------
        (float, float)
            Helioprojective cartesian (x, y) coordinates
    '''

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import frames

    # Calculate sun radius in meters
    rsun_obs = map_input.rsun_obs.to(u.rad)
    dsun_obs = map_input.dsun
    rsun = rsun_obs.value*dsun_obs

    # Calculate coordinate
    c = SkyCoord(
        crln, crlt, unit=u.deg,
        rsun=rsun,
        observer=map_input.observer_coordinate,
        frame=frames.HeliographicCarrington,
        obstime=map_input.reference_date,
    )

    hpx = c.helioprojective.Tx
    hpy = c.helioprojective.Ty

    return hpx, hpy


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            description=('Convert Carrington longitude/latitude to '
                         'helioprojective x/y.')
        )

        parser.add_argument(
            'input_fits',
            help='Input FITS file'
        )

        parser.add_argument(
            'crln',
            type=float,
            help='Carrington longitude (deg)'
        )

        parser.add_argument(
            'crlt',
            type=float,
            help='Carrington latitude (deg)'
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    # Check if file exists
    import os
    import sys
    if not os.path.exists(args.input_fits):
        sys.exit(f'File "{args.input_fits}" does not exist.')

    # ------------------------------------------------------------------

    from sunpy.map import Map

    # Load map
    map_input = Map(args.input_fits)

    # Convert lon/lat -> helioprojective x/y
    hpx, hpy = lonlat2hpxy(map_input, args.crln, args.crlt)

    print(hpx.value, hpy.value)
