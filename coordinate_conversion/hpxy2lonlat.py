#!/usr/bin/env python
'''
    Usage: xy2lonlat.py REF_FILE[FITS] hpx hpy

    Returns: (CRLN, CRLT)
'''


def hpxy2lonlat(map_input, hpx, hpy):
    '''
    Converts a coordinate pair from helioprojective to Heliographic Carrington

    Parameters
    ----------
        map_input : sunpy.map.GenericMap
        hpx : float or astropy.units.arcsec
            helioprojective cartesian x coordinate
        hpx : float or astropy.units.arcsec
            helioprojective cartesian y coordinate

    Returns
    -------
        (float, float)
            crln, crlt Heliographic Carrington location
    '''
    from sunpy.coordinates import frames
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    # Calculate sun radius in meters
    rsun_obs = map_input.rsun_obs.to(u.rad)
    dsun_obs = map_input.dsun
    rsun = rsun_obs.value*dsun_obs

    # _check_keywords(map_input)
    c = SkyCoord(
        hpx, hpy, unit=u.arcsec,
        rsun=rsun,
        observer=map_input.observer_coordinate,
        frame=frames.Helioprojective,
        obstime=map_input.reference_date,
    )

    lon = u.Quantity(c.heliographic_carrington.lon.value, u.deg)
    lat = u.Quantity(c.heliographic_carrington.lat.value, u.deg)
    return lon, lat


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            description=('Convert helioprojective x/y to Carrington longitude'
                         '/latitude.')
        )

        parser.add_argument(
            'input_fits',
            help='Input FITS file'
        )

        parser.add_argument(
            'hpx',
            type=float,
            help='Helioprojective X coordinate (arcsec)'
        )

        parser.add_argument(
            'hpy',
            type=float,
            help='Helioprojective Y coordinate (arcsec)'
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

    # Convert HP x/y to Carrington lon/lat
    crln, crlt = hpxy2lonlat(map_input, args.hpx, args.hpy)

    print(crln.value, crlt.value)
