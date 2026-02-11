#!/usr/bin/env python
'''
    Find x, y pixel location from Heliprojective cartesian

    Usage: hpxy2xy.py REF_FILE[FITS] hpx hpy

    Returns: (x, y)
'''


def hpxy2xy(map_input, hpx, hpy):
    '''
    Converts a coordinate pair from helioprojective to pixel coordinates

    Parameters
    ----------
        map_input : sunpy.map.GenericMap
        hpx : float or astropy.units.arcsec
            helioprojective cartesian x coordinate
        hpx : float or astropy.units.arcsec
            helioprojective cartesian y coordinate

    Returns
    -------
        (astropy.units.pix, astropy.units.pix)
            zero-based index (x, y) coordinates
    '''

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.coordinates import frames

    # Calculate sun radius in meters
    rsun_obs = map_input.rsun_obs.to(u.rad)
    dsun_obs = map_input.dsun
    rsun = rsun_obs.value*dsun_obs

    c = SkyCoord(
        hpx, hpy, unit=u.arcsec,
        rsun=rsun,
        observer=map_input.observer_coordinate,
        frame=frames.Helioprojective,
        obstime=map_input.reference_date,
    )

    coord = map_input.world_to_pixel(c)

    xpix = coord.x
    ypix = coord.y

    return xpix, ypix


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            description='Convert helioprojective x/y to pixel coordinates.'
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
    data = Map(args.input_fits)

    # Convert HP x/y to pixel x/y
    x, y = hpxy2xy(data, args.hpx, args.hpy)

    print(x.value, y.value)
