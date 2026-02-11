#!/usr/bin/env python
'''
    Find Carrington Longitude and Latitude (Heliographic coordinates) from
    a point on SDO images with Helioprojective Cartesian Coordinate System.
    Converts from pix to LonLat (Heliographic Carrington)

    Usage: xy2lonlat.py REF_FILE[FITS] x0 y0

    Returns: (CRLN, CRLT)
'''

import astropy.units as u


def xy2lonlat(map_input, x0, y0):
    '''
    Converts a coordinate pair from pixel to Heliographic Carrington

    Parameters
    ----------
        map_input : sunpy.map.GenericMap
        x0 : float or astropy.units.pixel
            zero-based index x-coordinate in pixels
        y0 : float or astropy.units.pixel
            zero-based index y-coordinate in pixels

    Returns
    -------
        (float, float)
            crln, crlt Heliographic Carrington location
    '''

    x0 = u.Quantity(x0, u.pix)
    y0 = u.Quantity(y0, u.pix)

    # Convert to heliographic
    coord = map_input.pixel_to_world(x0, y0)
    coord_car = coord.heliographic_carrington

    lon = u.Quantity(coord_car.lon.value, u.deg)
    lat = u.Quantity(coord_car.lat.value, u.deg)

    return lon, lat


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            description=(
                'Convert pixel coordinates (x, y) to Carrington lon/lat.'
            )
        )

        parser.add_argument(
            'input_fits',
            help='Input FITS file'
        )

        parser.add_argument(
            'x0',
            type=float,
            help='X pixel coordinate (zero-based index)'
        )

        parser.add_argument(
            'y0',
            type=float,
            help='Y pixel coordinate (zero-based index)'
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

    # Compute lon/lat
    crln, crlt = xy2lonlat(map_input, args.x0, args.y0)

    print(crln.value, crlt.value)
