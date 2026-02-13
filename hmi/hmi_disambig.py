#!/usr/bin/env python
'''
Resolve the 180 ambiguity of the magnetic field
'''

from sunpython.fits_work import load_map, writefits


def hmi_disambig(azimuth, disambig, method=2, output='azimuth_disambig.fits',
                 save=False):
    '''
    Resolves the ambiguity in the transverse component of the magnetic field.

    Parameters:
        azimuth : astropy.io.fits.HDUList or sunpy.map.GenericMap
            Azimuth map with unresolved 180° ambiguity.
        disambig : astropy.io.fits.HDUList or sunpy.map.GenericMap
            Disambiguation map.
        method : int, optional
            From IDL, integer from 0 to 2 indicating the bit used. 0 for
            potential acute, 1 for random, 2 for radial acute (default).
            Out-of-range values are ignored and set to 2.
        output : str
            Name of output
        save : bool
            Write fits file

    Returns:
        sunpy.map.Map: Ambiguity-resolved azimuth map
    '''
    import numpy as np

    # Load data and scaled scale values
    azimuth_data, azimuth_meta, bitpix, bscale, bzero, blank = load_map(
        azimuth, get_scale=True)
    disambig_data, _ = load_map(disambig)

    azimuth_data = np.asarray(azimuth_data, dtype=float)
    disambig = np.asarray(disambig_data, dtype=np.int32)

    if azimuth_data.shape != disambig.shape:
        raise ValueError('Dimension of two images do not agree')

    if method < 0 or method > 2:
        method = 2
        print('Invalid disambiguation method, set to default method = 2')

    disambig_shifted = disambig // (2 ** method)
    mask = (disambig_shifted % 2) != 0

    azimuth_data[mask] += 180

    azimuth_meta['_DATAMAX'] = 360.0
    azimuth_meta.pop('DATAMEAN', None)
    azimuth_meta.pop('DATAMEDN', None)
    azimuth_meta.pop('DATARMS', None)
    azimuth_meta.pop('DATASKEW', None)
    azimuth_meta.pop('DATAKURT', None)

    map_azimuth_dis = writefits(azimuth_data, azimuth_meta, output,
                                bitpix, bscale, bzero, blank, save=save)

    return map_azimuth_dis


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            prog='hmi_disambig.py',
            description=(
                'Resolve the 180° ambiguity of the HMI magnetic field.'
            )
        )

        parser.add_argument(
            'azimuth',
            type=str,
            help='Azimuth FITS file with unresolved 180° ambiguity'
        )

        parser.add_argument(
            'disambig',
            type=str,
            help='Disambiguation FITS file'
        )

        parser.add_argument(
            'output',
            type=str,
            help='Output FITS file'
        )

        parser.add_argument(
            '-m', '--method',
            type=int,
            choices=[0, 1, 2],
            default=1,
            help=(
                'Bit used for disambiguation: (default: %(default)s)\n'
                '  0 = potential acute.\n'
                '  1 = random.\n'
                '  2 = radial acute.'
            )
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    from astropy.io import fits

    azimuth_map = fits.open(args.azimuth)
    disambig_map = fits.open(args.disambig)

    hmi_disambig(
        azimuth_map,
        disambig_map,
        method=args.method,
        output=args.output,
        save=True
    )
