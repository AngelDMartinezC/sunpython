#!/usr/bin/env python
'''
Code to convert from a Line-Of-Sight projected map into Postel.
The output is a 1024
'''


def correct_solar_rotation(map_input):

    import numpy as np
    from sunpy.map import Map

    map_input = Map(map_input)
    data = map_input.data  # in m/s
    ny, nx = data.shape
    nan_pos = np.isnan(data)
    data[nan_pos] = 0

    y, x = np.mgrid[:ny, :nx]  # pixel coordinates
    x = x - nx/2  # center coordinates
    y = y - ny/2

    # Flatten arrays for fitting
    X = np.vstack([x.ravel(), y.ravel(), np.ones(x.size)]).T
    v = data.ravel()

    # Solve for plane coefficients
    coeffs, _, _, _ = np.linalg.lstsq(X, v, rcond=None)
    a, b, c = coeffs
    plane = a*x + b*y + c
    doppler_derot = data - plane

    # If NaNs where present, put them back
    doppler_derot[nan_pos] = np.nan
    doppler_clean = Map(doppler_derot, map_input.meta)
    return doppler_clean


def los2postel(map_input, crln, crlt, naxis=(1024, 1024),
               daxis=(0.0005, 0.0005), name_output=None, doppler=False,
               continuum=False, algorithm='interpolation',
               xy_ref=(None, None), save=False):
    '''
    Read info of a Line Of Sight map and reprojects it to Postel.

    Parameters
    ----------
    map_input : sunpy.map.GenericMap
        LOS projected map
    crln : float
        Carrington longitude of center of projection
    crlt : float
        Carrington latitude of center of projection
    naxis : (int, int)
        (NAXIS1, NAXIS2) of the projected frame
    daxis : (float, float)
        (CDELT1, CDELT2) in solar radii. Default to (0.0005, 0.0005)
    doppler : bool
        Apply solar rotation correction.
    continuum : bool
        Calculate limb darkening correction.
    algorithm : {'interpolation', 'adaptive', 'exact'}
        Projection algorithm (from fastest to slowest).
    name_output : str
        Name of output file.
    xy_ref : (float, float)
        (x, y) pixel coordinate pair of reference of projection. Default is
        (NAXIS1/2 + 0.5, NAXIS2/2 + 0.5).

    Returns
    -------
    sunpy.map.Map
        Reprojected map
    '''

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from sunpy.map import make_fitswcs_header
    from sunpy.coordinates import frames
    from numpy import pi
    from sunpython.fits_work import writefits, load_map
    from sunpython.hmi import limbcorrect
    import warnings
    from sunpy.util.exceptions import SunpyUserWarning
    from astropy.io.fits.verify import VerifyWarning

    # Filter some annoying warnings
    warnings.filterwarnings(
        'ignore',
        message='rsun mismatch detected',
        category=SunpyUserWarning
    )
    warnings.filterwarnings(
        'ignore',
        category=VerifyWarning,
        message='Keyword name.*greater than 8 characters'
    )

    # Get information on the bscale, bzero, blank keywords before runtime
    _, _, bitpix, bscale, bzero, blank = load_map(map_input, get_scale=True)

    # Calculate sun radius in meters
    rsun_obs = map_input.rsun_obs.to(u.rad)  # From arcsec to rad
    dsun = map_input.dsun
    rsun = rsun_obs.value*dsun

    # Temporal change of 'RSUN_REF'. It reverts to the other one later
    rsun_ref = map_input.meta['RSUN_REF']
    map_input.meta['RSUN_REF'] = rsun.value

    if 'RSUN_ARC' in map_input.meta:
        map_input.meta['RSUN_OBS'] = map_input.meta['RSUN_ARC']

    if continuum:
        map_input = limbcorrect(map_input, save=False)

    # Create center of projection
    origin_hpc = SkyCoord(
        crln, crlt, unit=u.deg,
        rsun=rsun,
        observer=map_input.observer_coordinate,
        frame=frames.HeliographicCarrington,
        obstime=map_input.reference_date,
    )

    origin = origin_hpc.heliographic_carrington
    out_shape = (naxis[0], naxis[1])
    scale1 = 180*daxis[0]/pi
    scale2 = 180*daxis[1]/pi

    # This ensures the data is well centered
    if xy_ref == (None, None):
        x0, y0 = (naxis[0]/2 + 0.5, naxis[1]/2 + 0.5)
    else:
        x0, y0 = xy_ref

    out_header = make_fitswcs_header(
        out_shape,
        origin,
        scale=[scale1, scale2]*u.deg/u.pix,
        reference_pixel=[x0, y0]*u.pix,
        projection_code='ARC'
    )

    map_out = map_input.reproject_to(
        out_header,
        algorithm=algorithm,
        parallel=False
    )

    # Scale data
    if 'BUNIT' in map_input.meta:
        if map_input.meta['BUNIT'] == 'DN':
            exptime = map_input.meta.get('EXPTIME',
                                         map_input.meta.get('XPOSURE', None))
            if exptime is not None:
                map_out.data[:, :] = map_out.data[:, :]/exptime
    if 'PIXLUNIT' in map_input.meta:
        if map_input.meta['PIXLUNIT'] == 'DN':
            exptime = map_input.meta.get('EXPTIME',
                                         map_input.meta.get('XPOSURE', None))
            if exptime is not None:
                map_out.data[:, :] = map_out.data[:, :]/exptime

    # Prefer CADENCE, else TRECSTEP
    if 'CADENCE' in map_input.meta:
        map_out.meta['CADENCE'] = map_input.meta['CADENCE']
    elif 'TRECSTEP' in map_input.meta:
        map_out.meta['CADENCE'] = map_input.meta['TRECSTEP']

    # Keyword condition added because of Solar Orbiter's deprecated DATE-OBS
    if 'DATE_EAR' in map_input.meta:
        map_out.meta['DATE-OBS'] = map_input.meta['DATE_EAR']
        map_out.meta['T_REC'] = map_input.meta['DATE_EAR']
    # Keyword condition for HMI
    else:
        map_out.meta['DATE-OBS'] = map_input.meta['DATE-OBS']
        if 'T_REC' in map_input.meta:
            map_out.meta['T_REC'] = map_input.meta['T_REC']
        else:
            map_out.meta['T_REC'] = map_input.meta['DATE-OBS']

    # Keyword added for SOLO/PHI
    if 'RSUN_OBS' in map_input.meta:
        map_out.meta['RSUN_OBS'] = map_input.meta['RSUN_OBS']
    else:
        map_out.meta['RSUN_OBS'] = map_input.meta['RSUN_ARC']

    # if 'WAVELNTH' in map_input.meta:
    #     map_out.meta['WAVELNTH'] = (
    #         map_input.meta['WAVELNTH'], '[angstrom] Wavelength')

    if 'WAVEUNIT' in map_input.meta:
        map_out.meta['WAVEUNIT'] = 'angstrom'

    if 'BUNIT' in map_input.meta:
        map_out.meta['BUNIT'] = map_input.meta['BUNIT']

    # Add new or update keywords
    keywords = {
        'WCSNAME': 'ARC',
        'DAXIS1': daxis[0],
        'DAXIS2': daxis[1],
        'RSUN_REF': rsun_ref,
        'CRPIX1': map_out.meta['CRPIX1'],
        'CRPIX2': map_out.meta['CRPIX2'],
        'CDELT1': map_out.meta['CDELT1'],
        'CDELT2': map_out.meta['CDELT2'],
        'CRVAL1': map_out.meta['CRVAL1'],
        'CRVAL2': map_out.meta['CRVAL2'],
        'HGLN_OBS': map_out.meta['HGLN_OBS'],
        'HGLT_OBS': map_out.meta['HGLT_OBS'],
        'DSUN_OBS': map_input.meta['DSUN_OBS'],
        'CRLN_OBS': map_input.meta['CRLN_OBS'],
        'CRLT_OBS': map_input.meta['CRLT_OBS'],
        'OBS_L0': map_input.meta['CRLN_OBS'],
        'OBS_B0': map_input.meta['CRLT_OBS'],
        'REF_L0': crln,
        'REF_B0': crlt,
    }

    map_out.meta['keycomments'] = {
        'WCSNAME':  'Coordinate system',
        'DAXIS1':   'solar radii per pixel',
        'DAXIS2':   'solar radii per pixel',
        'RSUN_REF': '[m] Reference radius of the Sun',
        'CRPIX1':   '[pixel]',
        'CRPIX2':   '[pixel]',
        'CDELT1':   '[deg/pixel] scale',
        'CDELT2':   '[deg/pixel] scale',
        'CRVAL1':   '[deg] origin at CRPIX1',
        'CRVAL2':   '[deg] origin at CRPIX2',
        'HGLN_OBS': '[deg] Stonyhurst longitude',
        'HGLT_OBS': '[deg] Stonyhurst latitude',
        'DSUN_OBS': '[m] Distance from Sun',
        'CRLN_OBS': '[deg] Carrington longitude of observer',
        'CRLT_OBS': '[deg] Carrington latitude of observer',
        'OBS_L0': '[deg] Same as {CRLN_OBS}',
        'OBS_B0': '[deg] Same as {CRLT_OBS}',
    }

    for key, value in keywords.items():
        map_out.meta[key] = value

    # Manually add units to have proper map.spatial_units
    map_out.meta['CUNIT1'] = 'deg'
    map_out.meta['CUNIT2'] = 'deg'
    map_out.meta['CUNIT3'] = 's'
    map_out.meta['CTYPE1'] = 'CRLN-ARC'
    map_out.meta['CTYPE2'] = 'CRLT-ARC'
    map_out.meta['CTYPE3'] = 'TIME'

    map_out.plot_settings = map_input.plot_settings

    if doppler:
        map_out = correct_solar_rotation(map_out)

    if save:
        map_out = writefits(map_out.data, map_out.meta, name_output,
                            bitpix, bscale, bzero, blank, save=save)

    return map_out


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            prog='los2postel.py',
            description=(
                'Read a Line-Of-Sight (LOS) map and reproject it to '
                'a Postel projection.'
            )
        )

        # Positional arguments
        parser.add_argument(
            'input',
            help='Input LOS FITS file'
        )

        parser.add_argument(
            'crln', type=float,
            help='Carrington longitude of the projection center'
        )

        parser.add_argument(
            'crlt', type=float,
            help='Carrington latitude of the projection center'
        )

        parser.add_argument(
            'output',
            help='Output FITS file'
        )

        # Optional arguments
        parser.add_argument(
            '--algorithm', choices=['interpolation', 'adaptive', 'exact'],
            default='interpolation',
            help=(
                'Projection algorithm (from fastest to slowest). '
                'Default: interpolation')
        )

        group = parser.add_mutually_exclusive_group()

        group.add_argument(
            '--int',
            action='store_true',
            help='Enable limb darkening correction (default: False)'
        )

        group.add_argument(
            '--dopp',
            action='store_true',
            help='Apply solar rotation (Doppler) correction'
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    from sunpy.map import Map

    # Load map
    map_input = Map(args.input)

    # Run projection
    map_out = los2postel(
        map_input,
        crln=args.crln,
        crlt=args.crlt,
        algorithm=args.algorithm,
        doppler=args.dopp,
        name_output=args.output,
        continuum=args.int,
        save=True
    )
