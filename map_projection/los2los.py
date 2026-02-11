#!/usr/bin/env python
'''
Code to convert from a Line-Of-Sight projected map into a cropped LOS.
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


def los2los(map_input, crln, crlt, naxis=(1024, 1024), name_output=None,
            doppler=False, continuum=False, save=False,
            algorithm='interpolation'):
    '''
    Read info of a Line Of Sight map and cuts it.

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
    doppler : bool
        Apply solar rotation correction.
    continuum : bool
        Calculate limb darkening correction.
    name_output : str
        Name of output file.
    algorithm : {'interpolation', 'adaptive', 'exact'}
        Projection algorithm (from fastest to slowest).

    Returns
    -------
    sunpy.map.Map
        Reprojected map
    '''

    import astropy.units as u
    from sunpython.fits_work import writefits, load_map
    from sunpython.hmi import limbcorrect
    import warnings
    from sunpy.util.exceptions import SunpyUserWarning
    from astropy.io.fits.verify import VerifyWarning
    from sunpython.coordinate_conversion import lonlat2hpxy
    from astropy.coordinates import SkyCoord
    from sunpy.map import make_fitswcs_header
    from sunpy.coordinates import frames
    from sunpy.map import Map

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

    hpx, hpy = lonlat2hpxy(map_input, crln, crlt)

    # Calculate sun radius in meters
    rsun_obs = map_input.rsun_obs.to(u.rad)  # From arcsec to rad
    dsun = map_input.dsun
    rsun = rsun_obs.value*dsun

    if 'RSUN_ARC' in map_input.meta:
        map_input.meta['RSUN_OBS'] = map_input.meta['RSUN_ARC']

    if continuum:
        map_input = limbcorrect(map_input, save=False)

    # # If CROTA angle is less than limit, then do not rotate
    # angle_lim = 0.05
    # if 'CROTA' in map_input.meta:
    #     angle_rot = map_input.meta['CROTA']
    #     map_input = map_input.rotate(angle_rot*u.deg)
    # elif 'CROTA2' in map_input.meta:
    #     angle_rot = map_input.meta['CROTA2']
    #     # TODO! Temporal solution to an AIA rotation issue
    #     if angle_rot > angle_lim:
    #         map_input = map_input.rotate(angle_rot*u.deg)

    # Temporal change of 'RSUN_REF'. It reverts to the other one later
    # rsun_ref = map_input.meta['RSUN_REF']
    # map_input.meta['RSUN_REF'] = rsun.value
    # Create center of projection
    origin_hpc = SkyCoord(
        hpx, hpy, unit=u.arcsec,
        rsun=rsun,
        observer=map_input.observer_coordinate,
        frame=frames.Helioprojective,
        obstime=map_input.reference_date,
    )

    origin = origin_hpc.helioprojective  # .heliographic_carrington
    out_shape = (naxis[0], naxis[1])
    scale1 = map_input.meta['CDELT1']
    scale2 = map_input.meta['CDELT2']

    xy_ref = (None, None)
    # This ensures the data is well centered
    if xy_ref == (None, None):
        x0, y0 = (naxis[0]/2 + 0.5, naxis[1]/2 + 0.5)
    else:
        x0, y0 = xy_ref

    out_header = make_fitswcs_header(
        out_shape,
        origin,
        scale=[scale1, scale2]*u.arcsec/u.pix,
        reference_pixel=[x0, y0]*u.pix,
        projection_code='TAN'
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
    #         map_input.meta['WAVELNTH'], '[Angstrom] Wavelength')

    if 'WAVEUNIT' in map_input.meta:
        map_out.meta['WAVEUNIT'] = 'angstrom'

    if 'BUNIT' in map_input.meta:
        map_out.meta['BUNIT'] = map_input.meta['BUNIT']

    keywords = {
        'CRVAL1': u.Quantity(map_out.meta['CRVAL1']*u.deg, u.arcsec).value,
        'CRVAL2': u.Quantity(map_out.meta['CRVAL2']*u.deg, u.arcsec).value,
        'CDELT1': u.Quantity(map_out.meta['CDELT1']*u.deg, u.arcsec).value,
        'CDELT2': u.Quantity(map_out.meta['CDELT2']*u.deg, u.arcsec).value,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'WCSNAME': 'Helioprojective-Cartesian',
        'CRLN_OBS': map_input.meta['CRLN_OBS'],
        'CRLT_OBS': map_input.meta['CRLT_OBS'],
        'OBS_L0': map_input.meta['CRLN_OBS'],
        'OBS_B0': map_input.meta['CRLT_OBS'],
        # 'CADENCE': trecstep,
        # 'TRECSTEP': trecstep,
    }

    map_out.meta['keycomments'] = {
        'CDELT1': '[arcsec/pixel] image scale in the x direction',
        'CDELT2': '[arcsec/pixel] image scale in the y direction',
        'CRVAL1': '[arcsec] CRVAL1: x origin',
        'CRVAL2': '[arcsec] CRVAL2: y origin',
        'CRPIX1': '[pixel]',
        'CRPIX2': '[pixel]',
        'CRLN_OBS': '[deg] Carrington longitude of the observer',
        'CRLT_OBS': '[deg] Carrington latitude of the observer',
        'OBS_L0': '[deg] Same as {CRLN_OBS}',
        'OBS_B0': '[deg] Same as {CRLT_OBS}',
        'DATE-OBS': r'[ISO] Observation date {DATE__OBS}',
        'RSUN_OBS': '[arcsec] angular radius of Sun',
        # 'CADENCE': '[seconds] Check {TRECSTEP}',
        # 'TRECSTEP': '[seconds] Time resolution',
    }

    for key, value in keywords.items():
        map_out.meta[key] = value

    # Manually add units to have proper map.spatial_units
    map_out.meta['CUNIT1'] = 'arcsec'
    map_out.meta['CUNIT2'] = 'arcsec'
    map_out.meta['CUNIT3'] = 's'
    map_out.meta['CTYPE3'] = 'TIME'

    map_out.plot_settings = map_input.plot_settings

    if doppler:
        map_out = correct_solar_rotation(map_out)

    if save:
        map_out = writefits(map_out.data, map_out.meta, name_output,
                            bitpix, bscale, bzero, blank, save=save)

    return Map(map_out.data, map_out.meta)


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            prog='los2los.py',
            description='Read anc crop a Line-Of-Sight (LOS).'
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
    map_out = los2los(
        map_input,
        crln=args.crln,
        crlt=args.crlt,
        doppler=args.dopp,
        name_output=args.output,
        continuum=args.int,
        save=True,
        algorithm=args.algorithm
    )
