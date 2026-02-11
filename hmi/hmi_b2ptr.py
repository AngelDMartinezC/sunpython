#!/usr/bin/env python
'''
Code to compute the spherical phi, theta and r components of the magnetic
field from B(inclination, azimuth_dis)
'''


def hmi_b2ptr(field, inclination, azimuth_dis, bp_output='hmi.Bp_720s.fits',
              bt_output='hmi.Bt_720s.fits', br_output='hmi.Br_720s.fits',
              save=False):
    '''
    Convert HMI vector field from native components (field, inclination,
    azimuth_dis) into spherical components (Bp, Bt, Br).

    Parameters
    ----------
    field : sunpy.map.GenericMap
        field map in Gauss
    inclination : sunpy.map.GenericMap
        inclination map in deg
    azimuth_dis : sunpy.map.GenericMap
        azimuth map in deg
    save : bool
        Write the fits file

    Returns
    -------
    bptr : ndarray, shape (3, nx, ny)
        (Bp, Bt, Br) in Gauss
    '''

    from sunpy.map import Map
    import numpy as np
    import astropy.units as u
    from sunpy.coordinates import frames
    from sunpython.fits_work import load_map, writefits

    field_data, header, _, bscale, bzero, blank = load_map(
        field, get_scale=True, astropy_header=True)

    # Check if shapes are the same for the maps
    if (
        field_data.shape != inclination.data.shape
        or field_data.shape != azimuth_dis.data.shape
    ):
        raise AssertionError('The sizes are not the same')

    b_vec = np.array([field_data, inclination.data, azimuth_dis.data])

    # Get dimmensions
    _, nx, ny = b_vec.shape

    # ----------------------------
    # Convert to B_xi, B_eta, B_zeta
    # Eq. (1) of Sun (2013)
    # ----------------------------
    b_abs = b_vec[0, :, :]
    gamma = np.deg2rad(b_vec[1, :, :])  # inclination
    psi = np.deg2rad(b_vec[2, :, :])   # azimuth_dis

    b_xi = -b_abs*np.sin(gamma)*np.sin(psi)
    b_eta = b_abs*np.sin(gamma)*np.cos(psi)
    b_zeta = b_abs*np.cos(gamma)

    # ----------------------------
    # WCS to Stonyhurst lon/lat
    # ----------------------------
    # Create a dummy map to access SunPy coordinate machinery
    dummy_data = np.zeros((nx, ny))
    smap = Map(dummy_data, header)

    # meshgrid order chosen to match SunPy pixel convention for HMI maps
    hpc = smap.pixel_to_world(
        *np.meshgrid(np.arange(ny), np.arange(nx))*u.pix)

    hgs = hpc.transform_to(frames.HeliographicStonyhurst)

    phi = hgs.lon.to(u.rad).value  # longitude
    lam = hgs.lat.to(u.rad).value  # latitude

    # if return_lonlat:
    #     lonlat = np.zeros((2, nx, ny))
    #     lonlat[0, :, :] = np.rad2deg(phi)
    #     lonlat[1, :, :] = np.rad2deg(lam)

    # ----------------------------
    # Angles from header
    # ----------------------------
    b = np.deg2rad(header['CRLT_OBS'])   # disk-center latitude
    p = -np.deg2rad(header['CROTA2'])    # p-angle (note minus!)

    sinb, cosb = np.sin(b), np.cos(b)
    sinp, cosp = np.sin(p), np.cos(p)

    sinphi, cosphi = np.sin(phi), np.cos(phi)
    sinlam, coslam = np.sin(lam), np.cos(lam)

    # ----------------------------
    # Transformation matrix
    # Eq. (7)(8) in Sun (2013)
    # ----------------------------
    k11 = coslam*(sinb*sinp*cosphi + cosp*sinphi) - sinlam*cosb*sinp
    k12 = -coslam*(sinb*cosp*cosphi - sinp*sinphi) + sinlam*cosb*cosp
    k13 = coslam*cosb*cosphi + sinlam*sinb

    k21 = sinlam*(sinb*sinp*cosphi + cosp*sinphi) + coslam*cosb*sinp
    k22 = -sinlam*(sinb*cosp*cosphi - sinp*sinphi) - coslam*cosb*cosp
    k23 = sinlam*cosb*cosphi - coslam*sinb

    k31 = -sinb*sinp*sinphi + cosp*cosphi
    k32 = sinb*cosp*sinphi + sinp*cosphi
    k33 = -cosb*sinphi

    # ----------------------------
    # Output: (Bp, Bt, Br)
    # ----------------------------
    # bptr = np.zeros((3, nx, ny))
    # bptr[0, :, :] = k31*b_xi + k32*b_eta + k33*b_zeta  # Bp (west)
    # bptr[1, :, :] = k21*b_xi + k22*b_eta + k23*b_zeta  # Bt (south)
    # bptr[2, :, :] = k11*b_xi + k12*b_eta + k13*b_zeta  # Br (radial)

    # Here we do: Bptr[i, x ,y] = K[i, j, x, y] * B[j, x, y]
    # Order: (Bp, Bt, Br) = (west, south, radial)
    B = np.stack((b_xi, b_eta, b_zeta), axis=0)
    K = np.stack((
        np.stack((k31, k32, k33), axis=0),
        np.stack((k21, k22, k23), axis=0),
        np.stack((k11, k12, k13), axis=0),
    ), axis=0)

    # Einstein summation notation
    bptr = np.einsum('ijxy,jxy->ixy', K, B)

    nan_mask = ~np.isfinite(bptr)
    bptr[nan_mask] = blank

    # Set to blank large values near the limb. The value 1e5 is empirical
    # TODO! This threshold is empirical and should be validated.
    LIMIT_THRESH = 1e4
    bptr[np.abs(bptr) > LIMIT_THRESH] = blank

    # HACK!! This function will probably be changing. Better implementation
    # required.
    # Save Bp, Bt, and Br
    writefits(bptr[0], header, bp_output, 64, bscale, bzero, blank, save=save)
    writefits(bptr[1], header, bt_output, 64, bscale, bzero, blank, save=save)
    writefits(bptr[2], header, br_output, 64, bscale, bzero, blank, save=save)

    # if return_lonlat:
    #     return bptr, lonlat

    return bptr


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    def parse_args():
        parser = argparse.ArgumentParser(
            prog='hmi_b2ptr.py',
            description=(
                'Compute spherical magnetic-field components (Bp, Bt, Br) '
                'from HMI field, inclination, and azimuth maps.\n\n'
                'Note: azimuth must be the disambiguated one.'
            )
        )

        # Required inputs
        parser.add_argument(
            'field',
            type=str,
            help='Magnetic field strength FITS file'
        )

        parser.add_argument(
            'inclination',
            type=str,
            help='Inclination FITS file (degrees)'
        )

        parser.add_argument(
            'azimuth_disambig',
            type=str,
            help='Disambiguated azimuth FITS file (degrees)'
        )

        # Optional outputs
        parser.add_argument(
            '--bp-output',
            default='hmi.Bp_720s.fits',
            help='Output filename for Bp (default: %(default)s)'
        )

        parser.add_argument(
            '--bt-output',
            default='hmi.Bt_720s.fits',
            help='Output filename for Bt (default: %(default)s)'
        )

        parser.add_argument(
            '--br-output',
            default='hmi.Br_720s.fits',
            help='Output filename for Br (default: %(default)s)'
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    import os
    import sys
    from sunpy.map import Map

    # ---- sanity checks ----
    for fname, label in [
        (args.field, 'field'),
        (args.inclination, 'inclination'),
        (args.azimuth_disambig, 'azimuth'),
    ]:
        if not os.path.isfile(fname):
            sys.exit(f'Error: {label} file "{fname}" does not exist.')

    # ---- load maps ----
    map_field = Map(args.field)
    map_inclination = Map(args.inclination)
    map_azimuth = Map(args.azimuth_disambig)

    # ---- compute Bp, Bt, Br ----
    bptr = hmi_b2ptr(
        map_field,
        map_inclination,
        map_azimuth,
        bp_output=args.bp_output,
        bt_output=args.bt_output,
        br_output=args.br_output,
        save=True
    )
