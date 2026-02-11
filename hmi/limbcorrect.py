#!/usr/bin/env python
"""
Program to correct limb darkening in python based on ILD/swwidl
Based on:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro

Coefficients taken from
    Cox, A. N.: Allen's Astrophysical Quantities, Springer, 2000 taken from IDL

Author:  AngelMartinezC
"""


# TODO! Make the same limbcorrect process outputs a BLANK
def limbcorrect(map_input, name_output=None, save=False):
    """
    Correct the limb darkening of an intensity continuum observation

    Parameters
    ----------
    map_input : astropy.io.fits.HDUList or sunpy.map.GenericMap
        Intensity map to correct
    name_output : str
        Name of output file.
    save : bool
        Write the fits file

    Returns
    -------
    sunpy.map.Map
        Corrected map
    """

    from sunpython.fits_work import load_map, writefits
    from sunpython.coordinate_conversion import hpxy2xy
    import numpy as np

    # -------------------------------------------------------------

    def darklimb_u(ll):
        pll = np.array([1.0, ll, ll**2, ll**3, ll**4, ll**5])
        au = -8.9829751
        bu = 0.0069093916
        cu = -1.8144591e-6
        du = 2.2540875e-10
        eu = -1.3389747e-14
        fu = 3.0453572e-19
        a = np.array([au, bu, cu, du, eu, fu])
        ul = sum(a*pll)
        return ul

    def darklimb_v(ll):
        pll = np.array([1.0, ll, ll**2, ll**3, ll**4, ll**5])
        av = 9.2891180
        bv = -0.0062212632
        cv = 1.5788029e-6
        dv = -1.9359644e-10
        ev = 1.1444469e-14
        fv = -2.599494e-19
        a = np.array([av, bv, cv, dv, ev, fv])
        vl = sum(a*pll)
        return vl

    # -------------------------------------------------------------

    # Read data files
    data, header, _, bscale, bzero, blank = load_map(map_input,
                                                     get_scale=True)

    # Parameters needed for the function
    wavelnth = header['WAVELNTH']  # Wavelenght
    # crval1 = header['CRVAL1']  # X center in arcsec
    # crval2 = header['CRVAL2']  # Y center in arcsec
    radius = header['RSUN_OBS']/header['CDELT1']  # Sun radius in pixels
    naxis1 = header['NAXIS1']  # X array size
    naxis2 = header['NAXIS2']  # X array size

    ll = 1.0*wavelnth
    crpix1, crpix2 = hpxy2xy(map_input, 0, 0)
    crpix1 = crpix1.value
    crpix2 = crpix2.value

    # NaNs = np.where(data < 1000)
    # print(NaNs)
    # data[NaNs] = -10000  # np.nan  # Make zero all NANs

    # Apply correction
    u_l = darklimb_u(ll)
    v_l = darklimb_v(ll)

    x_arr = np.arange(0, naxis1, 1)  # Make x array
    y_arr = np.arange(0, naxis2, 1)  # Make y array
    x_mesh, y_mesh = np.meshgrid(y_arr, x_arr)  # Make xy array

    # Make a circle centered at crpix(1, 2) with radius the Sun radius
    circle = np.sqrt((x_mesh - crpix1)**2 + (y_mesh - crpix2)**2)

    # Normalize the circle so that inside radius is the unity
    grid = circle/radius

    # Remove the outside domain of the circle
    out = np.where(grid >= 1.0)
    grid[out] = np.nan

    # Apply correction
    limbfilt = (1.0 - u_l - v_l + u_l*np.cos(np.arcsin(grid))
                + v_l*np.cos(np.arcsin(grid))**2)

    # Final image
    imgout = data/limbfilt

    # HACK! The value -1000 seems to work to produce blanks instead of NaNs.
    # This value is purely empirical.
    # HACK! Here I'm writing twice depending on wheather I want to save the
    # data or I want to use it -and therefore read blanks-.
    if save:
        imgout[out] = -1000  # Make zero outside arcsin domain
        map_out = writefits(imgout, header, name_output, 32, bscale, bzero,
                            blank, save=save)
        imgout[out] = np.nan
    else:
        imgout[out] = np.nan
        map_out = writefits(imgout, header, name_output, 32, bscale, bzero,
                            blank, save=False)

    return map_out


if __name__ == '__main__':

    import sys
    # Check number of arguments
    log_msg = \
        "Usage: limbcorrect.py input[FITS] output[FITS]"

    if len(sys.argv) != 3:
        sys.tracebacklimit = 0
        raise TypeError(log_msg)

    # Read input file and assign an output name
    name_input = sys.argv[1]
    name_output = sys.argv[2]

    from astropy.io import fits
    map_input = fits.open(name_input)

    limbcorrect(map_input, name_output, save=True)
