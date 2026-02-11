#!/usr/bin/env python

"""
Manage fits files using python
"""
import numpy as np

# Bitpix dynamic range constants
DYNAM_RANGE_BYTE = 255
DYNAM_RANGE_SHORT = 32767
DYNAM_RANGE_INT = 2147483647
DYNAM_RANGE_LL = 9223372036854775807

# Factor introduced to properly scale data and avoid overflow
FACTOR = 1


def _get_scale_keys(header):
    # Get scale keywords before data memory storage
    bitpix = header.get('BITPIX')
    bscale = header.get('BSCALE', None)
    bzero = header.get('BZERO', None)
    blank = header.get('BLANK', None)
    return bitpix, bscale, bzero, blank


# Return the maximum representable dynamic range for a given BITPIX value.
def dynamic_range(bitpix):
    if bitpix not in (8, 16, 32, 64, -32, -64, 128, -128):
        raise ValueError(f"Unsupported BITPIX value: {bitpix}")
    if bitpix == 8:
        dynam_range_type = (float(DYNAM_RANGE_BYTE), np.int8)
    elif bitpix == 16:
        dynam_range_type = (float(DYNAM_RANGE_SHORT), np.int16)
    elif abs(bitpix) == 32:
        dynam_range_type = (float(DYNAM_RANGE_INT), np.int32)
    elif abs(bitpix) == 64:
        dynam_range_type = (float(DYNAM_RANGE_LL), np.int64)
    elif abs(bitpix) == 128:
        dynam_range_type = (float(DYNAM_RANGE_LL), np.float64)
    dynam_range = -1*int(dynam_range_type[0] + 1)
    return dynam_range, dynam_range_type[1]


# Convert sunpy meta header into astropy header for proper scaling management
def sunpy_meta_to_fits_header(meta):
    from astropy.io import fits
    from astropy.io.fits.header import _HeaderCommentaryCards

    header = fits.Header()

    for key, value in meta.items():
        k = str(key).upper()

        # Blank key or comment-only section
        if k == "" or k.startswith("/"):
            if isinstance(value, _HeaderCommentaryCards):
                for line in value:
                    header.append(("COMMENT", line))
            else:
                header.append(("COMMENT", str(value)))
            continue

        # Skip keycomments dict
        if k == "KEYCOMMENTS" and isinstance(value, dict):
            continue

        # HISTORY / COMMENT with newlines
        if k in ("HISTORY", "COMMENT"):
            if isinstance(value, str):
                for line in value.splitlines():
                    header.append((k, line))
            elif isinstance(value, (list, tuple, _HeaderCommentaryCards)):
                for line in value:
                    header.append((k, str(line)))
            else:
                header.append((k, str(value)))
            continue

        # Only add if FITS accepts it
        try:
            fits.Card(k, value)
        except Exception:
            continue

        header[k] = value

    # Apply keycomments if present
    kc = meta.get("keycomments") or meta.get("KEYCOMMENTS")
    if isinstance(kc, dict):
        for kk, com in kc.items():
            if kk in header:
                header.comments[kk] = com

    return header


def load_map(map_input, bitpix=None, get_scale=False):
    '''
    Load data and metadata from a SunPy GenericMap, handling FITS scaling
    keywords and NaN values

    This function extracts the data array and header from the input map.
    If the BITPIX is of type INT, then NaN values in the data are replaced by
    the physical value corresponding to the FITS BLANK keyword.

    Optionally, the FITS scaling parameters (BITPIX, BSCALE, BZERO, BLANK)
    can be returned or recomputed.

    Parameters
    ----------
        map_input : sunpy.map.GenericMap
            Input FITS HDUList or SunPy map object from which to extract the
            data and header.
        get_scale : bool, optional
            If True, return the original FITS scaling keywords associated with
            the input data. Default is False.

    Returns
    -------
        data : ndarray
            Data array extracted from the input map.
        header : dict-like
            SunPy metadata associated with the data.
        bitpix : int, optional
            FITS BITPIX value. Returned if `get_scale` is True.
        bscale : float, optional
            FITS BSCALE value. Returned if `get_scale` is True.
        bzero : float, optional
            FITS BZERO value. Returned if `get_scale` is True.
        blank : int or float, optional
            FITS BLANK value. Returned if `get_scale` is True.
    '''

    meta = map_input.meta

    bitpix0, bscale0, bzero0, blank0 = _get_scale_keys(meta)

    # Condition added for Solar Orbiter's unscaled BITPIX keyword
    if bitpix is None:
        if 'FILENAME' in meta:
            if (
                ('icnt' in meta['FILENAME']) or
                ('vlos' in meta['FILENAME']) or
                ('bazi' in meta['FILENAME']) or
                ('binc' in meta['FILENAME']) or
                ('euv174' in meta['FILENAME'])
            ):
                bitpix = 16
            elif (
                ('blos' in meta['FILENAME']) or
                ('bmag' in meta['FILENAME'])
            ):
                bitpix = 32
            else:
                bitpix = 32
        else:
            bitpix = bitpix0

    blank_dyn, _ = dynamic_range(bitpix)

    data = map_input.data

    # CHECK: Scale data
    if 'BUNIT' in meta:
        if meta['BUNIT'] == 'DN':
            exptime = meta.get('EXPTIME', meta.get('XPOSURE', None))
            if exptime is not None:
                data = data/exptime
    if 'PIXLUNIT' in meta:
        if meta['PIXLUNIT'] == 'DN':
            exptime = meta.get('EXPTIME', meta.get('XPOSURE', None))
            if exptime is not None:
                data = data/exptime

    if abs(bitpix) in [8, 16, 32]:
        if blank0 is None:
            blank = blank_dyn
        else:
            blank = blank0

        if (bscale0 is None) or (bzero0 is None):
            # Get minmax range
            vmin0 = np.nanmin(data)
            vmax0 = np.nanmax(data)
            delta = (vmax0 - vmin0) * FACTOR
            vmin = np.nanmin(data - delta)
            vmax = np.nanmax(data + delta)
            # Calculate bscale, bzero
            imin, imax = -abs(blank), abs(blank)
            bscale = (vmax - vmin) / (imax - imin)
            bzero = vmin - imin * bscale
        else:
            bzero = bzero0
            bscale = bscale0

    else:
        blank = None
        bscale = None
        bzero = None

    if get_scale:
        return data, meta, bitpix, bscale, bzero, blank
    else:
        return data, meta


def writefits(data, meta, output_name, bitpix=64, bscale=1, bzero=0,
              blank=None, save=True):
    '''
    Save FITS files with an almost proper suitable header managing blanks and
    BSCALES. This results in size with smaller sizes (as much as half)

    Parameters
    ----------
        data : ndarray
            Input MxN shape data
        meta : dict-like
            SunPy metadata associated with the data.
        output_name : str
            Name of output map
        bitpix: int
            target BITPIX
        bscale: float
            BSCALE
        bzero: float
            BZERO
        blank: int
            Large numerical value to be interpreted as BLANK

    Returns
    -------
        Map : sunpy.map.GenericMap
            Output FITS HDUList or SunPy map object
    '''
    from astropy.io.fits import PrimaryHDU
    from sunpy.map import Map
    import warnings

    # Check for apparently invalid cast an silence the warningg
    warnings.filterwarnings(
        "ignore",
        category=RuntimeWarning,
        message="invalid value encountered in cast"
    )

    blank_dyn, dtype = dynamic_range(bitpix)

    if abs(bitpix) in [8, 16, 32]:
        if blank is None:
            blank = blank_dyn

        if (bscale is None) or (bzero is None):
            # Get minmax range
            vmin0 = np.nanmin(data)
            vmax0 = np.nanmax(data)
            delta = (vmax0 - vmin0) * FACTOR
            vmin = np.nanmin(data - delta)
            vmax = np.nanmax(data + delta)
            # Calculate bscale, bzero
            imin, imax = -abs(blank), abs(blank)
            bscale = (vmax - vmin) / (imax - imin)
            bzero = vmin - imin * bscale
        # elif (bscale < 1e-4):
        #     vmin = np.nanmin(data)
        #     vmax = np.nanmax(data)
        #     imin, imax = -abs(blank), abs(blank)
        #     bscale = 0.001
        #     bzero = vmin - imin * bscale

        data_scaled = np.rint((data - bzero) / bscale).astype(dtype)
        mask = ~np.isfinite(data)
        data_scaled[mask] = blank
    else:
        data_scaled = data
        blank = None
        bscale = None
        bzero = None

    meta = meta.copy()

    header_tmp = sunpy_meta_to_fits_header(meta)
    hdu = PrimaryHDU(data_scaled, header_tmp)

    header = hdu.header

    if abs(bitpix) in [8, 16, 32]:
        header['BSCALE'] = bscale
        header['BZERO'] = bzero
        header['BLANK'] = blank

    if save:
        hdu.writeto(output_name, overwrite=True)

    return Map(data, meta)
