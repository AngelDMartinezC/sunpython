#!/usr/bin/env python
'''
Code to convert from Postel to full-disk Helioprojective cartesian. The header
of the Postel map should include the info of the World Coordinate System (WCS)
and the T_REC keyword (same as the one specified in the filename).
The output is a 4096x4096 map
'''

# TODO: Update bitpix of output datacube to match the one of the input map.
import astropy.units as u
from astropy.io import fits


def _get_hmi_keywords(t_rec):
    from drms import Client
    '''
    Given a string date in the format YYYY.mm.dd_HH:MM:SS, extract keywords to
    use them as template in the helioprojective header. It returns a drms
    query response containting metadata info for reprojection
    '''

    # Selects initial and final time. Here I choose Doppler data but it is
    # only used for reading the header
    trecstep = 45  # TRECSTEP for HMI
    qstr = f'hmi.v_{int(trecstep)}s[{t_rec}_TAI]'

    # Get coordinate metadata information
    client = Client()
    query = client.query(qstr, key='T_REC, CRPIX1, CRPIX2, CDELT1, CDELT2,'
                         'CUNIT1, CUNIT2, CRVAL1, CRVAL2, CROTA2, CRLN_OBS,'
                         'CRLT_OBS, TELESCOP, DSUN_OBS, RSUN_OBS, DATE-OBS,'
                         'CAR_ROT, T_REC, T_OBS, WCSNAME',
                         convert_numeric=True)
    return query


def _update_map_header(map_out, keywords):
    '''
    Update the header of map_out with keywords from the Helioprojective
    template.

    Parameters
    ----------
    map_out : sunpy.map.Map
        Reprojected solar map header to be updated.

    keywords : drms.QueryResponse
        Query response containing metadata keywords.
    '''

    # Update header with relevant keywords
    for key in keywords:
        map_out.meta[key] = keywords[key].values[0]

    # Keyword comments
    map_out.meta['keycomments'] = {
        'CDELT1': '[arcsec/pixel] image scale in the x direction',
        'CDELT2': '[arcsec/pixel] image scale in the y direction',
        'CRPIX1': '[pixel] CRPIX1: location of the Sun center',
        'CRPIX2': '[pixel] CRPIX2: location of the Sun center',
        'CRVAL1': '[arcsec] CRVAL1: x origin',
        'CRVAL2': '[arcsec] CRVAL2: y origin',
        'CRLN_OBS': '[deg] Carrington longitude of the observer',
        'CRLT_OBS': '[deg] Carrington latitude of the observer',
        'CROTA2': '[deg]',
        'CAR_ROT': 'Carrington rotation number of CRLN_OB',
        'DSUN_OBS': '[m] Distance from SDO to Sun center',
        'RSUN_OBS': '[arcsec] angular radius of Sun',
        'T_OBS': '[TAI] Slot time',
        'T_REC': '[TAI] Nominal time',
        'DATE-OBS': r'[ISO] Observation date {DATE__OBS}',
    }

    return None


def _convert_to_helioprojective(map_input, naxis=4096, algorithm='exact'):
    from astropy.coordinates import SkyCoord
    from sunpy.map import make_fitswcs_header

    # # Calculate sun radius in meters
    # rsun_obs = map_input.rsun_obs.to(u.rad)
    # dsun = map_input.dsun
    # rsun = rsun_obs.value*dsun

    t_rec = map_input.date
    if t_rec is None:
        t_rec = map_input.meta['T_REC'][0:19]
    else:
        t_rec = str(t_rec)

    # Get keywords of interest
    crln_obs = map_input.meta['CRLN_OBS']
    crlt_obs = map_input.meta['CRLT_OBS']
    frame = map_input.coordinate_frame

    # Shift the center to align with the center of the reprojected map
    map_input.meta['CRPIX1'] -= 1
    map_input.meta['CRPIX2'] -= 1

    # Get keywords from the Helioprojective map
    keywords = _get_hmi_keywords(t_rec)
    naxis1 = naxis
    naxis2 = naxis
    crpix1 = float(keywords['CRPIX1'].values[0])
    crpix2 = float(keywords['CRPIX2'].values[0])
    cdelt1 = float(keywords['CDELT1'].values[0])
    cdelt2 = float(keywords['CDELT2'].values[0])
    crota2 = float(keywords['CROTA2'].values[0])

    # Retrieve data of center of map, then convert it to helioprojective
    origin_hgcar = SkyCoord(
        crln_obs,
        crlt_obs,
        unit=u.deg,
        frame=frame
    )
    origin = origin_hgcar.helioprojective

    # Create header of a helioprojective projection
    target_header = make_fitswcs_header(
        data=(naxis1, naxis2),  # size of output datacube
        scale=[cdelt1, cdelt2] * u.arcsec/u.pix,
        coordinate=origin,
        rotation_angle=crota2 * u.deg,
        projection_code="TAN",  # TAN: gnomonic projection
        reference_pixel=[crpix1, crpix2] * u.pix
    )

    # Reproject map from Postel into Helioprojective cartesian
    map_out = map_input.reproject_to(target_header, algorithm=algorithm,
                                     parallel=False)
    return map_out, keywords


def postel2los(input_map, naxis=4096, algorithm='interpolation',
               name_output='HP_map.fits', save=False):
    '''
    Read info of Postel map and reprojects it to Line Of Sight.

    Parameters
    ----------
    input_map : sunpy.map.GenericMap
        Postel projected map
    naxis : int
        NAXIS1 and NAXIS2 of the projected frame
    algorithm : {'interpolation', 'adaptive', 'exact'}
        Projection algorithm (from fastest to slowest).
    name_output : str
        Name of output file.

    Returns
    -------
    sunpy.map.Map
        Reprojected map
    '''

    map_out, keywords = _convert_to_helioprojective(input_map, naxis,
                                                    algorithm=algorithm)
    # Update header of reprojected map
    _update_map_header(map_out, keywords)
    out_hdu = fits.PrimaryHDU(data=map_out.data, header=map_out.fits_header)

    # TODO: Check BITPITX here
    # Change BITPIX of datacube
    # _check_and_change_bitpix(input_map, out_hdu)

    # Save file
    if save:
        out_hdu.writeto(name_output, overwrite=True)

    return out_hdu


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():
        parser = argparse.ArgumentParser(
            description=(
                "Read a Postel-projected map and reproject it back to "
                "Line-Of-Sight (LOS)."
            )
        )

        # Positional arguments
        parser.add_argument(
            "input",
            help="Input Postel FITS file"
        )

        parser.add_argument(
            "output",
            help="Output LOS FITS file"
        )

        # Optional arguments
        parser.add_argument(
            "--size",
            type=int,
            default=4096,
            help=(
                "X and Y dimensions (NAXIS1 = NAXIS2) of the output LOS map. "
                "Default: 4096"
            )
        )

        parser.add_argument(
            "--algorithm",
            choices=["interpolation", "adaptive", "exact"],
            default="interpolation",
            help=(
                "Reprojection algorithm (from fastest to slowest). "
                "Default: interpolation"
            )
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    from sunpy.map import Map

    # Read input map
    map_input = Map(args.input)

    # Run reprojection
    map_out = postel2los(
        input_map=map_input,
        naxis=args.size,
        algorithm=args.algorithm,
        name_output=args.output,
        save=True
    )
