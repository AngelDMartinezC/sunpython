#!/usr/bin/env python
'''
Resolves the ambiguation of the azimuth component of the magnetic field for
the azimuth and disambig files specified in the first and second argument to
be written in the directory of the third argument.
'''


def process_frame(args):
    (
        i, _, name_field, name_incl, name_azi_dis, dir_bp, dir_bt,
        dir_br) = args

    from sunpy.map import Map
    from sunpython.hmi import hmi_b2ptr
    from datetime import datetime

    map_field = Map(name_field)
    map_incl = Map(name_incl)
    map_azi_dis = Map(name_azi_dis)

    # Fix metadata for other instruments
    if 'T_REC' not in map_field.meta:
        map_field.meta['T_REC'] = map_field.meta['DATE-OBS']
        map_incl.meta['T_REC'] = map_incl.meta['DATE-OBS']
        map_azi_dis.meta['T_REC'] = map_azi_dis.meta['DATE-OBS']

    time_field = map_field.meta['T_REC']
    time_incl = map_incl.meta['T_REC']
    time_azi_dis = map_azi_dis.meta['T_REC']

    if time_field != time_incl or time_field != time_azi_dis:
        return i

    time_core = time_field[:23].replace('T', '_')
    dt = datetime.strptime(time_core, '%Y.%m.%d_%H:%M:%S.%f')
    time_out = dt.strftime('%Y%m%d_%H%M%S')

    name_bp = f'{dir_bp}/map.b_720s.{time_out}.Bp.fits'
    name_bt = f'{dir_bt}/map.b_720s.{time_out}.Bt.fits'
    name_br = f'{dir_br}/map.b_720s.{time_out}.Br.fits'

    hmi_b2ptr(
        map_field,
        map_incl,
        map_azi_dis,
        bp_output=name_bp,
        bt_output=name_bt,
        br_output=name_br,
        save=True
    )

    return i


def b2ptr(dir_field, dir_inclination, dir_azimuth_dis,
          dir_bp='./hmi_Bp', dir_bt='./hmi_Bt', dir_br='./hmi_Br',
          nproc=None):

    import glob
    import os
    from multiprocessing import Pool, cpu_count

    os.makedirs(dir_bp, exist_ok=True)
    os.makedirs(dir_bt, exist_ok=True)
    os.makedirs(dir_br, exist_ok=True)

    list_frames_field = sorted(glob.glob(f'{dir_field}/*.fits'))
    list_frames_incl = sorted(glob.glob(f'{dir_inclination}/*.fits'))
    list_frames_azi_dis = sorted(glob.glob(f'{dir_azimuth_dis}/*.fits'))

    if (
        len(list_frames_field) != len(list_frames_incl)
        or len(list_frames_field) != len(list_frames_azi_dis)
    ):
        raise AssertionError(
            'Directories do not contain the same number of files.'
        )

    n_files = len(list_frames_field)
    print(f'Number of files: {n_files}\n')

    if nproc is None:
        nproc = cpu_count()

    tasks = [
        (
            i,
            n_files,
            list_frames_field[i],
            list_frames_incl[i],
            list_frames_azi_dis[i],
            dir_bp,
            dir_bt,
            dir_br
        )
        for i in range(n_files)
    ]

    with Pool(nproc) as pool:
        for j, _ in enumerate(pool.imap(process_frame, tasks), 1):
            print(f'Calculating Bptr for file {j} of {n_files}')


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    import argparse

    # Bash checks
    def parse_args():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=(
                'Resolve the azimuth ambiguity of HMI magnetic field data and '
                'compute Bp, Bt, Br components.'
            )
        )

        # Required positional arguments
        parser.add_argument(
            'dir_field',
            help='Directory containing magnetic field strength FITS files'
        )
        parser.add_argument(
            'dir_inclination',
            help='Directory containing inclination FITS files'
        )
        parser.add_argument(
            'dir_azimuth_dis',
            help='Directory containing disambiguated azimuth FITS files'
        )

        # Optional arguments
        parser.add_argument(
            '--dir_bp',
            default='./hmi_Bp',
            help='Output directory for Bp maps (default: ./hmi_Bp)'
        )
        parser.add_argument(
            '--dir_bt',
            default='./hmi_Bt',
            help='Output directory for Bt maps (default: ./hmi_Bt)'
        )
        parser.add_argument(
            '--dir_br',
            default='./hmi_Br',
            help='Output directory for Br maps (default: ./hmi_Br)'
        )

        parser.add_argument(
            '-n', '--nproc',
            type=int,
            default=None,
            help='Number of CPU cores to use'
        )

        return parser.parse_args()

    # ------------------------------------------------------------------

    args = parse_args()

    b2ptr(
        args.dir_field,
        args.dir_inclination,
        args.dir_azimuth_dis,
        dir_bp=args.dir_bp,
        dir_bt=args.dir_bt,
        dir_br=args.dir_br,
        nproc=args.nproc
    )
