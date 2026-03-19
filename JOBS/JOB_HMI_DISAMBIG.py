#!/usr/bin/env python
'''
Resolves the ambiguation of the azimuth component of the magnetic field for
the azimuth and disambig files specified in the first and second argument to
be written in the directory of the third argument.
'''


def process_pair(args):
    i, n_files, name_azi, name_dis, dir_output, method = args

    from sunpy.map import Map
    from sunpython.hmi import hmi_disambig
    from datetime import datetime
    # import os

    map_azi_i = Map(name_azi)
    map_dis_i = Map(name_dis)

    time_azi = map_azi_i.meta['T_REC']
    time_dis = map_dis_i.meta['T_REC']

    if time_azi != time_dis:
        print(f'Time mismatch: {name_azi}')
        return

    time_core = time_azi[:23].replace('T', '_')
    dt = datetime.strptime(time_core, '%Y.%m.%d_%H:%M:%S.%f')
    time_out = dt.strftime('%Y%m%d_%H%M%S')

    observable = 'azimuth_disambig'
    output_name = f'{dir_output}/hmi.b_720s.{time_out}.{observable}.fits'

    # print(f'Calculating disambiguation for {time_azi}. '
    #       f'File {i+1} of {n_files}')

    hmi_disambig(
        map_azi_i,
        map_dis_i,
        method=method,
        output=output_name,
        save=True
    )

    return i


def disambiguate(dir_azimuth, dir_disambig,
                 dir_output='./hmi_azimuth_disambig/', method=2,
                 nproc=None):

    import glob
    import os
    from multiprocessing import Pool, cpu_count

    os.makedirs(dir_output, exist_ok=True)

    list_frames_azi = sorted(glob.glob(f'{dir_azimuth}/hmi*.fits'))
    list_frames_dis = sorted(glob.glob(f'{dir_disambig}/hmi*.fits'))

    n_files = len(list_frames_azi)
    print(f'Number of files: {n_files}\n')

    if nproc is None:
        nproc = cpu_count()

    tasks = [
        (
            i, n_files, list_frames_azi[i], list_frames_dis[i], dir_output,
            method)
        for i in range(n_files)
    ]

    with Pool(nproc) as pool:
        for j, _ in enumerate(pool.imap(process_pair, tasks), 1):
            print(f'Calculating disambiguation for file {j} of {n_files}')


# ---------------------------------------------------------------------------
# End of code
# ---------------------------------------------------------------------------


if __name__ == '__main__':

    import argparse
    import sys
    import os

    # Bash checks
    def parse_args():

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            prog='JOB_HMI_DISAMBIG.py',
            description=(
                'Resolve the azimuth ambiguity of HMI magnetic field data.\n'
                'The script matches azimuth and disambiguation FITS files '
                'by time and writes disambiguated azimuth maps.'
            )
        )

        parser.add_argument(
            'dir_azimuth',
            type=str,
            help='Directory containing azimuth FITS files'
        )

        parser.add_argument(
            'dir_disambig',
            type=str,
            help='Directory containing disambiguation FITS files'
        )

        parser.add_argument(
            '-o', '--output',
            dest='dir_output',
            type=str,
            default='./hmi_azimuth_disambig/',
            help='Output directory (default: %(default)s)'
        )

        parser.add_argument(
            '-m', '--method',
            type=int,
            choices=[0, 1, 2],
            default=2,
            help=(
                'Bit used for disambiguation: (default: %(default)s)\n'
                '  0 = potential acute.\n'
                '  1 = random.\n'
                '  2 = radial acute.'
            )
        )

        parser.add_argument(
            '-n', '--nproc',
            type=int,
            default=None,
            help='Number of CPU cores to use'
        )

        return parser.parse_args()

    args = parse_args()

    # ------------------------------------------------------------------

    if not os.path.isdir(args.dir_azimuth):
        sys.exit(f'Error: directory does not exist: {args.dir_azimuth}')

    if not os.path.isdir(args.dir_disambig):
        sys.exit(f'Error: directory does not exist: {args.dir_disambig}')

    disambiguate(
        args.dir_azimuth,
        args.dir_disambig,
        dir_output=args.dir_output,
        method=args.method,
        nproc=args.nproc
    )

