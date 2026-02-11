#!/usr/bin/env python
'''
Processes a set of solar FITS files and projects their data onto a
line-of-sight (LOS) map centered at a specified Carrington longitude and
latitude. The final output is a single stacked FITS datacube.
'''


def _los_worker(args):
    '''
    Worker function for parallel Postel projection.
    '''
    (
        i, name_input, frames_dir,
        crln, crlt, doppler, continuum, algorithm
    ) = args

    from sunpy.map import Map
    from sunpython.map_projection import los2los
    import os

    name_output_i = os.path.join(frames_dir, f'FRAME.{i:04d}.fits')

    input_map = Map(name_input)
    los2los(
        input_map,
        crln, crlt,
        name_output=name_output_i,
        doppler=doppler,
        continuum=continuum,
        save=True,
        algorithm=algorithm
    )

    return name_output_i


def los(directory, crln, crlt, name_output, doppler=False,
        continuum=False, prefix=None, nproc=1, algorithm='interpolation'):

    print(
        f'\n'
        f'Initial parameters:\n'
        f'    Projection:            Helioprojective-Cartesian\n'
        f'    CRLN:                  {crln} deg\n'
        f'    CRLT:                  {crlt} deg\n'
        f'    OUTPUT:                {name_output}\n'
        f'    algorithm:             {algorithm}\n'
        f'    doppler correction:    {doppler}\n'
        f'    limb dark correction:  {continuum}\n'
        f'    file prefix:           {prefix}\n'
    )

    # from sunpython.map_projection import los2los
    # from sunpy.map import Map
    # import multiprocessing
    import glob
    import subprocess
    import os
    import shutil
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from multiprocessing import set_start_method  # , Pool
    from astropy.io import fits

    set_start_method("spawn")

    # Write FRAMES directory to store temporal frames
    frames_dir = os.path.join(directory, 'FRAMES')

    # Cleanup to empty dir
    if os.path.exists(frames_dir):
        shutil.rmtree(frames_dir)

    os.makedirs(frames_dir)

    # Prefix name for files. By default it reads every fits file in dir
    if prefix is not None:
        list_frames = sorted(glob.glob(f'{directory}/{prefix}*.fits'))
    else:
        list_frames = sorted(glob.glob(f'{directory}/*.fits'))
    if not list_frames:
        print('No fits files found.')

    n_files = len(list_frames)
    print(f'Number of files: {n_files}\n')

    # # Old serial version
    # for i in range(n_files):
    #     # Read input file and assign an output name
    #     name_input = list_frames[i]
    #     print(f'Projection of {name_input}. File {i+1} of {n_files}')
    #     name_output_i = os.path.join(frames_dir, f'FRAME.{i:04d}.fits')
    #     input_map = Map(name_input)
    #     map_postel = los2los(
    #         input_map,
    #         crln, crlt,
    #         name_output=name_output_i,
    #         doppler=doppler,
    #         continuum=continuum,
    #         save=True,
    #     )

    # n_workers = min(multiprocessing.cpu_count(), n_files)
    n_workers = min(nproc, n_files)
    print(f'Running Postel projection using {n_workers} processes\n')
    tasks = [(
        i, list_frames[i], frames_dir,
        crln, crlt, doppler, continuum, algorithm
        )
        for i in range(n_files)
    ]

    # If nproc is 1, run sequentially
    if n_workers == 1:
        for i, fname in enumerate(list_frames):
            _los_worker((
                i, fname, frames_dir,
                crln, crlt, doppler, continuum, algorithm
            ))
            print(f'Projection of {fname}. File {i+1} of {n_files}')
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [executor.submit(_los_worker, t) for t in tasks]

            for k, future in enumerate(as_completed(futures), start=1):
                future.result()  # re-raise exceptions
                print(
                    f'Projection of {list_frames[k-1]}. File {k} of {n_files}'
                )

    # Read header from first projected frame
    first_frame = os.path.join(frames_dir, 'FRAME.0000.fits')
    hdr = fits.getheader(first_frame)

    if 'TRECSTEP' in hdr:
        trecstep = hdr['TRECSTEP']
    elif 'CADENCE' in hdr:
        trecstep = hdr['CADENCE']
    else:
        import numpy as np
        from astropy.time import Time
        from astropy.io import fits

        times = Time([fits.getheader(f)['DATE-OBS'] for f in list_frames])
        trecstep = np.median(np.diff(times.jd)) * 86400  # seconds

    num_frames = len(glob.glob(f'{frames_dir}/FRAME*.fits'))

    with open(f'{directory}/list.txt', 'w') as f:
        for i in range(num_frames):
            f.write(f'{frames_dir}/FRAME.{i:04d}.fits\n')

    if os.path.exists(name_output):
        os.remove(name_output)
        # print(f'File '{name_output}' has been removed.')

    # TODO! Implement stack_frames
    subprocess.run(f'stack_frames {name_output} -l {directory}/list.txt',
                   check=True, shell=True)

    # Add CDELT3 keyword
    subprocess.run(f'revise_head {name_output} -k CDELT3 -d',
                   check=True, shell=True)
    subprocess.run(f'revise_head {name_output} -k CDELT3 -v {trecstep} '
                   '-l 11 -c "/ [s]"', check=True, shell=True)

    # Add DAXIS3 keyword
    subprocess.run(f'revise_head {name_output} -k DAXIS3 -d',
                   check=True, shell=True)
    subprocess.run(f'revise_head {name_output} -k DAXIS3 -v {trecstep} '
                   '-l 12 -c "/ [s]"', check=True, shell=True)

    if 'TRECSTEP' not in hdr:
        subprocess.run(f'revise_head {name_output} -k TRECSTEP -v '
                       f'{trecstep} -l 16 -c "/ [s]"',
                       check=True, shell=True)

    print(f'File {name_output} written')
    # Cleanup
    if os.path.exists(frames_dir):
        shutil.rmtree(frames_dir)


if __name__ == '__main__':

    import argparse

    def parse_args():
        parser = argparse.ArgumentParser(
            prog='JOB_LOS.py',
            description='Crops a LOS map.'
        )

        parser.add_argument(
            'path',
            help='Input path (directory or FITS file)'
        )

        parser.add_argument(
            'crln',
            type=float,
            help='Carrington longitude of projection center (deg)'
        )

        parser.add_argument(
            'crlt',
            type=float,
            help='Carrington latitude of projection center (deg)'
        )

        parser.add_argument(
            'output',
            help='Output FITS filename'
        )

        parser.add_argument(
            '--algorithm',
            default='interpolation',
            choices=['interpolation', 'adaptive', 'exact'],
            help=(
                'Accuracy of the projected data. From fastest to slowest: '
                'interpolation, adaptive, exact (default: %(default)s).'
            )
        )

        parser.add_argument(
            '--nproc',
            type=int,
            default=1,
            help='Number of parallel processes (default: serieal).'
        )

        group = parser.add_mutually_exclusive_group()

        group.add_argument(
            '--dopp',
            action='store_true',
            help='Enable Doppler-related processing (default: False)'
        )

        group.add_argument(
            '--int',
            action='store_true',
            help='Enable limb darkening correction (default: False)'
        )

        parser.add_argument(
            '--prefix',
            type=str,
            default=None,
            help='Prefix of fits files to be projected'
        )

        return parser.parse_args()

    # ------------------------------------------------------------------

    args = parse_args()

    los(
        directory=args.path,
        crln=args.crln,
        crlt=args.crlt,
        name_output=args.output,
        doppler=args.dopp,
        continuum=args.int,
        prefix=args.prefix,
        nproc=args.nproc,
        algorithm=args.algorithm
    )
