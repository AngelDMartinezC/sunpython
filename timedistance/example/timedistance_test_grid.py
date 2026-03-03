from sunpython import TimeDistance as td
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm

# Import all sunpy colormaps and register them with matplotlib
cm

# Specify the file and the depth of the acoustic source
solar_radius = 696  # Mm
depth = 1000/(solar_radius*1000)

name = './AIA_0171.fits'
x0 = 269  # x-center
y0 = 496  # y-center
rad0 = 0  # Skip these pixels in the TD plot
rad1 = 57  # Final radius (distance) in pixels
time0 = 0  # Initial frame for calculations
time1 = 116  # 200  # Final frame for calculations
frame0 = None  # frame at which the acoustic progression begins
th0 = 150
th1 = 155

vmin = 100
vmax = 800


def grid(d_angle=12):
    th0 = 0
    count = 1
    plt.figure(figsize=(13, 8))
    for i in range(30):
        th1 = th0 + d_angle
        data = td(name, x0=x0, y0=y0, theta0=th0, theta1=th1, radius0=rad0,
                  radius=rad1, time0=time0, time1=time1)
        ax = plt.subplot(5, 6, count)
        _, _ = data.plot(colorbar=False, vmin=vmin, vmax=vmax,
                         interpolation='sinc', get_cfa=True,
                         cmap='sdoaia171')
        plt.title(' ')
        plt.xlabel(' ')
        plt.ylabel(' ')
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.text(
            0.05, 0.95,
            (f'{th0}\N{degree sign} - {th1}\N{degree sign}'),
            transform=ax.transAxes, ha='left', va='top', color='black',
            bbox=dict(boxstyle="round", ec='black', fc='white', alpha=0.85)
        )

        print(f'{count}: {th0}-{th1}')
        plt.grid(ls='--', c='white', alpha=0.35)
        count += 1
        th0 = th1
    plt.suptitle(f'{name}. Coords (x0, y0): ({x0}, {y0})')
    plt.subplots_adjust(wspace=0.03, hspace=0.03, left=0.04, right=0.97,
                        bottom=0.05, top=0.94)
    plt.show()
    exit()


grid(d_angle=12)


datatd = td(name, x0=x0, y0=y0, theta0=th0, theta1=th1, radius0=rad0,
            radius=rad1, time0=time0, time1=time1)

daxis1 = datatd.get_daxis
daxis3 = datatd.get_daxis3
fig = plt.figure(figsize=(6, 5))
ax = plt.subplot(111)
datatd.plot(vmin=vmin, vmax=vmax, cmap='sdoaia171', interpolation='sinc')
# plt.colorbar(ax.images[0])
plt.show()
exit()


def save(name='TD_diagram.fits', sq_int=0):
    from astropy.io import fits
    image, im = datatd.plot(colorbar=False, interpolation='sinc')
    hdu = fits.PrimaryHDU(image)
    hdr = hdu.header
    # play and change some of these variables
    hdr['CDELT1'] = datatd.header['CDELT1']
    hdr['CDELT2'] = datatd.header['DAXIS3']
    hdr['OBS_L0'] = (datatd.header['OBS_L0'], 'degree')
    hdr['OBS_B0'] = (datatd.header['OBS_B0'], 'degree')
    hdr['PANGLE'] = (datatd.header['PANGLE'], 'degree')
    hdr['REF_L0'] = (datatd.header['REF_L0'], 'Carrington Longitude (degree)')
    hdr['REF_B0'] = (datatd.header['REF_B0'], 'Carrington Latitude (degree)')
    hdr['DSUN_OBS'] = (datatd.header['DSUN_OBS'], 'm')
    hdr['X0'] = (x0, 'pix')
    hdr['Y0'] = (y0, 'pix')
    hdr['THETA1'] = (th0, 'degree')
    hdr['THETA2'] = (th1, 'degree')
    hdr['SQ_INT'] = (round(float(sq_int), 3), 'Normalized ripple '
                     'intensity (0: no ripple)')
    hdr['CUNIT1'] = ('arcsec  ', '[arcsec] CUNIT1: arcsec')
    hdr['CUNIT2'] = ('secs    ', 'same as T_REC_unit')
    hdr['TELESCOP'] = ('SDO/HMI ', 'Telescope')
    hdr['T_REC'] = datatd.header['T_REC']
    # save into file
    hdu.writeto(name, overwrite=True)
    exit()


def ray_path(xoff=0, toff=0):
    import subprocess
    import numpy as np

    fig = plt.figure(figsize=(9, 5))
    plt.subplots_adjust(wspace=0.1, left=0.08, right=0.85, top=0.8)

    plt.subplot(121)
    image, im = datatd.plot(colorbar=False, vmin=vmin, vmax=vmax,
                            interpolation='sinc')
    ext = im.get_extent()
    frame0_min = frame0*daxis3/60
    im.set_extent([ext[0], ext[1], ext[2]-frame0_min, ext[3]-frame0_min])
    plt.grid(ls='--', color='black', alpha=0.3)

    ax2 = plt.subplot(122)
    image, im = datatd.plot(colorbar=False, vmin=vmin, vmax=vmax,
                            interpolation='sinc')
    im.set_extent([ext[0], ext[1], ext[2]-frame0_min, ext[3]-frame0_min])
    subprocess.run('cp -rf $CHARLIE/tables/SOUND_SPEEDS .', shell=True,
                   check=True)
    subprocess.run('seq 0 0.001 0.2 > wkb_dist.txt', shell=True, check=True)
    subprocess.run(f'wkb SOUND_SPEEDS {depth} 0 1 0.0001 < wkb_dist.txt ' +
                   '> ray_path.txt', shell=True, check=True)
    dist, time, ampl, cos_a, cos_b = np.loadtxt('ray_path.txt', unpack=True)
    time /= 60  # Time min
    dist *= solar_radius
    plt.plot(dist+xoff, time+toff, lw=2, zorder=10)
    extent = im.get_extent()
    im.set_extent(extent)
    plt.ylabel(' ')
    plt.title(' ')
    plt.grid(ls='--', color='black', alpha=0.3)
    ax2.axes.get_yaxis().set_ticklabels([])
    cbar_ax = fig.add_axes([0.87, 0.14, 0.024, 0.64])

    fig.colorbar(im, cax=cbar_ax, label=r'$\Delta$ v$_{LOS}$ (m/s)',
                 orientation='vertical', location='right')

    # image_cart = data.toCartesian
    plt.show()
    exit()


ray_path(xoff=0, toff=0)

save('TD_diagram.fits', sq_int=0.40)
