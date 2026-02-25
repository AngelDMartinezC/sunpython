#!/usr/bin/env python
'''
Function to cut a sunpy map given the Carrington Longitude and Latitude

This works when there is a Postel projected map and the right WCS has been
added to the header. This is done by applying xlate_head2wcs_new to a
single unstacked frame.
'''

import astropy.units as u
from astropy.coordinates import SkyCoord


def map_cut(smap, x0, y0, width, height, coord=False):
    '''
    Returns a cropped map of size size `width x height` centered around
    `(x0, y0)`.

    Parameters
    ----------
    smap : sunpy.map.Map
        Input SunPy map
    x0, y0 : float or astropy.units.Quantity
        center of the field of view
    width : float or astropy.units.Quantity
        width of cropped map
    height : float or astropy.units.Quantity
        height of cropped map

    Returns
    -------
    sunpy.map.Map
        Cropped map
    '''
    unit_x = smap.spatial_units.axis1
    unit_y = smap.spatial_units.axis2
    x0 = u.Quantity(x0, unit_x)
    y0 = u.Quantity(y0, unit_y)
    x_size = u.Quantity(width, unit_x)
    y_size = u.Quantity(height, unit_y)
    c0 = x0 - x_size/2
    c1 = y0 - y_size/2
    c2 = x0 + x_size/2
    c3 = y0 + y_size/2
    bottom_left = SkyCoord(c0, c1, frame=smap.coordinate_frame)
    top_right = SkyCoord(c2, c3, frame=smap.coordinate_frame)
    submap = smap.submap(bottom_left, top_right=top_right)
    if coord:
        return submap, bottom_left, top_right
    else:
        return submap


def map_cut_box(smap, x0, x1, y0, y1, coord=False):
    '''
    Returns a cropped map inside a rectangular region befined by bottom-left
    and top-right coordinates.

    Parameters
    ----------
    smap : sunpy.map.Map
       Input SunPy map
    x0, y0 : float or astropy.units.Quantity
       Bottom-left coordinates of the rectangular region.
    x1, y1 : float or astropy.units.Quantity
        Top-right coordinates of the rectangular region.

    Returns
    -------
    sunpy.map.Map
        Cropped map
    '''
    if (x0.unit and y0.unit and x1.unit and y1.unit) is u.pixel:
        bottom_left = smap.pixel_to_world(x0, y0)
        top_right = smap.pixel_to_world(x1, y1)
        submap = smap.submap(bottom_left, top_right=top_right)
    else:
        unit_x = smap.spatial_units.axis1
        unit_y = smap.spatial_units.axis2
        x0 = u.Quantity(x0, unit_x)
        y0 = u.Quantity(y0, unit_y)
        x1 = u.Quantity(x1, unit_x)
        y1 = u.Quantity(y1, unit_y)
        top_right = SkyCoord(x1, y1, frame=smap.coordinate_frame)
        bottom_left = SkyCoord(x0, y0, frame=smap.coordinate_frame)
        submap = smap.submap(bottom_left, top_right=top_right)
    if coord:
        return submap, bottom_left, top_right
    else:
        return submap


if __name__ == "__main__":

    # Read map
    import matplotlib.pyplot as plt
    from sunpy.map import Map
    import sys

    map_input = Map(sys.argv[1])

    # Carrington coordinates
    crln = 178.5432
    crlt = -16.70895
    width = 500
    height = 500

    submap = map_cut(map_input, crln, crlt, width, height)

    plt.figure(figsize=(10, 4))
    ax1 = plt.subplot(121, projection=map_input)
    map_input.plot()
    plt.title('Original map')
    plt.colorbar()

    ax2 = plt.subplot(122, projection=submap)
    submap.plot()
    plt.colorbar()
    plt.title('cut map')
    plt.show()
