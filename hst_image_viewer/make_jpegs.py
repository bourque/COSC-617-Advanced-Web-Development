#! /usr/bin/env python

"""Create a JPEG and thumbnail file from a RAW, FLT, or FLC FITS file
from the ACS or WFC3 instrument.  Currently, other instruments (e.g.
STIS) and other filetyles (e.g. IMA) are not not supported, but likely
the code could be easily extended to support these.

Authors
-------

    Matthew Bourque
    Alex Viana
    Meredith Durbin

Use
---

    This module is intended to be used via the python environemnt, for
    example:

        import make_jpeg
        make_jpeg(<path_to_file>, <outfile>)
        make_thumbnail(<outfile>)

    Where <path_to_file> is the path to the FITS file to create a JPEG
    from and <outfile> is the path to where the JPEG (and subsequently
    the thumbnail) will be written.

Dependencies
------------

    This module is supported in both Python 2.7 and Python 3+.  This
    module also depends on astropy, numpy, and the Python Image Library
    (PIL), also known as 'Pillow'.  These can be installed via 'conda'
    or 'pip':

    conda/pip install astropy
    conda/pip install numpy
    conda/pip install Pillow
"""

import os

from astropy.io import fits
import numpy as np
from PIL import Image

from pyql.database.ql_database_interface import session, Master, UVIS_flt_0, IR_flt_0


def get_images_to_process(instrument):
    """Build a colletion of FITS files to create JPEGs from.

    Parameters
    ----------
    instrument : str
        The instrument. Can be 'WFC3', 'ACS', or 'STIS'.

    Returns
    -------
    images : list
        A list of paths to FITS files to create JPEGs for.
    """

    if instrument == 'WFC3':

        uvis_results = session.query(Master.ql_root, Master.dir)\
            .join(UVIS_flt_0)\
            .filter(UVIS_flt_0.pr_inv_l == 'Levay').all()

        ir_results = session.query(Master.ql_root, Master.dir)\
            .join(IR_flt_0)\
            .filter(IR_flt_0.pr_inv_l == 'Levay').all()

        results = uvis_results
        results.extend(ir_results)

        images = ['{}/{}q_flt.fits'.format(item[1], item[0]) for item in results]

        return images

    elif instrument == 'ACS':
        return []

    elif instrument == 'STIS':
        return []


def make_jpeg(filename, outfile):
    """Create a JPEG file from a raw, flt, or flc FITS file.  If the
    image is a full-frame WFC/UVIS or ACS/WFC image, the data from the
    1st and 4th extensions will be combined into one JPEG.

    Parameters
    ----------
    filename : str
        The path to the file.
    outfile : str
        The path to where the JPEG will be saved.
    """

    print('Creating JPEG for {}'.format(filename))

    try:
        with fits.open(filename) as hdulist:

            # Get the image data.
            data = hdulist[1].data

            # For full UVIS or WFC image.
            if (hdulist[0].header['DETECTOR'] in ['UVIS', 'WFC']
            and len(hdulist) > 4
            and hdulist[4].header['EXTNAME'] == 'SCI'):
                data2 = hdulist[4].data
                height = data.shape[0] + data2.shape[0]
                width = data.shape[1]
                temp = np.zeros((height, width))
                temp[0:int(height/2), :] = data
                temp[int(height/2):height, :] = data2
                data = temp

        # Clip the top and bottom 1% of pixels.
        top = np.percentile(data, 99)
        data[data > top] = top
        bottom = np.percentile(data, 1)
        data[data < bottom] = bottom

        # Scale the data.
        data = data - data.min()
        data = (data / data.max()) * 255.
        data = np.flipud(data)
        data = np.uint8(data)

        # Write the image to a JPEG
        image = Image.fromarray(data)
        image.save(outfile)

    except FileNotFoundError:
        pass


def make_thumbnail(filename):
    """Create a thumbnail file from a JPEG file.

    Parameters
    ----------
    filename : str
        The path to the JPEG file.
    outfile : str
        The path to where the thumbnail will be saved.
    """

    print('Creating thumbnail for {}'.format(filename))

    try:
        image = Image.open(filename)
        image.thumbnail((256, 256), Image.ANTIALIAS)
        image.save(filename.replace('.jpg', '.thumb'), 'JPEG')
    except FileNotFoundError:
        pass


if __name__ == '__main__':

    wfc3_images = get_images_to_process('WFC3')
    for filename in wfc3_images:
        outfile = os.path.join('/user/bourque/repositories/COSC-617-Advanced-Web-Development/hst_image_viewer/models/wfc3/', os.path.basename(filename).replace('_flt.fits', '.jpg'))
        make_jpeg(filename, outfile)
        make_thumbnail(outfile)
