#! /usr/bin/env python

"""Create a JSON file with image metadata, written to data.json.

Authors
-------

    Matthew Bourque

Use
---

    This module can be executed from the command line as such:

        >>> python make_json.py


Dependencies
------------

    - python 3
    - astropy
    - numpy
    - PIL
    - pyql
"""

import glob
import json
import os

from astropy.io import fits
import numpy as np
from PIL import Image

from acsql.database.database_interface import session as acs_session
from acsql.database.database_interface import Master as acs_Master
from acsql.database.database_interface import WFC_flt_0
from acsql.utils.utils import SETTINGS
from pyql.database.ql_database_interface import session as wfc3_session
from pyql.database.ql_database_interface import Master as wfc3_Master
from pyql.database.ql_database_interface import UVIS_flt_0
from pyql.database.ql_database_interface import IR_flt_0


def get_images_to_process(instrument):
    """Build a colletion of FITS files to create JSON data from.

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

        uvis_results = wfc3_session.query(wfc3_Master.ql_root, wfc3_Master.dir)\
            .join(UVIS_flt_0)\
            .filter(UVIS_flt_0.pr_inv_l == 'Levay').all()

        ir_results = wfc3_session.query(wfc3_Master.ql_root, wfc3_Master.dir)\
            .join(IR_flt_0)\
            .filter(IR_flt_0.pr_inv_l == 'Levay').all()

        results = uvis_results
        results.extend(ir_results)

        images = ['{}/{}q_flt.fits'.format(item[1], item[0]) for item in results]

        return images

    elif instrument == 'ACS':

        results = acs_session.query(acs_Master.path, WFC_flt_0.filename)\
            .join(WFC_flt_0)\
            .filter(WFC_flt_0.pr_inv_f == 'Zolt')\
            .filter(WFC_flt_0.pr_inv_l == 'Levay').all()

        images = ['{}{}/{}'.format(SETTINGS['filesystem'][:-1], item[0], item[1]) for item in results]

        return images

    elif instrument == 'STIS':

        images = glob.glob('/stis/data/*x2d.fits')
        return images


def make_file_dict(filename, instrument):
    """Gather image metadata and into a dictionary.

    Parameters
    ----------
    filename : str
        The path to the file to process.
    instrument : str
        Can be 'wfc3', 'acs', or 'stis'

    Returns
    -------
    file_dict : dict
        A dictionary containing the information to write to the JSON
        file.
    """

    print('Creating JSON for {}'.format(filename))

    # Build JSON for the filename
    file_dict = {}
    file_dict['instrument'] = instrument
    file_dict['rootname'] = os.path.basename(filename).split('_')[0]
    file_dict['jpeg_url'] = 'http://localhost:8080/data/{}/{}/jpg/image'.format(file_dict['instrument'], file_dict['rootname'])
    file_dict['thumb_url'] = 'http://localhost:8080/data/{}/{}/thumb/image'.format(file_dict['instrument'], file_dict['rootname'])

    # Add the metadata from the headers
    try:
        metadata = fits.getheader(filename, 0)
        exclude_list = ['HISTORY', 'COMMENT', ''] # keys to ignore
        file_dict['metadata'] = {}
        for key, value in metadata.items():
            key = key.strip()
            if key in exclude_list or value == '':
                continue
            file_dict['metadata'][key.lower()] = value
    except FileNotFoundError:
        file_dict['metadata'] = 'No metadata available.'

    return file_dict


if __name__ == '__main__':

    wfc3_images = get_images_to_process('WFC3')
    acs_images = get_images_to_process('ACS')
    stis_images = get_images_to_process('STIS')

    # Initialize the JSON file
    json_dict = {}
    json_dict['rootnames'] = []

    for filename in wfc3_images:
        file_dict = make_file_dict(filename, 'wfc3')
        json_dict['rootnames'].append(file_dict)

    for filename in acs_images:
        file_dict = make_file_dict(filename, 'acs')
        json_dict['rootnames'].append(file_dict)

    for filename in stis_images:
        file_dict = make_file_dict(filename, 'stis')
        json_dict['rootnames'].append(file_dict)

    # Write the JSON file
    with open('models/data.json', 'w') as outfile:
        json.dump(json_dict, outfile, indent=4)
