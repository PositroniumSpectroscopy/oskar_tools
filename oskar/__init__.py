#! python
from __future__ import print_function, division  # Require Python 2.6 or later
import os
from .oskar import *

colors = ['#4584CC', '#479657', '#C22F27', '#9354B0', '#ADB056', '#5FC7C2', '#ED8721']
markers = ['o', 's', '^', 'D', '*', 'h', 'x']

class Scripts(object):
    """ paths to scripts """
    def __init__(self):
        script_path = os.path.join(oskar.MOD_PATH, 'oskar')
        self.dset = os.path.join(script_path, 'dset.py')               # set defaults
        self.info = os.path.join(script_path, 'info.py')               # get run info
        self.average = os.path.join(script_path, 'average.py')         # average data
        self.vrange = os.path.join(script_path, 'vrange.py')           # find vertical ranges
        self.sspals = os.path.join(script_path, 'sspals_.py')          # reanalyse sspals
        self.count = os.path.join(script_path, 'count.py')             # count events
