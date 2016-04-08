#! python
""" oskar_tools: python tools for working with hdf5 structured datasets obtained using oskar.

    Copyright (c) 2014-2016, UNIVERSITY COLLEGE LONDON
    @author: Adam Deller
"""
from __future__ import print_function, division
import os
import ast
import json
import re
import glob
import h5py
import numpy as np
import pandas as pd

MOD_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

#    -------------
#    defaults.json
#    -------------

class Defaults(object):
    """ default settings for running analysis """
    def __init__(self, path=MOD_PATH):
        self.fil = os.path.abspath(os.path.join(path, 'defaults.json'))
        self.load()

    def load(self):
        """ load default.json from pyosq module """
        if os.path.isfile(self.fil):
            with open(self.fil, "r") as dfil:
                fstr = dfil.read()
            self.values = ast.literal_eval("".join(fstr.split()))
        else:
            print('Warning: default.json file not found.')
            self.values = dict()

    def save(self):
        """ save changes to default.json """
        fstw = json.dumps(self.values, sort_keys=False, indent=5,
                          separators=(',', ': '))
        with open(self.fil, 'w') as dfil:
            dfil.write(fstw)

    def drop(self, name):
        """ delete attribute [name] """
        if name in self.values:
            del self.values[name]
        else:
            pass

    def assign(self, name, value):
        """ assign attribute [name]: [value] """
        self.values[name] = value

    def rid(self, value):
        """ assign attribute 'rid': [value] """
        self.values['rid'] = value

    def dire(self, value):
        """ assign attribute 'dire': [value] """
        self.values['dire'] = value

    def pprint(self):
        """ print default settings """
        fstw = json.dumps(self.values, sort_keys=False, indent=5,
                          separators=(',', ': '))
        print('Default settings:\n', fstw)

#    -------------
#    h5 data files
#    -------------

class H5Data(object):
    """ information relating to hdf5 data """
    def __init__(self, rid=None, base=None):
        settings = Defaults()
        if rid is None:
            # get defaults
            if 'rid' in settings.values:
                self.rid = settings.values['rid'][0]
            else:
                raise Exception("rid not specified, nor found in defaults.json.")
        else:
            self.rid = rid
        if base is None:
            # get defaults
            if 'dire' in settings.values:
                self.base = settings.values['dire'][0]
            else:
                raise Exception("base directory not specified, nor found in defaults.json.")
        else:
            self.base = base
        self.build()
        if not os.path.isdir(self.dire):
            # no dire found for run id
            raise IOError(self.dire + " not found")
        if not os.path.isfile(self.fil):
            # no h5 file found for rid
            raise IOError(self.fil + " not found")
        self.log = None

    def build(self):
        """ build directory path to data file """
        year = self.rid[:4]
        month = self.rid[4:6]
        day = self.rid[6:8]
        self.dire = os.path.join(self.base, year, month, day, self.rid)
        self.fil = os.path.join(self.dire, str(self.rid + "_raw.h5"))
        self.log_file = os.path.join(self.dire, str(self.rid + "_log.pkl"))

    def info(self):
        """ read h5 attributes """
        with h5py.File(self.fil, 'r') as dfil:
            attrs = arr2dict(dfil.attrs.items())
        return attrs

    def squid_info(self, squid):
        """ read h5 attributes for squid entry """
        with h5py.File(self.fil, 'r') as dfil:
            data = dfil['.']
            if str(squid) in data:
                attrs = arr2dict(data[str(squid)].attrs.items())
                return attrs
            else:
                raise Exception("self.log is None.  Try method self.load_log().")

    def pprint(self):
        """ print author and description info """
        print(self.rid)
        inf = self.info()
        print(u'  ', u'Author:\t', inf['Author'])
        print(u'  ', u'Description: ', inf['Description'].replace('\n', '\n\t\t '))

    # H5 attributes
    def load_log(self, update=False):
        """ read from log_.pkl if exists unless update=True,
            otherwise read squid attributes, save as [rid]_log.pkl.
        """
        if self.log is None:
            if os.path.isfile(self.log_file) and not update:
                # if exists load log file
                self.log = pd.read_pickle(self.log_file)
            else:
                with h5py.File(self.fil, 'r') as dfil:
                    # read attributes from each squid
                    data = dfil['.']
                    squids = np.sort(np.array([int(key) for key in data.keys()]))
                    all_vars = []
                    for sq in squids:
                        # read info
                        info = arr2dict(data[str(sq)].attrs.items())
                        all_vars.append(pd.DataFrame([info], index=[sq]))
                    logDF = pd.concat(all_vars)
                    logDF.index.name = 'SQUID'
                    #remove duplicate squid column
                    logDF.drop('SQUID', axis=1, inplace=True)
                    if 'DATETIME' in logDF:
                        logDF.DATETIME = pd.to_datetime(logDF.DATETIME)
                    if 'ELAPSED' in logDF:
                        # legacy fix: rename Elapsed to Acquire
                        logDF.rename(columns={'ELAPSED': 'ACQUIRE'}, inplace=True)
                    self.log = logDF
                # save to pickle file
                logDF.to_pickle(self.log_file)

    def vDF(self):
        """ return just the VAR values from the log file """
        if self.log is not None:
            DF = self.log
            DF = DF.filter(regex="VAR:")
            DF.columns = [re.split('^VAR:', x)[1] for x in DF.columns.values]
            return DF
        else:
            raise Exception("self.log is None.  Try method self.load_log().")

    def rDF(self):
        """ return just the REC values from the log file """
        if self.log is not None:
            DF = self.log
            DF = DF.filter(regex="REC:")
            DF.columns = [re.split('^REC:', x)[1] for x in DF.columns.values]
            return DF
        else:
            raise Exception("self.log is None.  Try method self.load_log().")

    def uDF(self, lvars=None):
        """ return just the unique VAR values from the log file """
        if self.log is not None:
            DF = self.log
            DF = DF.filter(regex="VAR:")
            DF.columns = [re.split('^VAR:', x)[1] for x in DF.columns.values]
            return unqDF(DF, lvars=lvars)
        else:
            raise Exception("self.log is None.  Try method self.load_log().")

    # output
    def out_dire(self, folder="Analysis"):
        """ Build path to output directory.  Create if does not exist."""
        path = os.path.join(self.dire, folder)
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    ## array data (e.g., traces)
    def load_array(self, squids, name='CH_L0', ignore_missing=False):
        """ load HDF5 array data, default name=CH_L0.
            Returns 2D array of all arrays for squid(s) and hdf5 attributes as:
                numpy.array, dict
        """
        squids = np.array([squids]).flatten()
        trace = []
        with h5py.File(self.fil, 'r') as dfil:
            for sq in squids:
                squid = str(sq)
                if name in dfil[squid]:
                    trace.append(np.array(dfil[squid][name]))
                    info = dfil[squid][name].attrs.items()
                elif not ignore_missing:
                    raise Exception("Error: " + name + " not found for squid " \
                                    + squid + ".  Use ignore_missing=True if you don't care.")
        return np.concatenate(trace), arr2dict(info)

    ## DataFrame data (e.g., SSPALS)
    def load_DF(self, squids, name='SSPALS', cols=None, ignore_missing=False):
        """ load HDF5 DataFrame data, default name=SSPALS.
            Returns DataFrame of all arrays for squid(s) and hdf5 attributes as:
                pandas.DataFrame, dict
        """
        squids = np.array([squids]).flatten()
        all_data = []
        with h5py.File(self.fil, 'r') as dfil:
            data = dfil['.']
            for sq in squids:
                squid = str(sq)
                if str(sq) in data:
                    if name in data[squid]:
                        info = data[squid][name].attrs.items()
                        if cols is None:
                            tmp = pd.DataFrame(np.array(data[squid][name]))
                        else:
                            tmp = pd.DataFrame(np.array(data[squid][name]))[cols]
                        if len(tmp.index.values) > 0:
                            tmp['Repeat'] = tmp.index
                            tmp['SQUID'] = sq
                            all_data.append(tmp)
                    elif not ignore_missing:
                        raise Exception("Error: " + name + " not found for squid " \
                                    + squid + ".  Use ignore_missing=True if you don't care.")
                elif not ignore_missing:
                    raise Exception("Error: " + squid + " not found in data " \
                                + ".  Use ignore_missing=True if you don't care.")
        allDF = pd.concat(all_data, ignore_index=True)
        allDF = allDF.set_index(['SQUID', 'Repeat'])
        return allDF, arr2dict(info)

    ## load averaged data
    def load_av(self, **kwargs):
        """Load average data files.  Return pandas.DataFrame.
           default kwargs:
               loop=False     # data average over loops
               fils=[]        # file names to load from folder.  If empty loads *.dat
               ufil=False     # force use of unique_vars.dat, e.g., if analysed using lvars
               exclude=[]     # excluded file names.
        """
        # options
        loop = kwargs.get('loop', False)
        fils = kwargs.get('fils', [])
        ufil = kwargs.get('ufil', False)
        lvars = kwargs.get('lvars', None)
        exclude = kwargs.get('exclude', [])
        verbose = kwargs.get('verbose', False)
        #VAR files
        if loop:
            index = 'VID'
            ddir = 'Average_Loops'
            if ufil:
                unq_fil = os.path.join(self.dire, ddir, "unique_vars.dat")
                allDF = pd.read_csv(unq_fil, sep='\t')
            else:
                if self.log is not None:
                    allDF = self.uDF(lvars=lvars)
                else:
                    raise Exception("self.log is None.  Try method self.load_log().")
            exclude.append('unique_vars.dat')
        else:
            # all squids
            index = 'SQUID'
            ddir = 'Average'
            if self.log is not None:
                allDF = self.log
                allDF.columns = [re.split('^VAR:', x)[1] if len(re.split('^VAR:', x)) > 1 \
                                 else x for x in allDF.columns.values]
            else:
                raise Exception("self.log is None.  Try method self.load_log().")
        # data files
        if fils is None:
            return allDF
        else:
            if len(fils) == 0:
                fils = [x for x in glob.glob(os.path.join(self.dire, ddir, '*.dat')) \
                        if os.path.split(x)[1] not in exclude]
            if len(fils) == 0:
                raise IOError("no data files found")
            else:
                for fname in fils:
                    # read data
                    tmp = pd.read_csv(fname, sep="\t", index_col=index)
                    allDF = pd.merge(allDF, tmp, left_index=True, right_index=True, how='left')
                    if verbose:
                        print("Loaded: " + os.path.split(fname)[1])
                # convert columns to datetime
                if 'Date_Time' in allDF:
                    allDF['Date_Time'] = pd.to_datetime(allDF['Date_Time'])
                if 'DATETIME' in allDF:
                    allDF['DATETIME'] = pd.to_datetime(allDF['DATETIME'])
                return allDF

    ## load count data
    def load_count(self, names=['CH_A0'], folder='Count', **kwargs):
        """ Load count event data for all RIDS.  Return pandas DataFrame.
            Reads all files in run dires/ folder/ matching: name + '_triggers.pkl

            defaults:
                names = ['CH_A0']      # names of h5 data
                folder = Count         # folder where count info is stored
                include_vars = False   # merge with matching var information
        """
        include_vars = kwargs.get('include_vars', False)
        if include_vars:
            if self.log is None:
                raise Exception("self.log is None.  Try method self.load_log().")
            else:
                vvDF = self.vDF()
        all_dat = []
        for name in names:
            cfil = os.path.join(self.dire, folder, name + '_triggers.pkl')
            if os.path.exists(cfil):
                eventDF = pd.read_pickle(cfil)
                eventDF['FTYPE'] = name
                eventDF['EVENT'] = eventDF.index
                eventDF = eventDF.set_index(['FTYPE', 'EVENT'])
                all_dat.append(eventDF)
            else:
                raise IOError("pandas pickle file not found : " + name + '_triggers.pkl')
        allDF = pd.concat(all_dat)
        if include_vars:
            # combine VAR data and event data
            allDF = allDF.merge(vvDF, left_on='squid', right_index=True, how='left')
        return allDF

#    ----
#    data
#    ----
#    Functions for use with references to hdf5 files, e.g.
#
#        with h5py.File(h5.fil,'r') as dfil:
#            data = dfil['.']
#            arr = oskar.h5_array(data, 'CH_L1')

## attributes
def varDF(data):
    """ Read from an h5 object (data) the squid values and the VAR values
        associated with each. Return pandas DataFrame.
    """
    squids = np.sort(np.array([int(key) for key in data.keys()]))
    all_vars = []
    for sq in squids:
        # read info
        info = data[str(sq)].attrs.items()
        # filter VARs
        var = arr2dict([(re.split('^VAR:', x[0])[-1], x[1]) for x in info \
                        if re.match('^VAR:', x[0])])
        all_vars.append(pd.DataFrame([var], index=[sq]))
    vDF = pd.concat(all_vars)
    vDF.index.name = 'SQUID'
    return vDF

def recDF(data):
    """ Read from an h5 object (data) the squid values and the REC values
        associated with each. Return pandas DataFrame.
    """
    squids = np.sort(np.array([int(key) for key in data.keys()]))
    all_vars = []
    for sq in squids:
        # read info
        info = data[str(sq)].attrs.items()
        # filter REC values
        rec = arr2dict([(re.split('^REC:', x[0])[-1], x[1]) for x in info \
                        if re.match('^REC:', x[0])])
        all_vars.append(pd.DataFrame([rec], index=[sq]))
    rDF = pd.concat(all_vars)
    rDF.index.name = 'SQUID'
    return rDF

def unqDF(vDF, **kwargs):
    """ Read vDF (DataFrame) of VAR values, return uDF (DataFrame) containg
        unique combinations.  If kwarg 'lvars' is specified,  then use only
        those to find unique, otherwise use vDF.columns.
    """
    VARS = kwargs.get('lvars', None)
    if VARS is None:
        #if lvars not specified use ALL vars
        VARS = list(vDF.columns.values)
    elif not np.array([v in vDF.columns for v in VARS]).all():
        raise ValueError("at least one specified lvar was not found in the var list")
    # find unique
    if len(VARS) > 0:
        # sort and drop duplicate VAR combinations
        vDF = vDF.sort_values(VARS)
        uDF = vDF.drop_duplicates(VARS)[VARS]
        # re-index
        uDF['VID'] = np.arange(1, len(uDF.index) + 1, dtype='int32')
        uDF.set_index('VID', drop=True, inplace=True)
        return uDF
    else:
        # no VARS found
        raise Exception("No vars found for this rid.  Cannot average data over loops.")

## array data (e.g., traces)
def h5_array(data, squids, name='CH_L0', ignore_missing=False):
    """ Concatenate data[squid][name] traces for squid value(s) (integer dtype).
        ftype is expected to be a 2D array of repeat 1D measurements.
        Returns 2D array:
            np.array
    """
    squids = np.array([squids]).flatten()
    trace = []
    for sq in squids:
        squid = str(sq)
        if name in data[squid]:
            trace.append(np.array(data[squid][name]))
        elif not ignore_missing:
            raise Exception("Error: " + name + " not found for squid " \
                            + squid + ".  Use ignore_missing=True if you don't care.")
    return np.concatenate(trace)

## DataFrame data
def h5_DF(data, squids, name='SSPALS', cols=None, ignore_missing=False):
    """ Concatenate data[squid][name] values for squid value(s) (integer dtype).
        ftype is expected to be DataFrame-like (2D array with mixed type columns),
        e.g., SSPALS [t0, range, DF] (cols.) from repeat measurements (rows).
        If cols not spefied then return all.
    """
    squids = np.array([squids]).flatten()
    all_data = []
    for sq in squids:
        squid = str(sq)
        if squid in data:
            if name in data[squid]:
                if cols is None:
                    tmp = pd.DataFrame(np.array(data[squid][name]))
                else:
                    tmp = pd.DataFrame(np.array(data[squid][name]))[cols]
                if len(tmp.index.values) > 0:
                    tmp['Repeat'] = tmp.index
                    tmp['SQUID'] = sq
                    all_data.append(tmp)
            elif not ignore_missing:
                raise Exception("Error: " + name + " not found for squid " \
                            + squid + ".  Use ignore_missing=True if you don't care.")
        elif not ignore_missing:
            raise Exception("Error: " + squid + " not found in data " \
                        + ".  Use ignore_missing=True if you don't care.")
    allDF = pd.concat(all_data, ignore_index=True)
    allDF = allDF.set_index(['SQUID', 'Repeat'])
    return allDF

#    ----------------------
#    load multiple run data
#    ----------------------

def log_data(RIDS, **kwargs):
    """ Load and combine log files for multiple RIDS.

        defaults:
            update = False      # update log file (read hdf5 attributes -- slower)
    """
    update = kwargs.get('update', False)
    super_dat = []
    for rid in RIDS:
        h5 = H5Data(rid)
        h5.load_log(update=update)
        logDF = h5.log
        logDF['RID'] = rid
        logDF['SQUID'] = logDF.index
        super_dat.append(logDF)
    super_DF = pd.concat(super_dat)
    super_DF = super_DF.set_index(['RID', 'SQUID'])
    return super_DF

def average_data(RIDS, **kwargs):
    """Load average data files from all RIDS.  Return pandas DataFrame.

       default kwargs:
           update = False # update log file (read hdf5 attributes)
           loop=False     # data average over loops
           fils=[]        # file names to load from folder.  If empty loads *.dat
           ufil=False     # force use of unique_vars.dat, e.g., if analysed using lvars
           exclude=[]     # excluded file names.
    """
    # options
    loop = kwargs.get('loop', True)
    update = kwargs.get('update', False)
    ind = 'VID' if loop else 'SQUID'
    # get data
    super_dat = []
    for rid in RIDS:
        h5 = H5Data(rid)
        h5.load_log(update=update)
        tmp = h5.load_av(**kwargs)
        tmp['RID'] = rid
        tmp[ind] = tmp.index
        tmp = tmp.set_index(['RID', ind])
        super_dat.append(tmp)
    # concatenate
    allDF = pd.concat(super_dat)
    return allDF

def count_data(RIDS, names=['CH_A0'], folder='Count', include_vars=False, **kwargs):
    """ Load count event data for all RIDS.  Return pandas DataFrame.
        Reads all files in run dires/ folder/ matching:

            name + '_triggers.pkl

        defaults:
            NAMES = ['CH_A0']      # names of h5 data
            folder = Count         # folder where count info is stored
            include_vars = False   # merge with matching var information
            update = False         # update log file (read hdf5 attributes -- slower)
    """
    update = kwargs.get('update', False)
    super_dat = []
    for rid in RIDS:
        h5 = H5Data(rid)
        if include_vars:
            h5.load_log(update=update)
            vvDF = h5.vDF()
        all_dat = []
        for name in names:
            cfil = os.path.join(h5.dire, folder, name + '_triggers.pkl')
            if os.path.exists(cfil):
                eventDF = pd.read_pickle(cfil)
                eventDF['RID'] = rid
                eventDF['FTYPE'] = name
                eventDF['EVENT'] = eventDF.index
                eventDF = eventDF.set_index(['RID', 'FTYPE', 'EVENT'])
                all_dat.append(eventDF)
            else:
                raise IOError("pandas pickle file not found : " + name + '_triggers.pkl')
        allDF = pd.concat(all_dat)
        if include_vars:
            allDF = allDF.merge(vvDF, left_on='squid', right_index=True,
                                how='left')  # combine VAR data and event data
        super_dat.append(allDF)
    superDF = pd.concat(super_dat)
    return superDF

#    ---------------
#    misc. functions
#    ---------------

def out_dire(base, folder="Analysis"):
    """ Build path to output directory.  Create if does not exist."""
    path = os.path.join(base, folder)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def mean_stdev(DF, col):
    """ Calculate the mean standard deviation per loop using pandas"""
    return DF.groupby(['SQUID'])[col].std().mean()

def varcomp(row_a, row_b):
    """Compare two rows of values return True if all match, including NaNs."""
    test = ((row_a == row_b) | ((row_a != row_a) & (row_b != row_b))).all(1)
    return test

def arr2dict(arr):
    """Convert an array of [key, val] pairs into a dictionary."""
    dic = dict()
    for element in arr:
        dic[element[0]] = element[1]
    return dic
