#!python
""" info: read hdf5 attributes from a given run

    Copyright (c) 2014-2016, UNIVERSITY COLLEGE LONDON
    @author: Adam Deller
"""
from __future__ import print_function, division
import argparse
import pandas as pd
import oskar

def main():
    """ read hdf5 attributes and print result. """
    # defaults
    settings = oskar.Defaults()
    ddir = settings.values['dire'] if 'dire' in settings.values else None
    drid = settings.values['rid'] if 'rid' in settings.values else None

    # inputs
    # data options
    parser = argparse.ArgumentParser(description='Get HDF5 attributes for\
                                a given run.')
    dat_parser = parser.add_argument_group('HDF5 data')
    dat_parser.add_argument('-d', '--dire', nargs=1, default=ddir,
                            help='data directory, e.g. --dire \"Z:\\Data\"')
    dat_parser.add_argument('-r', '--rid', nargs=1, default=drid,
                            help='rid, e.g. --rid \"20160203_185633\"')

    # script options
    script_parser = parser.add_argument_group('script options')
    script_parser.add_argument('-v', '--verbose', action="store_true", default=False,
                               help="verbose script information")
    # defaults.json
    def_parser = parser.add_argument_group('defaults.json')
    def_parser.add_argument('-s', '--set', action="store_true", default=False,
                            help="save args. (dire/ rid/ ftype) as default values")
    args = parser.parse_args()
    # data dire
    if args.dire is not None:
        # overwrite default dire
        settings.assign('dire', args.dire)
    else:
        raise Exception("data directory not specified, nor found in defaults.json. Use flag --dire")
    # rid
    if args.rid is not None:
        # overwrite default rid
        settings.assign('rid', args.rid)
    else:
        raise Exception("rid not specified, nor found in defaults.json. Use flag --rid")

    # get info
    for rid in args.rid:
        if args.verbose:
            print('Loading hdf5 file.')
        h5 = oskar.H5Data(rid, args.dire[0])
        h5.pprint()
        if args.verbose:
            print('Loading squid, var, and rec info.')
        h5.load_log(update=True)
        vDF = pd.DataFrame(h5.vDF())
        print("SQUIDS:", str(len(vDF.index)))
        print("VARS:", str(vDF.columns.values))
        rDF = pd.DataFrame(h5.rDF())
        print("REC:", str(rDF.columns.values))
        if args.verbose:
            print('Searching for unique var permutations.')
        uDF = pd.DataFrame(h5.uDF())
        print("Unique VAR combinations:", str(len(uDF.index)))
        print("Number of loops:", len(vDF.index)/ len(uDF.index))
        if args.verbose:
            print("Done!")

if __name__ == "__main__":
    # execute only if run as a script
    main()
