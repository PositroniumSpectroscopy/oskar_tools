#!python
""" dset: assign default settings using the command line

    Copyright (c) 2014-2016, UNIVERSITY COLLEGE LONDON
    @author: Adam Deller
"""
from __future__ import print_function, division
import argparse
import oskar

def dset(settings, args):
    """ parse args to assign dire, rid or arb defaults. """
    if args.verbose:
        print(settings.fil)
    # data dire
    if args.dire is not None:
        # overwrite default dire
        if args.verbose:
            print('assign \'dire\' =', args.dire)
        settings.assign('dire', args.dire)
    # rid
    if args.rid is not None:
        # overwrite default rid
        if args.verbose:
            print('assign \'rid\' =', args.rid)
        settings.assign('rid', args.rid)
    # arb
    if args.arb is not None:
        if args.verbose:
            print('assign', args.arb[0], '=', [args.arb[1]])
        settings.assign(args.arb[0], [args.arb[1]])

    # append
    if args.append is not None:
        if args.verbose:
            print('append', args.arb[0], 'with', [args.arb[1]])
        val = settings.values[args.append[0]]
        val.append(args.append[1])
        settings.assign(args.append[0], val)

    # delete
    if args.drop is not None:
        if args.verbose:
            print('drop', args.drop[0])
        settings.drop(args.drop[0])

    # save
    if not args.dry_run:
        if args.verbose:
            print('save defaults.json')
        settings.save()

    # print
    if args.pprint:
        if args.verbose:
            print('pprint defaults')
        settings.pprint()

def main():
    """ assign values to default.json using the command line."""
    # inputs
    # data options
    parser = argparse.ArgumentParser(description='Set default.json values.')
    dat_parser = parser.add_argument_group('HDF5 data')
    dat_parser.add_argument('-d', '--dire', nargs=1, default=None,
                            help='data directory, e.g. --dire \"Z:\\Data\"')
    dat_parser.add_argument('-r', '--rid', nargs=1, default=None,
                            help='rid, e.g. --rid \"20160203_185633\"')
    dat_parser.add_argument('-a', '--arb', nargs=2, default=None,
                            help='arb, e.g. --arb \"average\" \"{\'SSPALS\': \
                            [\'t0\',\'DF\',\'Range\',\'FWHM\']}\"')
    dat_parser.add_argument('-x', '--append', nargs=2, default=None,
                            help='append value to defaults, name value')
    dat_parser.add_argument('-u', '--drop', nargs=1, default=None,
                            help='drop from defaults, e.g. --drop \"average\"')

    # script options
    script_parser = parser.add_argument_group('script options')
    script_parser.add_argument('-n', '--dry_run', action="store_true", default=False,
                               help="do not override defaults.json")
    script_parser.add_argument('-p', '--pprint', action="store_true", default=False,
                               help="print defaults")
    script_parser.add_argument('-v', '--verbose', action="store_true", default=False,
                               help="verbose script information")

    args = parser.parse_args()

    # defaults
    settings = oskar.Defaults()
    # set defaults
    dset(settings, args)

if __name__ == "__main__":
    # execute only if run as a script
    main()
