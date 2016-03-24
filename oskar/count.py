#!python
""" count: find triggers per trace. Output DataFrame to pandas .pkl file(s). ftypes
    are specfied in defaults.json (key: counts), or by using option -f.

    Copyright (c) 2015-2016, UNIVERSITY COLLEGE LONDON
    @author: Adam Deller
"""
from __future__ import print_function, division
import os
import sys
import argparse
import json
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import oskar

def find_triggers_1D(ydata, dt, min_x, min_y):
    """ Find the position, amplitude and duration of all trigger events in a
        1D array (ydata).  Trigger events are the leading edge of pulses in the
        array that surpase min_y for a period in excess of min_x.

        Returns each event in a mixed type array.
    """
    # initialise array
    all_trigs = np.array([], dtype=([('trigger', 'float64'),
                                     ('width', 'float64'),
                                     ('amp', 'float64')]))
    # find where trace exceeds threshold
    threshold = ydata > min_y
    # first and last element must be False
    threshold[0] = False
    threshold[-1] = False
    # find where True goes to False (pulse start/ stop)
    rngs = np.where(np.diff(threshold))[0].reshape(-1, 2)
    # find the pulses wider than the min_width
    pulse = rngs[rngs[:, 1]-rngs[:, 0] > min_x]
    # leading edge, width, and max value
    trigger = pulse[:, 0]
    width = pulse[:, 1] - pulse[:, 0]
    amp = [np.max(ydata[p[0]:p[1]]) for p in pulse]
    trig_arr = np.transpose([trigger, width, amp])
    # output
    for trig in trig_arr:
        newline = np.array([(trig[0]*dt, trig[1]*dt, trig[2])], dtype=all_trigs.dtype)
        all_trigs = np.append(all_trigs, newline)
    return all_trigs

def main():
    """ find triggers per trace. Output DataFrame to pandas .pkl file(s). """
    # defaults
    settings = oskar.Defaults()
    ddir = settings.values['dire'] if 'dire' in settings.values else None
    drid = settings.values['rid'] if 'rid' in settings.values else None
    dftype = settings.values['count'] if 'count' in settings.values else None
    # inputs
    # data options
    parser = argparse.ArgumentParser(description='Count the triggers per trace \
                                    in FTYPE. Output DataFrame to .pkl file.')
    dat_parser = parser.add_argument_group('HDF5 data')
    dat_parser.add_argument('-d', '--dire', nargs=1, default=ddir,
                            help='data directory, e.g. --dire \"Z:\\Data\"')
    dat_parser.add_argument('-r', '--rid', nargs=1, default=drid,
                            help='rid, e.g. --rid \"20160203_185633\"')
    dat_parser.add_argument('-f', '--ftype', nargs='+', default=dftype,
                            help="file type(s) to read, e.g. -f \"CH_A0\"")
    # trigger settings
    trigger_parser = parser.add_argument_group('trigger')
    trigger_parser.add_argument('-n', '--negative', action="store_true", default=False,
                                help="negative trigger [default: False]")
    trigger_parser.add_argument('--n_bsub', dest='n_bsub', type=int, default=50,
                               help="number of data points to use for background subtraction,\
                               default 50.")
    trigger_parser.add_argument('-l', '--min_level', nargs=1,
                                default=0.001, type=float,
                                help='trigger level (V) [default: 0.001]')
    trigger_parser.add_argument('-w', '--min_width', nargs=1,
                                default=1E-7, type=float,
                                help='minimum trigger pulse width (s) [default: 1E-7]')
    # script options
    script_parser = parser.add_argument_group('script options')
    script_parser.add_argument('-i', '--info', action="store_true", default=False,
                               help="pprint h5 file attributes")
    script_parser.add_argument('-b', '--progress', action="store_true", default=False,
                               help="display progress bar")
    script_parser.add_argument('-q', '--quiet', action="store_true", default=False,
                               help="surpress script output")
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

    # ftype
    if args.ftype is not None:
        settings.assign('count', args.ftype)
        ftypes = args.ftype
    else:
        raise Exception("ftype not specified nor found in defaults.json \
                        (key: count). Use flag --ftype")
    if args.set:
        settings.save()
    min_width = args.min_width
    min_level = args.min_level
    n_bsub = args.n_bsub
    # run counting analysis
    for rid in args.rid:
        if not (args.info or args.quiet):
            print(rid)
        if args.verbose:
            print('loading hdf5 file')
        h5 = oskar.H5Data(rid, args.dire[0])
        if args.info:
            h5.pprint()
        # output
        if args.verbose:
            print('building output path')
        out_dire = h5.out_dire('Count')
        if args.verbose:
            print(out_dire)
        # load hdf5 file
        with h5py.File(h5.fil, 'r') as dfil:
            data = dfil['.']
            # get squid values
            squids = np.sort(np.array([int(key) for key in data.keys()]))
            # count for each ftype
            for name in ftypes:
                sys.stdout.flush()
                if args.verbose:
                    print("counting trigger events")
                rtot = 0 # number of traces
                nsqs = 0 # number of squids
                all_events = []
                # loop over each squid
                for sq in tqdm(squids, smoothing=0.1, unit=' squids',
                               disable=bool(not args.progress), leave=True):
                    squid = str(sq)
                    if name in data[squid]:
                        nsqs = nsqs + 1
                        # load traces
                        trace = np.array(data[squid][name])
                        if np.shape(trace)[0] > 0:
                            rtot = rtot + np.shape(trace)[0]
                            # read oscilloscope attributes
                            osc = oskar.arr2dict(data[squid][name].attrs.items())
                            dt = osc['dt']
                            min_x = min_width / dt
                            min_y = min_level
                            if args.negative:
                                # invert if negative
                                trace = np.negative(trace)
                            for rep, ydata in enumerate(trace):
                                if n_bsub != 0:
                                    ydata = ydata - np.mean(ydata[:n_bsub])
                                events = pd.DataFrame(find_triggers_1D(ydata, dt, min_x, min_y))
                                events['squid'] = sq
                                events['rep'] = rep + 1
                                all_events.append(events)
                if len(all_events) == 0:
                    if not args.quiet:
                        print(name, ": no triggers found")
                else:
                    allDF = pd.concat(all_events, axis=0)
                    del all_events
                    allDF = allDF.reset_index(drop=True)
                    allDF.index.name = 'EVENT'
                    if len(allDF.index) > 0:
                        #save averaged data
                        out_fil = os.path.join(out_dire, name + "_triggers.pkl")
                        allDF.to_pickle(out_fil)
                        tot = len(allDF.index)
                        if not args.quiet:
                            print(name, ": Found", str(tot), " trigger events from", rtot,
                                  "traces, for", str(nsqs), "out of", str(len(squids)), "seqs.")
                    elif not args.quiet:
                        print(name, ": no triggers found")
            if args.verbose:
                print("saving trigger settings")
            # store settings as a dict
            vals = dict()
            vals['min_level'] = min_level
            vals['min_width'] = min_width
            fstw = json.dumps(vals, sort_keys=False, indent=5, separators=(',', ': '))
            # save settings to json file
            out_fil = os.path.join(out_dire, "settings.json")
            with open(out_fil, 'w') as out:
                out.write(fstw)
        if args.verbose:
            print("Done!")

if __name__ == "__main__":
    # execute only if run as a script
    main()
